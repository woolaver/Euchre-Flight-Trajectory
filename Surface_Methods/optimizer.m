
try
    swApp = actxGetRunningServer('SldWorks.Application');
catch
    swApp = actxserver('SldWorks.Application'); pause(5)
end
set(swApp, 'Visible', true);
spec = invoke(swApp, 'GetOpenDocSpec', 'N:\dotfiles\Desktop\euchre\blend_sharp.SLDPRT');
set(spec, 'Silent', true);
swModel = invoke(swApp, 'OpenDoc7', spec);




%%
%need to have flight regime info
M = 4.5;%Mach
alpha=deg2rad(15);
y = 1.15;%y as a function of mach

%these are hard coded because it will be very annoying to actually nfind
%these and they wont change that much
sref_lift = .42;
sref_drag = .03;



weight = 120*9.81;

rhoinf = 0.458; % Freestream density from standard atmosphere @30k ft
ainf =303; %speed of sound from standard atmosphere
qinf = 0.5*rhoinf*(M*ainf)^2; % 

%this gives the min required cl so it doesnt spaz out and give 0 cl and cd
min_cl = weight/qinf/sref_lift;


%%
%need to have a matrix with all the parameters I want to change
%              
%             'aft_heightFRAC','lowHT' 'bodyWD','bodyLG', 'humpLG %', 'humpWD','hmp_heightfrac', 'humpLowHT'];
parameters_0 = [.6              .074     .3      1.5      .5          .23        1.2           .09 ];

%%
%need to have lb and ub for the parameters

lb = [0   0 .15   1.25   .1    .1  .7  0 ];

ub = [.7 .15  .38   2   .7    .3   1.35  .15 ];

A = [0 0 -1 0 0 1 0 0];   % humpWD - bdyWD <= 0
B = 0;

%all of this for virginia

%             'bodyLG','noseLG' 'bodyWD','finWD', height fraction',];
%parameters_0 = [1.25     .74     .2        .075            .75 ];
% lb = [.5 .25 .15 .05 0 ];
% 
% ub = [1.5 1.5 .3 .125 1];
% 
% A = [1 1 0 0 0 ; %maxes out length
%      0 0 1 2 0; %maxes out width
% 
% 
% 
% ];
% 
% B= [2;.4
% 
% ];



%%
%need to have an optimizer surragote





options = optimoptions('surrogateopt', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 1000, ...
    'UseParallel', false, ...
    'PlotFcn', 'surrogateoptplot', ...
    'InitialPoints', parameters_0);     

fun = @(x) obj(x, swApp, swModel, alpha, M, sref_lift, sref_drag, y, min_cl);

% test_out = fun(parameters_0);
% disp(test_out.Ineq)
% return       

[optimal_shape, fval, exitflag, output, trials] = surrogateopt(fun, lb, ub, [], A,B, [], [], options);
%%






%need to have an objective function which takes inputs, creates cad file,
%runs the driver, and then spits out shit

function out = obj(x, swApp, swModel, alpha, minf, sref_lift, sref_drag, y, min_cl)
    out.Fval = 1e3;
    out.Ineq = [10;10];

    AX = 3;     % axial column in the CAD/STL frame 

    parameter_names = {'upHT','lowHT' 'bdyWD','bdyLG', 'humpLG', 'humpWD','humpUpHT', 'humpLowHT'};
    % parameter_names = {'bodyLG','noseLG','bodyWD','finWD','bodyHT'};
    R = .2;
    heightmax = R + sqrt(max(0, R^2 - (x(3)/2)^2));
    aft_height    = .15-x(2) + x(1)*(heightmax - (.15-x(2)));
    hmp_height = aft_height*x(7);



    parameter = [aft_height,x(2), x(3), x(4), x(5), x(6), hmp_height, x(8)];

    outputFile = 'N:\dotfiles\Desktop\euchre\blend_sharp.STL';
    if isfile(outputFile), delete(outputFile); end

    try
        eqMgr = invoke(swModel, 'GetEquationMgr');
        nEq   = double(invoke(eqMgr, 'GetCount'));
        nSet  = 0;
        for i = 0:(nEq-1)
            eqText = strtrim(get(eqMgr, 'Equation', i));
            for j = 1:numel(parameter_names)
                if startsWith(eqText, ['"' parameter_names{j} '"'])
                    set(eqMgr, 'Equation', i, sprintf('"%s"= %.10g', parameter_names{j}, parameter(j)));
                    nSet = nSet + 1;
                    break
                end
            end
        end
        if nSet ~= numel(parameter_names)
            fprintf('FAIL: updated %d of %d parameters\n', nSet, numel(parameter_names)); return
        end

        if ~invoke(swModel, 'ForceRebuild3', true)
            fprintf('FAIL: rebuild returned false\n'); return
        end

        mp = invoke(get(swModel, 'Extension'), 'CreateMassProperty');
        cg = get(mp, 'CenterOfMass');
        if isempty(cg), fprintf('FAIL: no center of mass\n'); return, end

        invoke(swModel, 'SaveAs3', outputFile, 0, 1);
        if ~isfile(outputFile), fprintf('FAIL: STL was not written\n'); return, end
    catch ME
        fprintf('FAIL (exception): %s\n', ME.message);
        return
    end

    TR = stlread(outputFile);
    ax = TR.Points(:, AX);
    cg_frac = (cg(AX)-min(ax) ) / (max(ax) - min(ax));     % 0 = nose, assuming nose at max end

    [Cl,Cd,~,x_cp_frac] = blended_driver(outputFile, alpha, minf, sref_lift, sref_drag, y);
    if ~isfinite(Cl) || ~isfinite(Cd) || Cd <= 0 || ~isfinite(x_cp_frac)
        fprintf('FAIL: bad solve  Cl=%g Cd=%g x_cp=%g\n', Cl, Cd, x_cp_frac); return
    end
    disp(parameter)
    fprintf('cg_frac=%.3f  x_cp_frac=%.3f  Cl=%.3g  Cd=%.3g  min_cl=%.4f\n', cg_frac, x_cp_frac, Cl, Cd, min_cl);
    fprintf('STL axial %.3f to %.3f m,  CG axial %.3f m\n', min(ax), max(ax), cg(AX));
    out.Fval = -Cl/Cd;
    out.Ineq = [cg_frac - x_cp_frac;
                min_cl - Cl];
end













% 
% 
% function out = obj(x,alpha,minf,sref_lift,sref_drag,y,min_cl)
% 
%         out.Fval = 1e3;
%         out.Ineq = [1;1];
% 
% 
% 
%     parameter_names = {'bodyLG','noseLG','bodyWD','finWD','bodyHT'};    
%     heightmax = .2+sqrt(max(0,.2^2-(x(3)/2+x(4))^2));
%     height = .15+x(3)*(heightmax-.15);
%     parameter= [x(1), x(2), x(3), x(4), height];
% 
%     %creates cad
%     inputFile  = 'N:\dotfiles\Desktop\euchre\virginia.SLDPRT'; 
%     outputFile = 'N:\dotfiles\Desktop\euchre\virginia.STL';   
%     try
%         % Hook into or launch SolidWorks
%         swApp = actxserver('SldWorks.Application');
%         set(swApp, 'Visible', false); % Set to false if you want it backgrounded
%     catch ME
%         error('Could not connect to SolidWorks: %s', ME.message);
%     end
%     docType = 1; % 1 for Part (.sldprt), 2 for Assembly (.sldasm)
%     options = 1; % 1 = Silent open (suppresses dialog boxes)
% 
%     fprintf('Opening file: %s...\n', inputFile);
%     swModel = swApp.OpenDoc6(inputFile, docType, options, '', 0, 0);
% 
%     if isempty(swModel)
%         error('Failed to open file. Verify the file path.');
%     end
% 
%     % 4. Access the Equation Manager to modify the global parameter
%     eqMgr = swModel.GetEquationMgr();
%     numEquations = eqMgr.GetCount();
% 
% 
% for i = 0:(numEquations - 1)
%     eqText = strtrim(eqMgr.Equation(i));          % e.g. "bodyLG"= 1.25
%     for j = 1:length(parameter_names)
%         if startsWith(eqText, ['"' parameter_names{j} '"'])
%             eqString = sprintf('"%s" = %f', parameter_names{j}, parameter(j));
%             eqMgr.Equation(i, eqString);
%             break
%         end
%     end
% end
% 
%     % 5. Force a full rebuild of the model geometry
%     % Rebuild options: 1 = Force rebuild specified, 2 = Force rebuild all
%     fprintf('Rebuilding geometry...\n');
%     swModel.ForceRebuild3(true); 
% 
%     % 1. Create a Mass Properties object from the active model
%     swMassProps = swModel.Extension.CreateMassProperty();
% 
%     % 2. Retrieve the Center of Gravity coordinates
%     % CenterOfMass returns a 3-element array of doubles: [X, Y, Z]
%     % IMPORTANT: SolidWorks API ALWAYS returns these coordinates in METERS.
%     cgCoordinates = swMassProps.CenterOfMass;
% 
%     if ~isempty(cgCoordinates)
%         % Extract the X coordinate (MATLAB is 1-indexed, so index 1 is X)
%         cgX = cgCoordinates(1);
% 
% 
%     else
%         warning('Could not calculate Center of Gravity. Ensure the part has a material assigned or density.');
%     end
% 
%     TR = stlread(outputFile);
%     axial_pts = TR.Points(:,1);                   
%     cg_frac = (cgX - min(axial_pts)) / (max(axial_pts) - min(axial_pts));
% 
%     % 6. Export to STL file
%     % SaveAs parameters: NewName, Version (0), Options (1 = Silent), Errors, Warnings
%     fprintf('Exporting to STL: %s...\n', outputFile);
%     saveSuccess = swModel.Extension.SaveAs(outputFile, 0, 1, [], 0, 0);
% 
%     if saveSuccess
%         fprintf('Export successful!\n');
%     else
%         warning('Export failed. Check folder write permissions.');
%     end
% 
%     % 7. Clean up connection (closes the document, keeps SolidWorks open)
%      swApp.CloseDoc(swModel.GetTitle());
%     swApp.ExitApp();
%     release(swApp);
% 
%     %now the actual stuff
%     %needs to have cp behind cg and needs to have a min cl to not fall out
%     %of sky
%     [Cl,Cd,~,x_cp]  = blended_driver(outputFile, alpha, minf, sref_lift, sref_drag,y);
%     if ~isfinite(Cl) || ~isfinite(Cd) || Cd <= 0 || ~isfinite(x_cp)
%         return          % bad solve -> keep the defaults
%     end
% 
%     out.Fval = -Cl/Cd;
%     out.Ineq = [cg_frac - x_cp;     % <= 0 when CP is aft of CG
%                 min_cl - Cl];       % <= 0 when lift is sufficient
% 
% end
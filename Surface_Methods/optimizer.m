clc;
clear
close all

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

%             'bodyLG','noseLG' 'bodyWD','finWD', height fraction',];
parameters_0 = [1.25     .74     .2        .1            .75 ];

%%
%need to have lb and ub for the parameters
lb = [.5 .25 .15 .05 0 ];

ub = [1.5 1.5 .3 .125 1];

A = [1 1 0 0 0 ; %maxes out length
     0 0 1 2 0; %maxes out width
   
    

];

B= [2;.4
    
];



%%
%need to have an optimizer surragote





options = optimoptions('surrogateopt', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 300, ...
    'UseParallel', false, ...
    'PlotFcn', 'surrogateoptplot', ...
    'InitialPoints', parameters_0);     

fun = @(x) obj(x, alpha, M, sref_lift, sref_drag, y, min_cl);

[optimal_shape, fval, exitflag, output, trials] = surrogateopt(fun, lb, ub, [], A, B, [], [], options);
%%






%need to have an objective function which takes inputs, creates cad file,
%runs the driver, and then spits out shit

function f = obj(x,alpha,minf,sref_lift,sref_drag,y,min_cl)
    parameter_names = {'bodyLG','noseLG','bodyWD','finWD','bodyHT'};    
    heightmax = 1.5+sqrt(max(0,.2^2-(x(2)/2+x(5))^2));
    height = .15+x(3)*(heightmax-.15);
    parameter= [x(1), x(2), x(3), x(4), height];

    %creates cad
    inputFile  = 'N:\dotfiles\Desktop\euchre\virginia.SLDPRT'; 
    outputFile = 'N:\dotfiles\Desktop\euchre\virginia.STL';   
    try
        % Hook into or launch SolidWorks
        swApp = actxserver('SldWorks.Application');
        set(swApp, 'Visible', false); % Set to false if you want it backgrounded
    catch ME
        error('Could not connect to SolidWorks: %s', ME.message);
    end
    docType = 1; % 1 for Part (.sldprt), 2 for Assembly (.sldasm)
    options = 1; % 1 = Silent open (suppresses dialog boxes)
    
    fprintf('Opening file: %s...\n', inputFile);
    swModel = swApp.OpenDoc6(inputFile, docType, options, '', 0, 0);
    
    if isempty(swModel)
        error('Failed to open file. Verify the file path.');
    end
    
    % 4. Access the Equation Manager to modify the global parameter
    eqMgr = swModel.GetEquationMgr();
    numEquations = eqMgr.GetCount();
 

for i = 0:(numEquations - 1)
    eqText = strtrim(eqMgr.Equation(i));          % e.g. "bodyLG"= 1.25
    for j = 1:length(parameter_names)
        if startsWith(eqText, ['"' parameter_names{j} '"'])
            eqString = sprintf('"%s" = %f', parameter_names{j}, parameter(j));
            eqMgr.Equation(i, eqString);
            break
        end
    end
end

    % 5. Force a full rebuild of the model geometry
    % Rebuild options: 1 = Force rebuild specified, 2 = Force rebuild all
    fprintf('Rebuilding geometry...\n');
    swModel.ForceRebuild3(true); 
    
    % 1. Create a Mass Properties object from the active model
    swMassProps = swModel.Extension.CreateMassProperty();
    
    % 2. Retrieve the Center of Gravity coordinates
    % CenterOfMass returns a 3-element array of doubles: [X, Y, Z]
    % IMPORTANT: SolidWorks API ALWAYS returns these coordinates in METERS.
    cgCoordinates = swMassProps.CenterOfMass;
    
    if ~isempty(cgCoordinates)
        % Extract the X coordinate (MATLAB is 1-indexed, so index 1 is X)
        cg = cgCoordinates(1);
        
        % Convert to millimeters for standard engineering use
        cgX_mm = cg * 1000; 
       
    else
        warning('Could not calculate Center of Gravity. Ensure the part has a material assigned or density.');
    end


    % 6. Export to STL file
    % SaveAs parameters: NewName, Version (0), Options (1 = Silent), Errors, Warnings
    fprintf('Exporting to STL: %s...\n', outputFile);
    saveSuccess = swModel.Extension.SaveAs(outputFile, 0, 1, [], 0, 0);
    
    if saveSuccess
        fprintf('Export successful!\n');
    else
        warning('Export failed. Check folder write permissions.');
    end
    
    % 7. Clean up connection (closes the document, keeps SolidWorks open)
     swApp.CloseDoc(swModel.GetTitle());
    swApp.ExitApp();
    release(swApp);




    %now the actual stuff
    %needs to have cp behind cg and needs to have a min cl to not fall out
    %of sky
    [Cl,Cd,~,x_cp]  = blended_driver(outputFile, alpha, minf, sref_lift, sref_drag,y);
    
    f = -Cl/Cd;
    if x_cp<cgX_mm
        f=1e3;
    end
    if Cl< min_cl
       f=1e3;
    end

end
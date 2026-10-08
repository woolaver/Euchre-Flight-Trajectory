function results = hypersonicHeatFlux(traj, geom, mat, opts)
%HYPERSONICHEATFLUX  Convective heat flux on a hypersonic body via Eckert's
%Reference Temperature Method, with 2-D wedge oblique-shock edge conditions
%(exact theta-beta-M solution).
%
%   results = HYPERSONICHEATFLUX(traj, geom, mat, opts)
%
%   TWO ANALYSIS MODES (select with opts.mode)
%   ---------------------------------------------------------------
%   'time'      (default) Heat flux, wall temperature, etc. vs. TIME at a
%               single downstream station geom.x.
%   'position'  Heat flux, wall temperature, etc. vs. DOWNSTREAM POSITION
%               (geom.xStations) at a single snapshot time
%               (opts.snapshotTime). The lumped wall temperature at each
%               station is found by integrating that station's wall energy
%               balance from the start of the trajectory up to the
%               snapshot time, so Tw(x) reflects each station's own
%               heating history.
%
%   METHOD
%   ---------------------------------------------------------------
%   1. Freestream state from the 1976 U.S. Standard Atmosphere.
%   2. Boundary-layer-edge state (Me, Te, pe) from an attached oblique
%      shock on a 2-D wedge of half-angle geom.halfAngle (+|alpha| as a
%      crude AoA placeholder). Leave geom.halfAngle empty ([]) to use
%      freestream conditions (flat-plate limit).
%   3. Eckert reference temperature T*, local Reynolds number at T*,
%      laminar (Blasius) or turbulent (Schlichting) skin friction.
%   4. Reynolds-Colburn analogy gives Stanton number, then
%      q_w = rho* Ue Cp St (Taw - Tw).
%   5. Wall is a single lumped node: rho*t*cp*dTw/dt = q_conv - q_rad,
%      integrated with ode45.
%
%   REQUIRED INPUT: traj (struct)
%   ---------------------------------------------------------------
%     traj.time      [s]     time vector, strictly increasing
%     traj.altitude  [m]     altitude vector, same length as time
%     traj.velocity  [m/s]   velocity vector, same length as time
%     traj.alpha     [deg]   (optional) angle of attack, scalar or same
%                            length as time. Default: 0.
%
%   OPTIONAL INPUT: geom (struct)
%   ---------------------------------------------------------------
%     geom.x           [m]   station for mode 'time'. Default: 0.5
%     geom.xStations   [m]   vector of stations for mode 'position'.
%                            Default: linspace(0.1, 1, 10)
%     geom.halfAngle   [deg] wedge half-angle, [] for flat-plate limit.
%                            Default: 10
%     geom.emissivity  [-]   wall emissivity. Default: 0.80
%
%   OPTIONAL INPUT: mat (struct)
%   ---------------------------------------------------------------
%     mat.thickness  [m]        Default: 0.003
%     mat.density    [kg/m^3]   Default: 8000
%     mat.cp         [J/kg-K]   Default: 500
%     mat.Tw0        [K]        initial wall temperature. Default: 290
%
%   OPTIONAL INPUT: opts (struct)
%   ---------------------------------------------------------------
%     opts.mode             'time' (default) or 'position'
%     opts.snapshotTime     [s] time at which the position profile is
%                           evaluated (mode 'position'). Must lie within
%                           traj.time. Default: last trajectory time.
%     opts.gamma            Default: 1.4
%     opts.R                [J/kg-K] Default: 287
%     opts.Pr               Default: 0.71
%     opts.Re_trans         reference-condition Reynolds number above which
%                           the boundary layer is turbulent. Default: 5e5
%     opts.includeRadiation include radiative wall cooling. Default: true
%     opts.makePlots        generate plots. Default: true
%
%   OUTPUT: results (struct)
%   ---------------------------------------------------------------
%   mode 'time'  -- fields are [N x 1] time histories:
%     time, altitude, velocity, Mach_inf, Me, Te, Taw, Tstar, Tw, qw,
%     qrad, Cf, regime
%
%   mode 'position' -- fields are [Nx x 1] profiles along the body:
%     x, Taw, Tstar, Tw, qw, qrad, Cf, regime
%   plus scalars (independent of x for the wedge edge model):
%     time, altitude, velocity, Mach_inf, Me, Te
%
%   Units of qw and qrad are W/m^2 (plots are shown in W/cm^2).
%
%   NOTES / SIMPLIFICATIONS
%   ---------------------------------------------------------------
%   - A 2-D wedge shock is an approximation for an axisymmetric cone.
%   - Atmosphere model is valid to 86 km; above that it is approximate.
%   - Because the edge state is the same at every station, only the
%     boundary-layer quantities (T*, Cf, qw, Tw, regime) vary with x.
%
%   Examples:
%       % Time history at x = 0.4 m
%       results = hypersonicHeatFlux(traj, geom, mat, struct('mode','time'));
%
%       % Profile from 0.1 to 1 m at t = 25 s
%       geom.xStations = linspace(0.1, 1, 30);
%       o.mode = 'position'; o.snapshotTime = 25;
%       results = hypersonicHeatFlux(traj, geom, mat, o);

% =========================================================================
% 1. INPUT HANDLING / DEFAULTS
% =========================================================================
if nargin < 2 || isempty(geom), geom = struct(); end
if nargin < 3 || isempty(mat),  mat  = struct(); end
if nargin < 4 || isempty(opts), opts = struct(); end

t_traj = traj.time(:);
h_traj = traj.altitude(:);
V_traj = traj.velocity(:);
if isfield(traj, 'alpha') && ~isempty(traj.alpha)
    if isscalar(traj.alpha)
        alpha_traj = repmat(traj.alpha, size(t_traj));
    else
        alpha_traj = traj.alpha(:);
    end
else
    alpha_traj = zeros(size(t_traj));
end

if ~isequal(numel(t_traj), numel(h_traj), numel(V_traj), numel(alpha_traj))
    error('hypersonicHeatFlux:sizeMismatch', ...
        'traj.time, traj.altitude, traj.velocity, and traj.alpha must be the same length.');
end
if any(diff(t_traj) <= 0)
    error('hypersonicHeatFlux:timeNotMonotonic', 'traj.time must be strictly increasing.');
end

geom = setDefault(geom, 'x',           0.5);
geom = setDefault(geom, 'xStations',   linspace(0.1, 1, 10));
geom = setDefault(geom, 'halfAngle',   10);
geom = setDefault(geom, 'emissivity',  0.80);

mat  = setDefault(mat, 'thickness', 0.003);
mat  = setDefault(mat, 'density',   8000);
mat  = setDefault(mat, 'cp',        500);
mat  = setDefault(mat, 'Tw0',       290);

opts = setDefault(opts, 'mode',              'time');
opts = setDefault(opts, 'snapshotTime',      t_traj(end));
opts = setDefault(opts, 'gamma',             1.4);
opts = setDefault(opts, 'R',                 287);
opts = setDefault(opts, 'Pr',                0.71);
opts = setDefault(opts, 'Re_trans',          5e5);
opts = setDefault(opts, 'includeRadiation',  true);
opts = setDefault(opts, 'makePlots',         true);

% Shared context passed to the helper functions
ctx.geom    = geom;
ctx.mat     = mat;
ctx.opts    = opts;
ctx.sigmaSB = 5.670374419e-8; % Stefan-Boltzmann constant, W/m^2-K^4
ctx.hOf     = @(tq) interp1(t_traj, h_traj,     tq, 'linear', 'extrap');
ctx.VOf     = @(tq) interp1(t_traj, V_traj,     tq, 'linear', 'extrap');
ctx.alphaOf = @(tq) interp1(t_traj, alpha_traj, tq, 'linear', 'extrap');
ctx.odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-3);

% =========================================================================
% 2. SOLVE + PLOT, by mode
% =========================================================================
switch lower(opts.mode)
    case 'time'
        results = solveTimeHistory(ctx, t_traj);
        if opts.makePlots
            plotTimeHistory(results, geom);
        end

    case 'position'
        if opts.snapshotTime < t_traj(1) || opts.snapshotTime > t_traj(end)
            error('hypersonicHeatFlux:badSnapshotTime', ...
                'opts.snapshotTime (%.4g s) must lie within traj.time [%.4g, %.4g] s.', ...
                opts.snapshotTime, t_traj(1), t_traj(end));
        end
        if any(geom.xStations <= 0)
            error('hypersonicHeatFlux:badStations', 'geom.xStations must all be > 0.');
        end
        results = solvePositionProfile(ctx, t_traj(1));
        if opts.makePlots
            plotPositionProfile(results);
        end

    otherwise
        error('hypersonicHeatFlux:badMode', ...
            'opts.mode must be ''time'' or ''position'' (got ''%s'').', opts.mode);
end

end % ================== END MAIN FUNCTION ==================


% =========================================================================
% LOCAL FUNCTION: MODE 'time' -- histories at one station geom.x
% =========================================================================
function results = solveTimeHistory(ctx, t_traj)

x = ctx.geom.x;
opts = ctx.opts;

odefun = @(t, Tw) wallEnergyBalance(t, Tw, x, ctx);
[tSol, TwSol] = ode45(odefun, t_traj, ctx.mat.Tw0, ctx.odeOpts);

N = numel(tSol);
altOut = zeros(N,1); velOut = zeros(N,1); Mach_inf = zeros(N,1);
Me = zeros(N,1); Te = zeros(N,1);
qw = zeros(N,1); Cf = zeros(N,1); Taw = zeros(N,1); Tstar = zeros(N,1);
qrad = zeros(N,1);
regimeCell = cell(N,1);

for k = 1:N
    st  = flightState(tSol(k), ctx);
    Twk = TwSol(k);

    altOut(k)   = st.h;
    velOut(k)   = st.V;
    Mach_inf(k) = st.Minf;
    Me(k)       = st.Me;
    Te(k)       = st.Te;

    [qw(k), Cf(k), Taw(k), Tstar(k), regimeCell{k}] = ...
        referenceTempHeating(st.Me, st.Te, st.pe, st.Ue, Twk, x, opts);

    if opts.includeRadiation
        qrad(k) = ctx.geom.emissivity * ctx.sigmaSB * (Twk^4 - st.Tinf^4);
    end
end

results.time     = tSol;
results.altitude = altOut;
results.velocity = velOut;
results.Mach_inf = Mach_inf;
results.Me       = Me;
results.Te       = Te;
results.Taw      = Taw;
results.Tstar    = Tstar;
results.Tw       = TwSol;
results.qw       = qw;
results.qrad     = qrad;
results.Cf       = Cf;
results.regime   = regimeCell;

end


% =========================================================================
% LOCAL FUNCTION: MODE 'position' -- profile along the body at one time
% =========================================================================
function results = solvePositionProfile(ctx, t0)

opts  = ctx.opts;
xs    = ctx.geom.xStations(:);
tSnap = opts.snapshotTime;
Nx    = numel(xs);

% Edge state is independent of x (2-D wedge), so compute it once
st = flightState(tSnap, ctx);

Tw = zeros(Nx,1); qw = zeros(Nx,1); Cf = zeros(Nx,1);
Taw = zeros(Nx,1); Tstar = zeros(Nx,1); qrad = zeros(Nx,1);
regimeCell = cell(Nx,1);

for i = 1:Nx
    % Wall temperature at this station at the snapshot time
    if tSnap > t0
        odefun = @(t, Twv) wallEnergyBalance(t, Twv, xs(i), ctx);
        [~, TwSol] = ode45(odefun, [t0 tSnap], ctx.mat.Tw0, ctx.odeOpts);
        Tw(i) = TwSol(end);
    else
        Tw(i) = ctx.mat.Tw0;
    end

    [qw(i), Cf(i), Taw(i), Tstar(i), regimeCell{i}] = ...
        referenceTempHeating(st.Me, st.Te, st.pe, st.Ue, Tw(i), xs(i), opts);

    if opts.includeRadiation
        qrad(i) = ctx.geom.emissivity * ctx.sigmaSB * (Tw(i)^4 - st.Tinf^4);
    end
end

results.x        = xs;
results.time     = tSnap;
results.altitude = st.h;
results.velocity = st.V;
results.Mach_inf = st.Minf;
results.Me       = st.Me;
results.Te       = st.Te;
results.Taw      = Taw;
results.Tstar    = Tstar;
results.Tw       = Tw;
results.qw       = qw;
results.qrad     = qrad;
results.Cf       = Cf;
results.regime   = regimeCell;

end


% =========================================================================
% LOCAL FUNCTION: plots for mode 'time'
% =========================================================================
function plotTimeHistory(res, geom)

t = res.time;

figure('Name', 'Trajectory - Altitude vs Time');
plot(t, res.altitude/1000, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Altitude [km]');
title('Trajectory: Altitude vs. Time');
grid on;

figure('Name', 'Trajectory - Velocity vs Time');
plot(t, res.velocity, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Trajectory: Velocity vs. Time');
grid on;

figure('Name', 'Mach Number vs Time');
plot(t, res.Mach_inf, 'LineWidth', 1.6); hold on;
plot(t, res.Me, '--', 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Freestream vs. Boundary-Layer-Edge Mach Number vs. Time');
legend('Freestream, M_\infty', 'Boundary-layer edge, M_e', 'Location', 'best');
grid on; hold off;

figure('Name', 'Convective Heat Flux vs Time');
plot(t, res.qw/1e4, 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Convective Heat Flux [W/cm^2]');
title(sprintf('Convective Heat Flux vs. Time (x = %.3g m)', geom.x));
grid on;

figure('Name', 'Temperatures vs Time');
plot(t, res.Tw,    'LineWidth', 1.9); hold on;
plot(t, res.Taw,   '--', 'LineWidth', 1.6);
plot(t, res.Te,    ':',  'LineWidth', 1.6);
plot(t, res.Tstar, '-.', 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Temperature [K]');
title(sprintf('Wall, Recovery, Edge, and Reference Temperature vs. Time (x = %.3g m)', geom.x));
legend('Wall, T_w', 'Adiabatic wall, T_{aw}', 'Edge, T_e', ...
       'Eckert reference, T^*', 'Location', 'best');
grid on; hold off;

figure('Name', 'Skin Friction Coefficient vs Time');
plot(t, res.Cf, 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Skin Friction Coefficient, C_f [-]');
title(sprintf('Local Skin Friction Coefficient vs. Time (x = %.3g m)', geom.x));
grid on;

end


% =========================================================================
% LOCAL FUNCTION: plots for mode 'position'
% =========================================================================
function plotPositionProfile(res)

x = res.x;
tag = sprintf('t = %.3g s, alt = %.3g km, M_\\infty = %.3g', ...
    res.time, res.altitude/1000, res.Mach_inf);

figure('Name', 'Convective Heat Flux vs Position');
plot(x, res.qw/1e4, '-o', 'LineWidth', 1.8);
xlabel('Downstream Position, x [m]'); ylabel('Convective Heat Flux [W/cm^2]');
title(sprintf('Convective Heat Flux vs. Position (%s)', tag));
grid on;

figure('Name', 'Temperatures vs Position');
plot(x, res.Tw,    '-o', 'LineWidth', 1.9); hold on;
plot(x, res.Taw,   '--', 'LineWidth', 1.6);
plot(x, res.Tstar, '-.', 'LineWidth', 1.6);
plot(x, res.Te*ones(size(x)), ':', 'LineWidth', 1.6);
xlabel('Downstream Position, x [m]'); ylabel('Temperature [K]');
title(sprintf('Wall, Recovery, Reference, and Edge Temperature vs. Position (%s)', tag));
legend('Wall, T_w', 'Adiabatic wall, T_{aw}', 'Eckert reference, T^*', ...
       'Edge, T_e', 'Location', 'best');
grid on; hold off;

figure('Name', 'Skin Friction Coefficient vs Position');
plot(x, res.Cf, '-o', 'LineWidth', 1.8);
xlabel('Downstream Position, x [m]'); ylabel('Skin Friction Coefficient, C_f [-]');
title(sprintf('Local Skin Friction Coefficient vs. Position (%s)', tag));
grid on;

end


% =========================================================================
% LOCAL FUNCTION: freestream + boundary-layer-edge state at time t
% =========================================================================
function st = flightState(t, ctx)

gam = ctx.opts.gamma;
R   = ctx.opts.R;

h     = ctx.hOf(t);
V     = ctx.VOf(t);
alpha = ctx.alphaOf(t);

[~, T_inf, p_inf, ~] = standardAtmosphere(h);
M_inf = V / sqrt(gam * R * T_inf);

if ~isempty(ctx.geom.halfAngle)
    thetaEff = ctx.geom.halfAngle + abs(alpha); % crude AoA placeholder, deg
    [Me, Te, pe, ~] = obliqueShockEdge(M_inf, T_inf, p_inf, thetaEff, gam);
else
    Me = M_inf; Te = T_inf; pe = p_inf;
end

st.h    = h;
st.V    = V;
st.Tinf = T_inf;
st.Minf = M_inf;
st.Me   = Me;
st.Te   = Te;
st.pe   = pe;
st.Ue   = Me * sqrt(gam * R * Te);

end


% =========================================================================
% LOCAL FUNCTION: single-node wall energy balance (ODE right-hand side)
% =========================================================================
function dTwdt = wallEnergyBalance(t, Tw, x, ctx)

st = flightState(t, ctx);

qw = referenceTempHeating(st.Me, st.Te, st.pe, st.Ue, Tw, x, ctx.opts);

if ctx.opts.includeRadiation
    qrad = ctx.geom.emissivity * ctx.sigmaSB * (Tw^4 - st.Tinf^4);
else
    qrad = 0;
end

massArea = ctx.mat.density * ctx.mat.thickness; % kg/m^2
dTwdt = (qw - qrad) / (massArea * ctx.mat.cp);

end


% =========================================================================
% LOCAL FUNCTION: Eckert reference temperature method
% =========================================================================
function [qw, Cf, Taw, Tstar, regime] = referenceTempHeating(Me, Te, pe, Ue, Tw, x, opts)

gamma = opts.gamma;
R     = opts.R;
Pr    = opts.Pr;

% Eckert reference temperature: T*/Te = 1 + 0.032*Me^2 + 0.58*(Tw/Te - 1)
Tstar = Te * (1 + 0.032*Me^2 + 0.58*(Tw/Te - 1));
[rho_star, mu_star] = referenceGasProps(Tstar, pe, R);
Rex_star = rho_star * Ue * x / mu_star;

%%%% hardcoding turb
%if Rex_star <= opts.Re_trans
%    regime = 'laminar';
%    r  = sqrt(Pr);                 % laminar recovery factor
%    Cf = 0.664 / sqrt(Rex_star);   % Blasius, evaluated at T*
%else
    regime = 'turbulent';
    r  = Pr^(1/3);                 % turbulent recovery factor
    Cf = 0.0592 / Rex_star^0.2;    % Schlichting flat-plate, at T*
%end

Taw = Te * (1 + r * (gamma - 1)/2 * Me^2);

% Reynolds-Colburn analogy
Cp_star = gamma * R / (gamma - 1);
St_star = (Cf/2) / Pr^(2/3);

qw = rho_star * Ue * Cp_star * St_star * (Taw - Tw);

end


% =========================================================================
% LOCAL FUNCTION: reference density & viscosity at T*
% =========================================================================
function [rho_star, mu_star] = referenceGasProps(Tstar, pe, R)

rho_star = pe / (R * Tstar);

mu0 = 1.716e-5;   % kg/(m-s), Sutherland reference viscosity
T0  = 273.15;     % K
S   = 110.4;      % K
mu_star = mu0 * (Tstar/T0)^1.5 * (T0 + S) / (Tstar + S);

end


% =========================================================================
% LOCAL FUNCTION: exact oblique-shock boundary-layer-edge conditions
% =========================================================================
function [M2, T2, p2, rho2] = obliqueShockEdge(M1, T1, p1, thetaDeg, gamma)

R_air = 287;
theta = deg2rad(thetaDeg);

if theta <= 0 || M1 <= 1
    M2 = M1; T2 = T1; p2 = p1; rho2 = p1/(R_air*T1);
    return;
end

muMach = asin(1/M1); % Mach angle, lower bound for beta

% Exact theta-beta-M relation, solved for the weak-shock angle beta
thetaBetaM = @(beta) tan(theta) - 2*cot(beta).* ...
    (M1^2 .* sin(beta).^2 - 1) ./ (M1^2 .* (gamma + cos(2*beta)) + 2);

betaLow  = muMach*1.0001;
betaHigh = pi/2 - 1e-6;

if thetaBetaM(betaLow) * thetaBetaM(betaHigh) > 0
    % Detached shock: fall back to freestream
    M2 = M1; T2 = T1; p2 = p1; rho2 = p1/(R_air*T1);
    return;
end

beta = fzero(thetaBetaM, [betaLow, betaHigh]);

M1n = M1 * sin(beta);
p2_p1     = 1 + 2*gamma/(gamma+1) * (M1n^2 - 1);
rho2_rho1 = (gamma+1) * M1n^2 / ((gamma-1) * M1n^2 + 2);
T2_T1     = p2_p1 / rho2_rho1;

M2n = sqrt( (1 + (gamma-1)/2 * M1n^2) / (gamma*M1n^2 - (gamma-1)/2) );
M2  = M2n / sin(beta - theta);

p2   = p1 * p2_p1;
T2   = T1 * T2_T1;
rho2 = rho2_rho1 * (p1/(R_air*T1));

end


% =========================================================================
% LOCAL FUNCTION: 1976 U.S. Standard Atmosphere (valid to 86 km geometric)
% =========================================================================
function [Hgp, T, p, rho] = standardAtmosphere(h_geometric)

Re = 6356766;    % m
g0 = 9.80665;    % m/s^2
Rgas = 287.0528; % J/kg-K

Hgp = Re * h_geometric / (Re + h_geometric);

baseH = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];
baseT = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946];
lapse = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, 0];
baseP = zeros(1, numel(baseH));
baseP(1) = 101325;

for i = 2:numel(baseH)
    dH = baseH(i) - baseH(i-1);
    L  = lapse(i-1);
    Tb = baseT(i-1);
    Pb = baseP(i-1);
    if abs(L) < 1e-12
        baseP(i) = Pb * exp(-g0*dH/(Rgas*Tb));
    else
        baseP(i) = Pb * (1 + L*dH/Tb)^(-g0/(Rgas*L));
    end
end

if Hgp > baseH(end)
    Tref = baseT(end);
    Pref = baseP(end);
    scaleHeight = Rgas * Tref / g0;
    T = Tref;
    p = Pref * exp(-(Hgp - baseH(end)) / scaleHeight);
    rho = p / (Rgas * T);
    return;
end

idx = find(Hgp >= baseH, 1, 'last');
dH = Hgp - baseH(idx);
L  = lapse(idx);
Tb = baseT(idx);
Pb = baseP(idx);

T = Tb + L*dH;
if abs(L) < 1e-12
    p = Pb * exp(-g0*dH/(Rgas*Tb));
else
    p = Pb * (1 + L*dH/Tb)^(-g0/(Rgas*L));
end
rho = p / (Rgas * T);

end


% =========================================================================
% LOCAL FUNCTION: fill in default struct fields
% =========================================================================
function s = setDefault(s, field, value)
if ~isfield(s, field)
    s.(field) = value;
end
end

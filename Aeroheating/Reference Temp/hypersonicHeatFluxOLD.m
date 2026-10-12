function results = hypersonicHeatFlux(traj, geom, mat, opts)
%HYPERSONICHEATFLUX  Convective heat flux at a downstream body station
%over a hypersonic flight trajectory, via Eckert's Reference Temperature
%Method, with the body approximated as a 2-D wedge and the oblique-shock
%boundary-layer-edge conditions computed from the EXACT theta-beta-M
%relation (solved numerically, not a small-angle/linearized approximation).
%
%   results = HYPERSONICHEATFLUX(traj, geom, mat, opts)
%
%   This is a deliberately simplified tool focused on one question:
%   "What is the convective heat flux at a given downstream location on
%   the body, as a function of time, over this flight?" Everything else
%   (wall temperature, skin friction, recovery temperature) is carried
%   along only because the reference-temperature method needs them.
%
%   METHOD
%   ---------------------------------------------------------------
%   1. Freestream state comes from the 1976 U.S. Standard Atmosphere.
%   2. The body is approximated as a 2-D wedge of half-angle
%      geom.halfAngle. The boundary-layer-edge state (Me, Te, pe) is the
%      state behind an attached oblique shock at that wedge angle,
%      computed by solving the exact theta-beta-M relation for the shock
%      angle (via fzero), then applying the exact oblique-shock jump
%      relations -- no weak-shock/small-angle simplification is used. If
%      geom.halfAngle is left empty, the edge state is simply set equal
%      to the freestream (flat-plate limit).
%   3. At the station geom.x (running length from the nose/leading
%      edge), the Eckert reference temperature T* is computed, and the
%      local Reynolds number at T* is used to decide laminar vs.
%      turbulent and to evaluate the skin-friction coefficient (Blasius
%      or Schlichting flat-plate correlations, evaluated with reference-
%      temperature properties).
%   4. The Reynolds(-Colburn) analogy converts skin friction to a Stanton
%      number, giving the convective heat flux q_w = rho* Ue Cp St (Taw - Tw).
%   5. The wall is a single lumped-capacitance node (uniform through the
%      thickness): rho*t*cp*dTw/dt = q_conv - q_rad, integrated with
%      ode45 over the supplied trajectory. This gives Tw(t) self-
%      consistently rather than assuming a fixed/cold wall.
%
%   ---------------------------------------------------------------
%   REQUIRED INPUT: traj (struct)
%   ---------------------------------------------------------------
%     traj.time      [s]     time vector, strictly increasing
%     traj.altitude  [m]     altitude vector, same length as time
%     traj.velocity  [m/s]   velocity magnitude vector, same length
%     traj.alpha     [deg]   (optional) angle of attack, scalar or same
%                            length as time. Default: 0. Added directly
%                            to geom.halfAngle as a crude effective wedge
%                            angle (theta_eff = halfAngle + |alpha|) --
%                            a placeholder, not a true windward-ray
%                            analysis.
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: geom (struct) -- reference/placeholder values
%   ---------------------------------------------------------------
%     geom.x           [m]   downstream station (running length from the
%                              nose/leading edge) at which the heat flux
%                              is evaluated. Default: 0.5
%     geom.halfAngle   [deg] wedge half-angle used for the oblique-shock
%                              edge-condition estimate. Leave empty ([])
%                              to skip the shock calculation and use
%                              freestream conditions directly (flat-plate
%                              limit). Default: 10
%     geom.emissivity  [-]   wall emissivity. Default: 0.80
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: mat (struct) -- lumped wall thermal properties
%   ---------------------------------------------------------------
%     mat.thickness  [m]        wall/skin thickness. Default: 0.003
%     mat.density    [kg/m^3]   wall material density. Default: 8000
%     mat.cp         [J/kg-K]   wall material specific heat. Default: 500
%     mat.Tw0        [K]        initial wall temperature. Default: 290
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: opts (struct) -- gas/model parameters
%   ---------------------------------------------------------------
%     opts.gamma            [-]  ratio of specific heats. Default: 1.4
%     opts.R                [J/kg-K] specific gas constant, air. Default: 287
%     opts.Pr               [-]  Prandtl number (assumed constant). Default: 0.71
%     opts.Re_trans         [-]  local (reference-condition) Reynolds
%                                 number above which the boundary layer
%                                 is treated as turbulent. Default: 5e5
%     opts.includeRadiation (logical) include radiative wall cooling in
%                                 the transient wall energy balance.
%                                 Default: true
%     opts.makePlots        (logical) generate the plot suite described
%                                 below. Default: true
%
%   ---------------------------------------------------------------
%   OUTPUT: results (struct), all fields are time histories [N x 1]
%   ---------------------------------------------------------------
%     results.time     [s]
%     results.altitude [m]
%     results.velocity [m/s]
%     results.Mach_inf [-]      freestream Mach number
%     results.Me       [-]      boundary-layer edge Mach number
%     results.Te       [K]      boundary-layer edge static temperature
%     results.Taw      [K]      adiabatic wall (recovery) temperature
%     results.Tstar    [K]      Eckert reference temperature
%     results.Tw       [K]      wall temperature (transient solution)
%     results.qw       [W/m^2]  convective heat flux to the wall
%     results.qrad      [W/m^2]  radiative heat flux leaving the wall
%     results.Cf        [-]      local skin-friction coefficient
%     results.regime    {'laminar'|'turbulent'} cell array, per time step
%
%   ---------------------------------------------------------------
%   NOTES / SIMPLIFICATIONS
%   ---------------------------------------------------------------
%   - The 2-D wedge/oblique-shock model is an approximation for a true
%     axisymmetric (conical) body -- a conical shock runs at a shallower
%     angle than a 2-D wedge shock of the same deflection, so this will
%     tend to over-predict edge Mach loss / over-predict heating somewhat
%     for a true cone. It is a reasonable, simple first pass; swap in a
%     Taylor-Maccoll cone solution later if higher fidelity is needed.
%   - The oblique-shock angle itself IS solved exactly (via fzero on the
%     full theta-beta-M relation), so within the 2-D wedge idealization
%     there is no shock-strength approximation being made.
%   - The atmosphere model (1976 U.S. Standard Atmosphere) is valid to
%     86 km geometric altitude; above that it falls back to a simple
%     exponential extrapolation and should be treated as approximate.
%   - The wall is a single lumped node (uniform through the thickness) --
%     this captures the transient response of a thin skin but does not
%     distinguish an outer vs. inner surface temperature.
%
%   Example:
%       traj.time     = linspace(0, 60, 300);
%       traj.altitude = linspace(40000, 2000, 300);   % m
%       traj.velocity = linspace(3000, 1200, 300);    % m/s
%       results = hypersonicHeatFlux(traj, [], [], []);

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
geom = setDefault(geom, 'halfAngle',   10);
geom = setDefault(geom, 'emissivity',  0.80);

mat  = setDefault(mat, 'thickness', 0.003);
mat  = setDefault(mat, 'density',   8000);
mat  = setDefault(mat, 'cp',        500);
mat  = setDefault(mat, 'Tw0',       290);

opts = setDefault(opts, 'gamma',             1.4);
opts = setDefault(opts, 'R',                 287);
opts = setDefault(opts, 'Pr',                0.71);
opts = setDefault(opts, 'Re_trans',          5e5);
opts = setDefault(opts, 'includeRadiation',  true);
opts = setDefault(opts, 'makePlots',         true);

sigmaSB = 5.670374419e-8; % Stefan-Boltzmann constant, W/m^2-K^4

% Interpolants used inside the ODE integration (queried at arbitrary t)
hOf     = @(tq) interp1(t_traj, h_traj,     tq, 'linear', 'extrap');
VOf     = @(tq) interp1(t_traj, V_traj,     tq, 'linear', 'extrap');
alphaOf = @(tq) interp1(t_traj, alpha_traj, tq, 'linear', 'extrap');

% =========================================================================
% 2. TRANSIENT SOLVE: single-node wall temperature at station geom.x
% =========================================================================
odefun = @(t, Tw) wallEnergyBalance(t, Tw, hOf, VOf, alphaOf, geom, mat, opts, sigmaSB);
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-3);
[tSol, TwSol] = ode45(odefun, t_traj, mat.Tw0, odeOptions);

% =========================================================================
% 3. POST-PROCESS: recompute all quantities of interest on the solution grid
% =========================================================================
N = numel(tSol);

Mach_inf = zeros(N,1);
Me       = zeros(N,1);
Te       = zeros(N,1);
altOut   = zeros(N,1);
velOut   = zeros(N,1);
qw       = zeros(N,1);
Cf       = zeros(N,1);
Taw      = zeros(N,1);
Tstar    = zeros(N,1);
qrad     = zeros(N,1);
regimeCell = cell(N,1);

for k = 1:N
    tk = tSol(k);
    Twk = TwSol(k);

    hk     = hOf(tk);
    Vk     = VOf(tk);
    alphak = alphaOf(tk);

    altOut(k) = hk;
    velOut(k) = Vk;

    [~, T_inf, p_inf, ~] = standardAtmosphere(hk);
    a_inf = sqrt(opts.gamma * opts.R * T_inf);
    M_inf = Vk / a_inf;
    Mach_inf(k) = M_inf;

    if ~isempty(geom.halfAngle)
        thetaEff = geom.halfAngle + abs(alphak); % crude AoA placeholder, deg
        [Me_k, Te_k, pe_k, ~] = obliqueShockEdge(M_inf, T_inf, p_inf, thetaEff, opts.gamma);
    else
        Me_k = M_inf;
        Te_k = T_inf;
        pe_k = p_inf;
    end
    Me(k) = Me_k;
    Te(k) = Te_k;

    Ue_k = Me_k * sqrt(opts.gamma * opts.R * Te_k);

    [qw_k, Cf_k, Taw_k, Tstar_k, regime_k] = referenceTempHeating( ...
        Me_k, Te_k, pe_k, Ue_k, Twk, geom.x, opts);

    qw(k)    = qw_k;
    Cf(k)    = Cf_k;
    Taw(k)   = Taw_k;
    Tstar(k) = Tstar_k;
    regimeCell{k} = regime_k;

    if opts.includeRadiation
        qrad(k) = geom.emissivity * sigmaSB * (Twk^4 - T_inf^4);
    else
        qrad(k) = 0;
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

% =========================================================================
% 4. PLOTS (driver section) -- one figure per plot, no subplots
% =========================================================================
if opts.makePlots
    plotHeatFluxResults(results, geom);
end

end % ================== END MAIN FUNCTION ==================


% =========================================================================
% LOCAL FUNCTION: plot suite
% =========================================================================
function plotHeatFluxResults(res, geom)

t = res.time;

% ---- Trajectory: altitude vs time ----
figure('Name', 'Trajectory - Altitude vs Time');
plot(t, res.altitude/1000, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Altitude [km]');
title('Trajectory: Altitude vs. Time');
grid on;

% ---- Trajectory: velocity vs time ----
figure('Name', 'Trajectory - Velocity vs Time');
plot(t, res.velocity, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Trajectory: Velocity vs. Time');
grid on;

% ---- Mach number: freestream vs boundary-layer edge ----
figure('Name', 'Mach Number vs Time');
plot(t, res.Mach_inf, 'LineWidth', 1.6); hold on;
plot(t, res.Me, '--', 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Freestream vs. Boundary-Layer-Edge Mach Number vs. Time');
legend('Freestream, M_\infty', 'Boundary-layer edge, M_e', 'Location', 'best');
grid on; hold off;

% ---- Convective heat flux vs time ----
figure('Name', 'Convective Heat Flux vs Time');
plot(t, res.qw/1e4, 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Convective Heat Flux [W/cm^2]');
title(sprintf('Convective Heat Flux vs. Time (station x = %.3g m)', geom.x));
grid on;

% ---- Wall / recovery / edge / reference temperature ----
figure('Name', 'Wall, Recovery, Edge, and Reference Temperature vs Time');
plot(t, res.Tw,    'LineWidth', 1.9); hold on;
plot(t, res.Taw,   '--',  'LineWidth', 1.6);
plot(t, res.Te,    ':',   'LineWidth', 1.6);
plot(t, res.Tstar, '-.',  'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Temperature [K]');
title(sprintf('Wall Temperature vs. Adiabatic Wall, Edge, and Reference Temperature (x = %.3g m)', geom.x));
legend('Wall temperature, T_w', 'Adiabatic wall (recovery), T_{aw}', ...
       'Boundary-layer edge, T_e', 'Eckert reference temp, T^*', 'Location', 'best');
grid on; hold off;

% ---- Skin friction coefficient vs time ----
figure('Name', 'Skin Friction Coefficient vs Time');
plot(t, res.Cf, 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Skin Friction Coefficient, C_f [-]');
title(sprintf('Local Skin Friction Coefficient vs. Time (station x = %.3g m)', geom.x));
grid on;

end


% =========================================================================
% LOCAL FUNCTION: single-node wall energy balance (ODE right-hand side)
% =========================================================================
function dTwdt = wallEnergyBalance(t, Tw, hOf, VOf, alphaOf, geom, mat, opts, sigmaSB)

h     = hOf(t);
V     = VOf(t);
alpha = alphaOf(t);

[~, T_inf, p_inf, ~] = standardAtmosphere(h);
a_inf = sqrt(opts.gamma * opts.R * T_inf);
M_inf = V / a_inf;

if ~isempty(geom.halfAngle)
    thetaEff = geom.halfAngle + abs(alpha);
    [Me, Te, pe, ~] = obliqueShockEdge(M_inf, T_inf, p_inf, thetaEff, opts.gamma);
else
    Me = M_inf;
    Te = T_inf;
    pe = p_inf;
end

Ue = Me * sqrt(opts.gamma * opts.R * Te);

qw = referenceTempHeating(Me, Te, pe, Ue, Tw, geom.x, opts);

if opts.includeRadiation
    qrad = geom.emissivity * sigmaSB * (Tw^4 - T_inf^4);
else
    qrad = 0;
end

massArea = mat.density * mat.thickness; % kg/m^2 (thin-skin, per unit area)
dTwdt = (qw - qrad) / (massArea * mat.cp);

end


% =========================================================================
% LOCAL FUNCTION: Eckert reference temperature method
%   Returns local convective heat flux, skin friction coefficient,
%   adiabatic wall temperature, and reference temperature at a station
%   located a running length x from the leading edge/stagnation point.
% =========================================================================
function [qw, Cf, Taw, Tstar, regime] = referenceTempHeating(Me, Te, pe, Ue, Tw, x, opts)

gamma = opts.gamma;
R     = opts.R;
Pr    = opts.Pr;

% --- Eckert reference temperature (Anderson, Ref. Temp. Method) ---
% T*/Te = 1 + 0.032*Me^2 + 0.58*(Tw/Te - 1)
Tstar = Te * (1 + 0.032*Me^2 + 0.58*(Tw/Te - 1));
[rho_star, mu_star] = referenceGasProps(Tstar, pe, R);
Rex_star = rho_star * Ue * x / mu_star;

if Rex_star <= opts.Re_trans
    % ---- Laminar ----
    regime = 'laminar';
    r  = sqrt(Pr);                 % laminar recovery factor
    Cf = 0.664 / sqrt(Rex_star);   % Blasius, evaluated at T*
else
    % ---- Turbulent ----
    regime = 'turbulent';
    r  = Pr^(1/3);                 % turbulent recovery factor
    Cf = 0.0592 / Rex_star^0.2;    % Schlichting flat-plate, at T*
end

Taw = Te * (1 + r * (gamma - 1)/2 * Me^2);

% --- Reynolds-Colburn analogy for Stanton number at reference conditions ---
Cp_star = gamma * R / (gamma - 1); % calorically perfect gas
St_star = (Cf/2) / Pr^(2/3);

% --- Convective heat flux ---
qw = rho_star * Ue * Cp_star * St_star * (Taw - Tw);

end


% =========================================================================
% LOCAL FUNCTION: reference density & viscosity at T*
% =========================================================================
function [rho_star, mu_star] = referenceGasProps(Tstar, pe, R)

rho_star = pe / (R * Tstar); % perfect gas, static pressure ~const across BL

mu0 = 1.716e-5;   % kg/(m-s), Sutherland reference viscosity
T0  = 273.15;     % K
S   = 110.4;      % K
mu_star = mu0 * (Tstar/T0)^1.5 * (T0 + S) / (Tstar + S);

end


% =========================================================================
% LOCAL FUNCTION: exact oblique-shock boundary-layer-edge conditions
%   2-D wedge relation. The shock angle beta is the EXACT root of the
%   full theta-beta-M relation (solved via fzero) -- no small-angle or
%   weak-shock linearization is used.
% =========================================================================
function [M2, T2, p2, rho2] = obliqueShockEdge(M1, T1, p1, thetaDeg, gamma)

R_air = 287;
theta = deg2rad(thetaDeg);

if theta <= 0 || M1 <= 1
    M2 = M1; T2 = T1; p2 = p1; rho2 = p1/(R_air*T1);
    return;
end

muMach = asin(1/M1); % Mach angle, lower bound for beta

% Exact theta-beta-M relation, solved for the weak-shock angle beta:
%   tan(theta) = 2*cot(beta) * (M1^2 sin^2(beta) - 1) / (M1^2 (gamma+cos(2*beta)) + 2)
thetaBetaM = @(beta) tan(theta) - 2*cot(beta).* ...
    (M1^2 .* sin(beta).^2 - 1) ./ (M1^2 .* (gamma + cos(2*beta)) + 2);

betaLow  = muMach*1.0001;
betaHigh = pi/2 - 1e-6;

if thetaBetaM(betaLow) * thetaBetaM(betaHigh) > 0
    % Detached shock at this Mach/angle combination -- the wedge
    % relation has no attached-shock solution. Fall back to freestream
    % (a normal-shock/blunt-body model would be needed for a rigorous
    % treatment of this regime).
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
%   Returns geopotential height, temperature, pressure, density. Above
%   86 km, falls back to a simple exponential extrapolation (approximate
%   only).
% =========================================================================
function [Hgp, T, p, rho] = standardAtmosphere(h_geometric)

Re = 6356766;    % m, effective Earth radius for geopotential conversion
g0 = 9.80665;    % m/s^2
Rgas = 287.0528; % J/kg-K, air

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
% LOCAL FUNCTION: small helper to fill in default struct fields
% =========================================================================
function s = setDefault(s, field, value)
if ~isfield(s, field)
    s.(field) = value;
end
end
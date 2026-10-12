function results = hypersonicAeroheating(traj, geom, mat, opts)
%HYPERSONICAEROHEATING  Transient aeroheating analysis for a hypersonic
%projectile using Eckert's Reference Temperature Method.
%
%   results = HYPERSONICAEROHEATING(traj, geom, mat, opts)
%
%   Computes a time-resolved (transient) estimate of convective heat
%   flux, wall temperature, adiabatic (recovery) wall temperature,
%   boundary-layer edge temperature, and local skin-friction coefficient
%   along a flight trajectory, using the Reference Temperature Method
%   described in:
%
%       Anderson, J.D., "Hypersonic and High-Temperature Gas Dynamics,"
%       (relevant chapters on laminar/turbulent boundary layers and the
%       reference temperature method for compressible skin friction and
%       heat transfer).
%
%   The wall temperature is not assumed known a priori -- it is solved
%   as a state variable via a lumped-thermal-capacitance energy balance
%   (convective heating in, radiative loss out), integrated with ode45
%   across the supplied trajectory. This captures the transient thermal
%   response of the structure, not just the instantaneous heating rate.
%
%   ---------------------------------------------------------------
%   REQUIRED INPUT: traj (struct)
%   ---------------------------------------------------------------
%     traj.time      [s]     time vector, monotonically increasing
%     traj.altitude  [m]     altitude vector, same length as time
%     traj.velocity  [m/s]   velocity magnitude vector, same length
%     traj.alpha     [deg]   (optional) angle of attack vector, same
%                            length as time. Defaults to zero. Used only
%                            as a simple cosine modifier on the local
%                            edge Mach number for the windward ray; see
%                            NOTES below for the approximation made.
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: geom (struct) -- reference/placeholder values below
%   ---------------------------------------------------------------
%     geom.x            [m]   running length from stagnation point/nose
%                              to the surface station of interest, used
%                              as the reference length in the local
%                              Reynolds number. Default: 0.5 m
%     geom.Rn           [m]   nose radius, used only for the auxiliary
%                              stagnation-point heating estimate.
%                              Default: 0.03 m
%     geom.halfAngle    [deg] body half-angle (wedge-equivalent) used to
%                              estimate boundary-layer-edge conditions
%                              downstream of an attached oblique shock.
%                              Leave empty ([]) to skip the shock
%                              calculation and take edge conditions equal
%                              to freestream (flat-plate/slender-body
%                              approximation). Default: [] (freestream)
%     geom.emissivity   [-]   wall total hemispherical emissivity.
%                              Default: 0.80
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: mat (struct) -- lumped wall thermal properties
%   ---------------------------------------------------------------
%     mat.thickness   [m]        wall/skin thickness. Default: 0.003 m
%     mat.density     [kg/m^3]   wall material density. Default: 8000
%                                 (generic steel-like placeholder)
%     mat.cp          [J/kg-K]   wall material specific heat.
%                                 Default: 500
%     mat.Tw0         [K]        initial wall temperature. Default: 290
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: opts (struct) -- gas/model parameters
%   ---------------------------------------------------------------
%     opts.gamma          [-]  ratio of specific heats. Default: 1.4
%     opts.R              [J/kg-K] specific gas constant (air). Default: 287
%     opts.Pr             [-]  Prandtl number (assumed const). Default: 0.71
%     opts.Re_trans       [-]  transition Reynolds number (local, based
%                               on reference conditions) above which the
%                               boundary layer is treated as turbulent.
%                               Default: 5e5
%     opts.includeRadiation (logical) include radiative wall cooling.
%                               Default: true
%     opts.includeStagnation (logical) also compute an auxiliary
%                               stagnation-point convective heating
%                               estimate (Sutton-Graves correlation) for
%                               comparison. Default: true
%     opts.makePlots      (logical) generate the summary figure.
%                               Default: true
%
%   ---------------------------------------------------------------
%   OUTPUT: results (struct)
%   ---------------------------------------------------------------
%     results.time     [s]      time vector (solver output grid)
%     results.altitude [m]
%     results.velocity [m/s]
%     results.Mach_inf [-]      freestream Mach number
%     results.Me       [-]      boundary-layer edge Mach number
%     results.Te       [K]      boundary-layer edge static temperature
%     results.Tw       [K]      wall temperature (solved, transient)
%     results.Taw      [K]      adiabatic wall (recovery) temperature
%     results.Tstar    [K]      Eckert reference temperature
%     results.qw       [W/m^2]  convective heat flux to the wall
%     results.qrad     [W/m^2]  radiative heat flux leaving the wall
%     results.Cf       [-]      local skin-friction coefficient
%     results.regime   {'laminar'|'turbulent'} cell array per time step
%     results.qstag    [W/m^2]  (if enabled) stagnation-point convective
%                                heat flux, Sutton-Graves estimate
%
%   ---------------------------------------------------------------
%   NOTES / SIMPLIFICATIONS (fill in / refine as needed)
%   ---------------------------------------------------------------
%   1. Boundary-layer edge conditions: if geom.halfAngle is supplied, a
%      2-D oblique-shock (wedge) relation is used to estimate edge
%      conditions. This is only an approximation for an axisymmetric
%      body (a true conical analysis requires the Taylor-Maccoll
%      equations); it is adequate for a first-pass reference-temperature
%      heating estimate and can be replaced with a better edge-condition
%      model later without changing the rest of the pipeline.
%   2. Angle of attack is folded in only as a crude cosine-type modifier
%      on the effective flow deflection seen by the windward ray
%      (theta_eff = geom.halfAngle + alpha). This is a placeholder --
%      replace with a proper windward/leeward ray analysis if AoA
%      effects matter for your case.
%   3. The atmosphere model is a standard 1976 U.S. Standard Atmosphere
%      implementation valid to 86 km geometric altitude; above that it
%      falls back to an exponential extrapolation and should be treated
%      as approximate only.
%   4. The lumped-capacitance wall model assumes a thermally-thin skin
%      with no through-thickness gradient and no conduction to
%      surrounding structure. Replace mat.* with actual vehicle values
%      when available; consider a finite-difference through-thickness
%      model if the thermally-thin assumption is not valid (a large
%      Biot number).
%
%   Example:
%       traj.time     = linspace(0,60,300);
%       traj.altitude = linspace(40000,2000,300);   % m
%       traj.velocity = linspace(3000,1200,300);    % m/s
%       results = hypersonicAeroheating(traj, [], [], []);

% =========================================================================
% 1. INPUT HANDLING / DEFAULTS
% =========================================================================
if nargin < 2, geom = struct(); end
if nargin < 3, mat  = struct(); end
if nargin < 4, opts = struct(); end

t_traj = traj.time(:);
h_traj = traj.altitude(:);
V_traj = traj.velocity(:);
if isfield(traj,'alpha') && ~isempty(traj.alpha)
    alpha_traj = traj.alpha(:);
else
    alpha_traj = zeros(size(t_traj));
end

if ~isequal(numel(t_traj), numel(h_traj), numel(V_traj), numel(alpha_traj))
    error('hypersonicAeroheating:sizeMismatch', ...
        'traj.time, traj.altitude, traj.velocity, and traj.alpha must be the same length.');
end
if any(diff(t_traj) <= 0)
    error('hypersonicAeroheating:timeNotMonotonic', ...
        'traj.time must be strictly increasing.');
end

geom = setDefault(geom, 'x',           0.5);     % m, reference station
geom = setDefault(geom, 'Rn',          0.03);    % m, nose radius (placeholder)
geom = setDefault(geom, 'halfAngle',   []);       % deg, [] = freestream edge
geom = setDefault(geom, 'emissivity',  0.80);

mat  = setDefault(mat, 'thickness',    0.003);   % m
mat  = setDefault(mat, 'density',      8000);    % kg/m^3
mat  = setDefault(mat, 'cp',           500);     % J/kg-K
mat  = setDefault(mat, 'Tw0',          290);     % K

opts = setDefault(opts, 'gamma',              1.4);
opts = setDefault(opts, 'R',                  287);
opts = setDefault(opts, 'Pr',                 0.71);
opts = setDefault(opts, 'Re_trans',           5e5);
opts = setDefault(opts, 'includeRadiation',   true);
opts = setDefault(opts, 'includeStagnation',  true);
opts = setDefault(opts, 'makePlots',          true);

sigmaSB = 5.670374419e-8; % Stefan-Boltzmann constant, W/m^2-K^4

% Interpolants used inside the ODE integration (queried at arbitrary t)
hOf     = @(tq) interp1(t_traj, h_traj,     tq, 'linear', 'extrap');
VOf     = @(tq) interp1(t_traj, V_traj,     tq, 'linear', 'extrap');
alphaOf = @(tq) interp1(t_traj, alpha_traj, tq, 'linear', 'extrap');

% =========================================================================
% 2. TRANSIENT SOLVE: wall temperature via lumped energy balance
% =========================================================================
odefun = @(t, Tw) wallEnergyBalance(t, Tw, hOf, VOf, alphaOf, geom, mat, ...
                                     opts, sigmaSB);

odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-3);
[tSol, TwSol] = ode45(odefun, t_traj, mat.Tw0, odeOptions);

% =========================================================================
% 3. POST-PROCESS: recompute all quantities of interest on the solution grid
% =========================================================================
N = numel(tSol);
Mach_inf = zeros(N,1);
Me       = zeros(N,1);
Te       = zeros(N,1);
Taw      = zeros(N,1);
Tstar    = zeros(N,1);
qw       = zeros(N,1);
qrad     = zeros(N,1);
Cf       = zeros(N,1);
regime   = cell(N,1);
qstag    = zeros(N,1);
altOut   = zeros(N,1);
velOut   = zeros(N,1);

for k = 1:N
    tk  = tSol(k);
    Twk = TwSol(k);

    hk     = hOf(tk);
    Vk     = VOf(tk);
    alphak = alphaOf(tk);

    altOut(k) = hk;
    velOut(k) = Vk;

    [rho_inf, p_inf, T_inf] = standardAtmosphere(hk);
    a_inf = sqrt(opts.gamma * opts.R * T_inf);
    M_inf = Vk / a_inf;
    Mach_inf(k) = M_inf;

    % --- Boundary-layer edge conditions ---
    if ~isempty(geom.halfAngle)
        thetaEff = geom.halfAngle + abs(alphak); % crude AoA placeholder, deg
        [Me_k, Te_k, pe_k, ~] = obliqueShockEdge(M_inf, T_inf, p_inf, ...
                                                  thetaEff, opts.gamma);
    else
        Me_k = M_inf;
        Te_k = T_inf;
        pe_k = p_inf;
    end
    Me(k) = Me_k;
    Te(k) = Te_k;

    Ue_k = Me_k * sqrt(opts.gamma * opts.R * Te_k);

    % --- Recovery / reference temperature / local Cf, qw ---
    [qw_k, Cf_k, Taw_k, Tstar_k, regime_k] = referenceTempHeating( ...
        Me_k, Te_k, pe_k, Ue_k, Twk, geom.x, opts);

    Taw(k)    = Taw_k;
    Tstar(k)  = Tstar_k;
    qw(k)     = qw_k;
    Cf(k)     = Cf_k;
    regime{k} = regime_k;

    if opts.includeRadiation
        qrad(k) = geom.emissivity * sigmaSB * (Twk^4 - T_inf^4);
    else
        qrad(k) = 0;
    end

    if opts.includeStagnation
        % Sutton-Graves correlation, SI units (approximate engineering
        % estimate of convective stagnation-point heating, independent
        % of the reference-temperature flat-plate result above; useful
        % as a sanity check / comparison for the nose region).
        C_sg = 1.7415e-4;
        qstag(k) = C_sg * sqrt(rho_inf / geom.Rn) * Vk^3;
    else
        qstag(k) = NaN;
    end
end

results.time     = tSol;
results.altitude = altOut;
results.velocity = velOut;
results.Mach_inf = Mach_inf;
results.Me       = Me;
results.Te       = Te;
results.Tw       = TwSol;
results.Taw      = Taw;
results.Tstar    = Tstar;
results.qw       = qw;
results.qrad     = qrad;
results.Cf       = Cf;
results.regime   = regime;
if opts.includeStagnation
    results.qstag = qstag;
end

% =========================================================================
% 4. PLOTS
% =========================================================================
if opts.makePlots
    figure('Name', 'Hypersonic Aeroheating - Reference Temperature Method', ...
           'Color', 'w');

    subplot(3,1,1);
    plot(results.time, results.qw/1e3, 'LineWidth', 1.6); hold on;
    if opts.includeStagnation
        plot(results.time, results.qstag/1e3, '--', 'LineWidth', 1.2);
        legend('Local (ref. temp method)', 'Stagnation point (Sutton-Graves)', ...
               'Location', 'best');
    end
    xlabel('Time [s]');
    ylabel('Heat flux [kW/m^2]');
    title('Convective Heat Flux vs Time');
    grid on;

    subplot(3,1,2);
    plot(results.time, results.Tw,  'LineWidth', 1.6); hold on;
    plot(results.time, results.Taw, '--', 'LineWidth', 1.4);
    plot(results.time, results.Te,  ':', 'LineWidth', 1.4);
    xlabel('Time [s]');
    ylabel('Temperature [K]');
    title('Wall / Recovery / Boundary-Layer-Edge Temperature vs Time');
    legend('Wall temperature T_w', 'Adiabatic wall (recovery) T_{aw}', ...
           'Boundary-layer edge T_e', 'Location', 'best');
    grid on;

    subplot(3,1,3);
    plot(results.time, results.Cf, 'LineWidth', 1.6);
    xlabel('Time [s]');
    ylabel('Skin friction coefficient C_f [-]');
    title('Local Skin Friction Coefficient vs Time');
    grid on;
end

end % ================== END MAIN FUNCTION ==================


% =========================================================================
% LOCAL FUNCTION: wall energy balance (RHS of the transient ODE)
% =========================================================================
function dTwdt = wallEnergyBalance(t, Tw, hOf, VOf, alphaOf, geom, mat, ...
                                    opts, sigmaSB)

h     = hOf(t);
V     = VOf(t);
alpha = alphaOf(t);

[~, p_inf, T_inf] = standardAtmosphere(h);
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

% Lumped thermal capacitance: rho*t*cp*dTw/dt = q_convective_in - q_radiative_out
massPerArea = mat.density * mat.thickness;
dTwdt = (qw - qrad) / (massPerArea * mat.cp);

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
% Note: this does NOT depend on the recovery factor, so it can be
% evaluated once, before the laminar/turbulent regime is even known.
Tstar = Te * (1 + 0.032*Me^2 + 0.58*(Tw/Te - 1));
[rho_star, mu_star] = referenceGasProps(Tstar, pe, R);
Rex_star = rho_star * Ue * x / mu_star;

% --- Regime check (based on reference-condition local Reynolds number) ---
if Rex_star <= opts.Re_trans
    % ---- Laminar ----
    regime = 'laminar';
    r  = sqrt(Pr);                            % laminar recovery factor
    Cf = 0.664 / sqrt(Rex_star);              % Blasius, evaluated at T*
else
    % ---- Turbulent ----
    regime = 'turbulent';
    r  = Pr^(1/3);                            % turbulent recovery factor
    Cf = 0.0592 / Rex_star^0.2;               % Schlichting flat-plate, at T*
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

% Static pressure is ~constant across the boundary layer -> perfect gas
rho_star = pe / (R * Tstar);

% Sutherland's law for air
mu0 = 1.716e-5;   % kg/(m-s) at T0
T0  = 273.15;     % K
S   = 110.4;      % K
mu_star = mu0 * (Tstar/T0)^1.5 * (T0 + S) / (Tstar + S);

end


% =========================================================================
% LOCAL FUNCTION: approximate oblique-shock boundary-layer-edge conditions
%   2-D wedge relation, weak-shock root. Approximation for axisymmetric
%   bodies -- replace with Taylor-Maccoll cone solution for higher
%   fidelity if needed.
% =========================================================================
function [M2, T2, p2, rho2] = obliqueShockEdge(M1, T1, p1, thetaDeg, gamma)

theta = deg2rad(thetaDeg);

if theta <= 0 || M1 <= 1
    M2 = M1; T2 = T1; p2 = p1; rho2 = p1/(287*T1);
    return;
end

muMach = asin(1/M1); % Mach angle, lower bound for beta

% theta-beta-M relation, solved for weak shock angle beta
thetaBetaM = @(beta) tan(theta) - 2*cot(beta).* ...
    (M1^2 .* sin(beta).^2 - 1) ./ (M1^2 .* (gamma + cos(2*beta)) + 2);

betaLow  = muMach*1.0001;
betaHigh = pi/2 - 1e-6;

% Guard against detached shock (no real solution for the given theta)
if thetaBetaM(betaLow) * thetaBetaM(betaHigh) > 0
    % Detached / too blunt for the wedge relation at this Mach number;
    % fall back to freestream (a normal-shock/blunt-body model should
    % replace this branch for a more rigorous high-alpha treatment).
    M2 = M1; T2 = T1; p2 = p1; rho2 = p1/(287*T1);
    return;
end

beta = fzero(thetaBetaM, [betaLow, betaHigh]);

M1n = M1 * sin(beta);
p2_p1   = 1 + 2*gamma/(gamma+1) * (M1n^2 - 1);
rho2_rho1 = (gamma+1) * M1n^2 / ((gamma-1) * M1n^2 + 2);
T2_T1   = p2_p1 / rho2_rho1;

M2n = sqrt( (1 + (gamma-1)/2 * M1n^2) / (gamma*M1n^2 - (gamma-1)/2) );
M2  = M2n / sin(beta - theta);

p2   = p1 * p2_p1;
T2   = T1 * T2_T1;
rho2 = rho2_rho1 * (p1/(287*T1));

end


function [rho,p,T] = standardAtmosphere(h)
% PURPOSE: Computes density (rho [kg/m^3]) at a specified altitude (h [m])
%          Also returns pressure (Pa) and temperature (K)

g  = 9.8;  % gravitational acceleration [m/s^2]
R  = 287;  % gas constant for air [J/kg.K]

hl = 1000*[0, 11, 20, 32, 47, 51, 71, 85, 90, 140, 500];  % altitudes at layer endpoints [m]
Tl = [288, 217, 217, 229, 271, 271, 215, 187, 187, 320, 700];
Ll = (Tl(2:end)-Tl(1:end-1))./(hl(2:end)-hl(1:end-1)); % lapse rate for each layer [K/m]
pl = [101325, 2.26656e+04, 5.49942e+03, 8.75227e+02, 1.12265e+02, ...
      6.78201e+01, 4.03063e+00, 3.72214e-01, 1.49376e-01, 1.51115e-04];% pressures at endpoints [Pa]

% identify the layer in which h is found
h = max(h,0); I = find(hl-h>0);
if (isempty(I)), rho = 0.; p=0; T=288; return; end;
l = I(1)-1; T = Tl(l) + Ll(l)*(h-hl(l)); % temperature at the point
p0 = pl(l);                  % pressure at the start of the layer

% pressure, density at the point
if (Ll(l) == 0), % isothermal layer
  p = p0*exp(-g*(h-hl(l))./(R*T));
else             % constant-temperature-gradient layer
  p = p0*(T/Tl(l)).^(-g/(R*Ll(l)));
end
rho = p/(R*T);
end

% =========================================================================
% LOCAL FUNCTION: small helper to fill in default struct fields
% =========================================================================
function s = setDefault(s, field, value)
if ~isfield(s, field) || isempty(s.(field))
    if isfield(s, field) && isempty(s.(field)) && ...
            (strcmp(field,'halfAngle')) % allow explicit [] to mean "skip"
        return;
    end
    s.(field) = value;
end
end

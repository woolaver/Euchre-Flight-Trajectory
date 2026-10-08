function results = hypersonicAeroheating(traj, geom, mat, opts)
%HYPERSONICAEROHEATING  Transient aeroheating analysis for a hypersonic
%projectile using Eckert's Reference Temperature Method, with a
%through-thickness (outer/inner) wall thermal model and a spatial survey
%of heating along the body length.
%
%   results = HYPERSONICAEROHEATING(traj, geom, mat, opts)
%
%   Computes a time-resolved (transient) estimate of convective heat
%   flux, OUTER (aero-heated) and INNER (payload-facing) wall
%   temperature, adiabatic (recovery) wall temperature, boundary-layer
%   edge temperature, and local skin-friction coefficient along a flight
%   trajectory, using the Reference Temperature Method described in:
%
%       Anderson, J.D., "Hypersonic and High-Temperature Gas Dynamics,"
%       (relevant chapters on laminar/turbulent boundary layers and the
%       reference temperature method for compressible skin friction and
%       heat transfer).
%
%   The calculation is performed simultaneously at an array of body
%   stations (geom.xStations) running from near the nose to the aft
%   body, so heating and wall temperature can be examined both as a
%   function of TIME (at a representative "primary" station) and as a
%   function of POSITION along the body (at selected snapshots in time).
%
%   WALL THERMAL MODEL (per station)
%   ---------------------------------------------------------------
%   Each station is modeled as a 2-node through-thickness lumped-
%   capacitance system:
%     - Outer node (T_outer): the aero-heated surface. Receives
%       convective heating from the boundary layer, loses heat by
%       radiation to the freestream, and conducts heat inward.
%     - Inner node (T_inner): the payload-facing surface. Receives
%       conduction from the outer node and (optionally) loses heat to an
%       internal payload cavity by convection.
%   Each node is treated as a lumped half-thickness slice of the wall;
%   conduction between the two node centers uses the material's thermal
%   conductivity mat.k. This is the standard 2-node finite-difference
%   idealization of a 1-D slab and is the minimum level of fidelity
%   needed to separate "outer skin" temperature from "payload-side"
%   temperature. If very large through-thickness gradients or very thick
%   walls are of interest, this can be extended to more nodes without
%   changing the rest of the pipeline (see wallEnergyBalanceMulti).
%
%   The wall temperatures are NOT assumed known a priori -- they are
%   solved as state variables via the coupled energy balances above,
%   integrated with ode45 across the supplied trajectory. This captures
%   the transient thermal response of the structure, not just the
%   instantaneous heating rate.
%
%   ---------------------------------------------------------------
%   REQUIRED INPUT: traj (struct)
%   ---------------------------------------------------------------
%     traj.time      [s]     time vector, monotonically increasing
%     traj.altitude  [m]     altitude vector, same length as time
%     traj.velocity  [m/s]   velocity magnitude vector, same length
%     traj.alpha     [deg]   (optional) angle of attack vector, same
%                            length as time. Defaults to zero. Used only
%                            as a simple cosine-type modifier on the
%                            effective flow deflection seen by the
%                            windward ray; see NOTES below.
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: geom (struct) -- reference/placeholder values below
%   ---------------------------------------------------------------
%     geom.xStations   [m]   vector of running-length stations from the
%                              nose/stagnation point to survey along the
%                              body. Default: linspace(0.05, 1.0, 15)
%     geom.xPrimary    [m]   the station (nearest entry in xStations)
%                              used for all "vs. time" plots and for the
%                              primary-station output fields (Touter,
%                              Tinner, qw, Cf, etc.). Default: 0.5 m
%     geom.Rn          [m]   nose radius, used only for the auxiliary
%                              stagnation-point heating estimate.
%                              Default: 0.03 m
%     geom.halfAngle   [deg] body half-angle (wedge-equivalent) used to
%                              estimate boundary-layer-edge conditions
%                              downstream of an attached oblique shock.
%                              Leave empty ([]) to skip the shock
%                              calculation and take edge conditions equal
%                              to freestream (flat-plate/slender-body
%                              approximation). Default: [] (freestream)
%     geom.emissivity  [-]   outer wall total hemispherical emissivity.
%                              Default: 0.80
%
%   ---------------------------------------------------------------
%   OPTIONAL INPUT: mat (struct) -- lumped wall thermal properties
%   ---------------------------------------------------------------
%     mat.thickness    [m]        total wall/skin thickness (outer node
%                                   and inner node each represent one
%                                   half of this). Default: 0.003 m
%     mat.density      [kg/m^3]   wall material density. Default: 8000
%                                   (generic steel-like placeholder)
%     mat.cp           [J/kg-K]   wall material specific heat.
%                                   Default: 500
%     mat.k            [W/m-K]    wall material thermal conductivity,
%                                   used for outer-to-inner conduction.
%                                   Default: 16 (generic steel-like)
%     mat.h_internal   [W/m^2-K]  convective heat transfer coefficient
%                                   from the inner wall surface to the
%                                   payload cavity air. Default: 0
%                                   (adiabatic inner surface -- the most
%                                   conservative/worst-case assumption
%                                   for payload thermal margin; increase
%                                   to model natural/forced convection
%                                   inside the payload bay)
%     mat.T_payload    [K]        effective payload cavity air
%                                   temperature used with h_internal.
%                                   Default: 290
%     mat.Tw0          [K]        initial temperature of BOTH the outer
%                                   and inner nodes (uniform IC).
%                                   Default: 290
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
%     opts.nSnapshots     (integer) number of time snapshots used for the
%                               "vs. body station" family-of-curves plots.
%                               Default: 6
%     opts.snapshotSpacing ('log'|'linear') how the snapshot times are
%                               distributed across the trajectory.
%                               'log' (default) clusters snapshots densely
%                               just after the first time point and
%                               spreads them out logarithmically toward
%                               the end -- useful for resolving the sharp
%                               thermal transient right after launch,
%                               where Touter/Tinner/qw change fastest.
%                               'linear' spaces them evenly in time
%                               (the previous default behavior).
%     opts.snapshotLogFrac [-]  only used when snapshotSpacing = 'log'.
%                               Fraction of the total trajectory duration
%                               used to place the *second* snapshot (the
%                               first is always t = traj.time(1)). E.g.
%                               0.001 places the second snapshot 0.1% of
%                               the way through the flight, so with a
%                               60-second flight that is 0.06 s after
%                               launch. Smaller values push the earliest
%                               post-launch snapshots even closer to t0.
%                               Default: 1e-3
%     opts.makePlots      (logical) generate the full driver plot suite.
%                               Default: true
%
%   ---------------------------------------------------------------
%   OUTPUT: results (struct)
%   ---------------------------------------------------------------
%   Primary-station (geom.xPrimary) time-history fields, each [N x 1]:
%     results.time      [s]
%     results.altitude  [m]
%     results.velocity  [m/s]
%     results.Mach_inf  [-]      freestream Mach number
%     results.Me        [-]      boundary-layer edge Mach number
%     results.Te        [K]      boundary-layer edge static temperature
%     results.Touter     [K]      OUTER (aero-heated) wall temperature
%     results.Tinner     [K]      INNER (payload-facing) wall temperature
%     results.Tw         [K]      alias of results.Touter (back-compat)
%     results.Taw        [K]      adiabatic wall (recovery) temperature
%     results.Tstar      [K]      Eckert reference temperature
%     results.qw         [W/m^2]  convective heat flux to the outer wall
%     results.qrad       [W/m^2]  radiative heat flux leaving the outer wall
%     results.qcond      [W/m^2]  conduction flux from outer to inner node
%     results.Cf         [-]      local skin-friction coefficient
%     results.regime     {'laminar'|'turbulent'} cell array per time step
%     results.qstag      [W/m^2]  (if enabled) stagnation-point convective
%                                  heat flux, Sutton-Graves estimate
%
%   Full spatio-temporal survey fields, each [N x Nx] (rows = time,
%   columns = body station):
%     results.xStations   [m]      the surveyed body stations (1 x Nx)
%     results.primaryIdx  [-]      column index nearest geom.xPrimary
%     results.TouterField [K]
%     results.TinnerField [K]
%     results.qwField     [W/m^2]
%     results.CfField      [-]
%     results.TawField     [K]
%     results.TstarField   [K]
%     results.isTurbField  [logical] true where turbulent
%     results.snapshotIdx  [-]      time indices used for the "vs x" plots
%
%   ---------------------------------------------------------------
%   NOTES / SIMPLIFICATIONS (fill in / refine as needed)
%   ---------------------------------------------------------------
%   1. Boundary-layer edge conditions: if geom.halfAngle is supplied, a
%      2-D oblique-shock (wedge) relation is used to estimate edge
%      conditions, applied uniformly along the whole body (only the local
%      Reynolds number, and hence Cf/heat flux, varies with x). This is
%      only an approximation for an axisymmetric body (a true conical
%      analysis requires the Taylor-Maccoll equations) and does not
%      capture expansion/compression corners; adequate for a first-pass
%      estimate and can be refined later without changing the rest of
%      the pipeline.
%   2. Angle of attack is folded in only as a crude cosine-type modifier
%      on the effective flow deflection seen by the windward ray
%      (theta_eff = geom.halfAngle + alpha). This is a placeholder --
%      replace with a proper windward/leeward ray analysis if AoA
%      effects matter for your case.
%   3. The atmosphere model is a standard 1976 U.S. Standard Atmosphere
%      implementation valid to 86 km geometric altitude; above that it
%      falls back to an exponential extrapolation and should be treated
%      as approximate only.
%   4. The 2-node through-thickness wall model assumes 1-D conduction
%      only (no in-plane/streamwise conduction between stations, no
%      circumferential conduction) and lumps each half-thickness as an
%      isothermal node. Replace mat.* with actual vehicle values when
%      available, and add more nodes through the thickness if the
%      through-thickness gradient needs finer resolution.
%   5. The default mat.h_internal = 0 (adiabatic inner wall) is the
%      conservative bounding case for payload temperature margin -- it
%      gives the highest possible inner-wall temperature for a given
%      outer heating history. Set it to a representative internal
%      convection coefficient once the payload bay is defined.
%
%   Example:
%       traj.time     = linspace(0,60,300);
%       traj.altitude = linspace(40000,2000,300);   % m
%       traj.velocity = linspace(3000,1200,300);    % m/s
%       results = hypersonicAeroheating(traj, [], [], []);

% =========================================================================
% 1. INPUT HANDLING / DEFAULTS
% =========================================================================
if nargin < 2 || isempty(geom), geom = struct(); end
if nargin < 3 || isempty(mat),  mat  = struct(); end
if nargin < 4 || isempty(opts), opts = struct(); end

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

geom = setDefault(geom, 'xStations',   linspace(0.05, 1.0, 15));
geom = setDefault(geom, 'xPrimary',    0.5);
geom = setDefault(geom, 'Rn',          0.03);    % m, nose radius (placeholder)
geom = setDefault(geom, 'halfAngle',   []);       % deg, [] = freestream edge
geom = setDefault(geom, 'emissivity',  0.80);

mat  = setDefault(mat, 'thickness',    0.003);   % m
mat  = setDefault(mat, 'density',      8000);    % kg/m^3
mat  = setDefault(mat, 'cp',           500);     % J/kg-K
mat  = setDefault(mat, 'k',            16);      % W/m-K
mat  = setDefault(mat, 'h_internal',   0);       % W/m^2-K (0 = adiabatic inner wall)
mat  = setDefault(mat, 'T_payload',    290);     % K
mat  = setDefault(mat, 'Tw0',          290);     % K, initial temp, both nodes

opts = setDefault(opts, 'gamma',              1.4);
opts = setDefault(opts, 'R',                  287);
opts = setDefault(opts, 'Pr',                 0.71);
opts = setDefault(opts, 'Re_trans',           5e5);
opts = setDefault(opts, 'includeRadiation',   true);
opts = setDefault(opts, 'includeStagnation',  true);
opts = setDefault(opts, 'nSnapshots',         6);
opts = setDefault(opts, 'snapshotSpacing',    'log');
opts = setDefault(opts, 'snapshotLogFrac',    1e-3);
opts = setDefault(opts, 'makePlots',          true);

sigmaSB = 5.670374419e-8; % Stefan-Boltzmann constant, W/m^2-K^4

xStations = geom.xStations(:)';  % row vector, 1 x Nx
Nx = numel(xStations);
[~, primaryIdx] = min(abs(xStations - geom.xPrimary));

% Interpolants used inside the ODE integration (queried at arbitrary t)
hOf     = @(tq) interp1(t_traj, h_traj,     tq, 'linear', 'extrap');
VOf     = @(tq) interp1(t_traj, V_traj,     tq, 'linear', 'extrap');
alphaOf = @(tq) interp1(t_traj, alpha_traj, tq, 'linear', 'extrap');

% =========================================================================
% 2. TRANSIENT SOLVE: outer/inner wall temperature at every station
% =========================================================================
% State vector layout: y = [Touter(1..Nx); Tinner(1..Nx)]
y0 = [repmat(mat.Tw0, Nx, 1); repmat(mat.Tw0, Nx, 1)];

odefun = @(t, y) wallEnergyBalanceMulti(t, y, hOf, VOf, alphaOf, xStations, ...
                                         Nx, geom, mat, opts, sigmaSB);

odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-3);
[tSol, Ysol] = ode45(odefun, t_traj, y0, odeOptions);

% =========================================================================
% 3. POST-PROCESS: recompute all quantities of interest on the solution grid
% =========================================================================
N = numel(tSol);

Mach_inf = zeros(N,1);
Me       = zeros(N,1);
Te       = zeros(N,1);
altOut   = zeros(N,1);
velOut   = zeros(N,1);
qstag    = zeros(N,1);

TouterField = zeros(N, Nx);
TinnerField = zeros(N, Nx);
qwField     = zeros(N, Nx);
CfField     = zeros(N, Nx);
TawField    = zeros(N, Nx);
TstarField  = zeros(N, Nx);
isTurbField = false(N, Nx);
qradField   = zeros(N, Nx);
qcondField  = zeros(N, Nx);

for k = 1:N
    tk = tSol(k);
    Touter_k = Ysol(k, 1:Nx);
    Tinner_k = Ysol(k, Nx+1:2*Nx);

    hk     = hOf(tk);
    Vk     = VOf(tk);
    alphak = alphaOf(tk);

    altOut(k) = hk;
    velOut(k) = Vk;

    [~, T_inf, p_inf, rho_inf] = standardAtmosphere(hk);
    a_inf = sqrt(opts.gamma * opts.R * T_inf);
    M_inf = Vk / a_inf;
    Mach_inf(k) = M_inf;

    % --- Boundary-layer edge conditions (uniform along body) ---
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

    for i = 1:Nx
        [qw_i, Cf_i, Taw_i, Tstar_i, regime_i] = referenceTempHeating( ...
            Me_k, Te_k, pe_k, Ue_k, Touter_k(i), xStations(i), opts);

        TouterField(k,i) = Touter_k(i);
        TinnerField(k,i) = Tinner_k(i);
        qwField(k,i)     = qw_i;
        CfField(k,i)     = Cf_i;
        TawField(k,i)    = Taw_i;
        TstarField(k,i)  = Tstar_i;
        isTurbField(k,i) = strcmp(regime_i, 'turbulent');

        if opts.includeRadiation
            qradField(k,i) = geom.emissivity * sigmaSB * (Touter_k(i)^4 - T_inf^4);
        else
            qradField(k,i) = 0;
        end
        qcondField(k,i) = 2*mat.k*(Touter_k(i) - Tinner_k(i)) / mat.thickness;
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

regimeCell = cell(N,1);
regimeCell(isTurbField(:,primaryIdx))  = {'turbulent'};
regimeCell(~isTurbField(:,primaryIdx)) = {'laminar'};

results.time     = tSol;
results.altitude = altOut;
results.velocity = velOut;
results.Mach_inf = Mach_inf;
results.Me       = Me;
results.Te       = Te;

results.Touter   = TouterField(:, primaryIdx);
results.Tinner   = TinnerField(:, primaryIdx);
results.Tw       = results.Touter;  % back-compat alias
results.Taw      = TawField(:, primaryIdx);
results.Tstar    = TstarField(:, primaryIdx);
results.qw       = qwField(:, primaryIdx);
results.qrad     = qradField(:, primaryIdx);
results.qcond    = qcondField(:, primaryIdx);
results.Cf       = CfField(:, primaryIdx);
results.regime   = regimeCell;
if opts.includeStagnation
    results.qstag = qstag;
end

results.xStations   = xStations;
results.primaryIdx  = primaryIdx;
results.TouterField = TouterField;
results.TinnerField = TinnerField;
results.qwField     = qwField;
results.CfField     = CfField;
results.TawField    = TawField;
results.TstarField  = TstarField;
results.isTurbField = isTurbField;

nSnap = max(2, min(opts.nSnapshots, N));
t0      = tSol(1);
tEnd    = tSol(end);
totalDur = tEnd - t0;

if strcmpi(opts.snapshotSpacing, 'log') && totalDur > 0
    % First snapshot is always the initial time point (t0). The
    % remaining (nSnap-1) snapshots are log-spaced across the rest of
    % the flight, so they cluster tightly just after launch -- where
    % Touter/Tinner/qw change fastest -- and spread out at later times.
    tFirstOffset = max(opts.snapshotLogFrac * totalDur, eps(totalDur)*10);
    snapshotTimes = [t0, t0 + logspace(log10(tFirstOffset), log10(totalDur), nSnap-1)];
else
    % Evenly spaced in time (previous default behavior).
    snapshotTimes = linspace(t0, tEnd, nSnap);
end

snapshotIdx = zeros(1, numel(snapshotTimes));
for s = 1:numel(snapshotTimes)
    [~, snapshotIdx(s)] = min(abs(tSol - snapshotTimes(s)));
end
snapshotIdx = unique(sort(snapshotIdx), 'stable');
results.snapshotIdx = snapshotIdx;

% =========================================================================
% 4. PLOTS (driver section)
% =========================================================================
% One figure per plot -- no subplots -- so each result can be inspected,
% zoomed, and saved independently during troubleshooting/verification.
% Multiple related lines are still overlaid on a single plot where that
% aids comparison (e.g., Touter/Taw/Te/T* on one plot, or several time
% snapshots on one "vs. body station" plot).
if opts.makePlots
    plotAeroheatingDriver(results, geom, mat, opts);
end

end % ================== END MAIN FUNCTION ==================


% =========================================================================
% LOCAL FUNCTION: full driver plotting suite
%   Separate figure per plot (no subplots). Covers the trajectory itself,
%   the primary-station time histories, and the along-body spatial survey
%   at several flight-time snapshots.
% =========================================================================
function plotAeroheatingDriver(res, geom, mat, opts)

t = res.time;
x = res.xStations;
snapIdx = res.snapshotIdx;
nSnap = numel(snapIdx);
snapColors = lines(nSnap);

% ---- Trajectory: altitude vs time ----
figure('Name', 'Trajectory - Altitude vs Time', 'Color', 'w');
plot(t, res.altitude/1000, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Altitude [km]');
title('Trajectory: Altitude vs. Time');
grid on;

% ---- Trajectory: velocity vs time ----
figure('Name', 'Trajectory - Velocity vs Time', 'Color', 'w');
plot(t, res.velocity, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Trajectory: Velocity vs. Time');
grid on;

% ---- Trajectory map: altitude vs velocity ----
figure('Name', 'Trajectory - Altitude vs Velocity', 'Color', 'w');
plot(res.velocity, res.altitude/1000, 'LineWidth', 1.6);
xlabel('Velocity [m/s]'); ylabel('Altitude [km]');
title('Trajectory Map: Altitude vs. Velocity');
grid on;

% ---- Freestream Mach number vs time ----
figure('Name', 'Freestream Mach Number vs Time', 'Color', 'w');
plot(t, res.Mach_inf, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Freestream Mach Number [-]');
title('Freestream Mach Number vs. Time');
grid on;

% ---- Boundary-layer edge Mach number vs time (freestream overlay) ----
figure('Name', 'Boundary-Layer Edge vs Freestream Mach Number', 'Color', 'w');
plot(t, res.Mach_inf, 'LineWidth', 1.6); hold on;
plot(t, res.Me, '--', 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Mach Number [-]');
title('Boundary-Layer Edge Mach Number vs. Freestream Mach Number');
legend('Freestream, M_\infty', 'Boundary-layer edge, M_e', 'Location', 'best');
grid on; hold off;

% ---- Boundary-layer edge temperature vs time ----
figure('Name', 'Boundary-Layer Edge Temperature vs Time', 'Color', 'w');
plot(t, res.Te, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Boundary-Layer Edge Temperature, T_e [K]');
title('Boundary-Layer Edge Static Temperature vs. Time');
grid on;

% ---- Heat flux: local (reference-temperature method) vs stagnation ----
figure('Name', 'Convective Heat Flux vs Time', 'Color', 'w');
plot(t, res.qw/1e4, 'LineWidth', 1.8); hold on;
if isfield(res, 'qstag')
    plot(t, res.qstag/1e4, '--', 'LineWidth', 1.6);
    legend(sprintf('Primary station, x = %.3g m (ref. temp. method)', geom.xPrimary), ...
           'Stagnation point (Sutton-Graves)', 'Location', 'best');
end
xlabel('Time [s]'); ylabel('Convective Heat Flux [W/cm^2]');
title('Convective Heat Flux vs. Time');
grid on; hold off;

% ---- Radiative heat flux leaving the outer wall vs time ----
figure('Name', 'Radiative Heat Flux vs Time', 'Color', 'w');
plot(t, res.qrad/1e4, 'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Radiative Heat Flux [W/cm^2]');
title('Radiative Heat Flux Leaving the Outer Wall vs. Time');
grid on;

% ---- Net heat flux into the outer wall (convective minus radiative) ----
figure('Name', 'Net Heat Flux Into Outer Wall vs Time', 'Color', 'w');
plot(t, (res.qw - res.qrad)/1e4, 'LineWidth', 1.6); hold on;
yline(0, 'k:', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Net Heat Flux [W/cm^2]');
title('Net Heat Flux Into the Outer Wall vs. Time (Convective - Radiative)');
grid on; hold off;

% ---- Outer vs inner wall temperature vs time (payload margin view) ----
figure('Name', 'Outer vs Inner Wall Temperature vs Time', 'Color', 'w');
plot(t, res.Touter, 'LineWidth', 1.9); hold on;
plot(t, res.Tinner, 'LineWidth', 1.9);
xlabel('Time [s]'); ylabel('Temperature [K]');
title(sprintf('Outer (Aero) vs. Inner (Payload-Facing) Wall Temperature, x = %.3g m', geom.xPrimary));
legend('Outer wall, T_{outer}', 'Inner wall, T_{inner}', 'Location', 'best');
grid on; hold off;

% ---- Full temperature picture: outer, inner, recovery, edge, reference ----
figure('Name', 'Wall, Recovery, Edge, and Reference Temperature vs Time', 'Color', 'w');
plot(t, res.Touter, 'LineWidth', 1.9); hold on;
plot(t, res.Tinner, 'LineWidth', 1.9);
plot(t, res.Taw,   '--',  'LineWidth', 1.6);
plot(t, res.Te,    ':',   'LineWidth', 1.6);
plot(t, res.Tstar, '-.',  'LineWidth', 1.6);
xlabel('Time [s]'); ylabel('Temperature [K]');
title('Outer/Inner Wall Temperature vs. Adiabatic Wall, Boundary-Layer Edge, and Reference Temperature');
legend('Outer wall, T_{outer}', 'Inner wall, T_{inner}', ...
       'Adiabatic wall (recovery), T_{aw}', 'Boundary-layer edge, T_e', ...
       'Eckert reference temp, T^*', 'Location', 'best');
grid on; hold off;

% ---- Skin friction coefficient vs time (primary station) ----
figure('Name', 'Skin Friction Coefficient vs Time', 'Color', 'w');
plot(t, res.Cf, 'LineWidth', 1.8);
xlabel('Time [s]'); ylabel('Skin Friction Coefficient, C_f [-]');
title(sprintf('Local Skin Friction Coefficient vs. Time (station x = %.3g m)', geom.xPrimary));
grid on;

% ---- Boundary-layer regime indicator (laminar/turbulent), primary station ----
isTurb = strcmp(res.regime, 'turbulent');
figure('Name', 'Boundary-Layer Regime vs Time', 'Color', 'w');
stairs(t, double(isTurb), 'LineWidth', 1.8);
ylim([-0.2, 1.2]);
yticks([0 1]); yticklabels({'Laminar', 'Turbulent'});
xlabel('Time [s]'); ylabel('Boundary-Layer Regime');
title(sprintf('Boundary-Layer Regime vs. Time (station x = %.3g m)', geom.xPrimary));
grid on;

% =====================================================================
% SPATIAL SURVEY: family-of-curves plots vs. body station, one line per
% flight-time snapshot.
% =====================================================================

legendStrs = cell(1, nSnap);
for s = 1:nSnap
    legendStrs{s} = sprintf('t = %.3g s', t(snapIdx(s)));
end

% ---- Outer wall temperature vs body station ----
figure('Name', 'Outer Wall Temperature vs Body Station', 'Color', 'w');
hold on;
for s = 1:nSnap
    plot(x, res.TouterField(snapIdx(s), :), 'LineWidth', 1.8, 'Color', snapColors(s,:));
end
xlabel('Body Station, x [m]'); ylabel('Outer Wall Temperature, T_{outer} [K]');
title('Outer (Aero-Heated) Wall Temperature vs. Body Station, at Several Flight Times');
legend(legendStrs, 'Location', 'best');
grid on; hold off;

% ---- Inner wall temperature vs body station ----
figure('Name', 'Inner Wall Temperature vs Body Station', 'Color', 'w');
hold on;
for s = 1:nSnap
    plot(x, res.TinnerField(snapIdx(s), :), 'LineWidth', 1.8, 'Color', snapColors(s,:));
end
xlabel('Body Station, x [m]'); ylabel('Inner Wall Temperature, T_{inner} [K]');
title('Inner (Payload-Facing) Wall Temperature vs. Body Station, at Several Flight Times');
legend(legendStrs, 'Location', 'best');
grid on; hold off;

% ---- Outer minus inner temperature (through-thickness gradient) vs x ----
figure('Name', 'Through-Thickness Temperature Difference vs Body Station', 'Color', 'w');
hold on;
for s = 1:nSnap
    dT = res.TouterField(snapIdx(s), :) - res.TinnerField(snapIdx(s), :);
    plot(x, dT, 'LineWidth', 1.8, 'Color', snapColors(s,:));
end
xlabel('Body Station, x [m]'); ylabel('T_{outer} - T_{inner} [K]');
title('Through-Thickness Temperature Difference vs. Body Station, at Several Flight Times');
legend(legendStrs, 'Location', 'best');
grid on; hold off;

% ---- Convective heat flux vs body station ----
figure('Name', 'Convective Heat Flux vs Body Station', 'Color', 'w');
hold on;
for s = 1:nSnap
    plot(x, res.qwField(snapIdx(s), :)/1e4, 'LineWidth', 1.8, 'Color', snapColors(s,:));
end
xlabel('Body Station, x [m]'); ylabel('Convective Heat Flux [W/cm^2]');
title('Convective Heat Flux vs. Body Station, at Several Flight Times');
legend(legendStrs, 'Location', 'best');
grid on; hold off;

% ---- Skin friction coefficient vs body station ----
figure('Name', 'Skin Friction Coefficient vs Body Station', 'Color', 'w');
hold on;
for s = 1:nSnap
    plot(x, res.CfField(snapIdx(s), :), 'LineWidth', 1.8, 'Color', snapColors(s,:));
end
xlabel('Body Station, x [m]'); ylabel('Skin Friction Coefficient, C_f [-]');
title('Skin Friction Coefficient vs. Body Station, at Several Flight Times');
legend(legendStrs, 'Location', 'best');
grid on; hold off;

end


% =========================================================================
% LOCAL FUNCTION: multi-station wall energy balance (RHS of the transient
%   ODE). State vector y = [Touter(1..Nx); Tinner(1..Nx)].
% =========================================================================
function dydt = wallEnergyBalanceMulti(t, y, hOf, VOf, alphaOf, xStations, ...
                                        Nx, geom, mat, opts, sigmaSB)

Touter = y(1:Nx);
Tinner = y(Nx+1:2*Nx);

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

dTouter = zeros(Nx,1);
dTinner = zeros(Nx,1);

massAreaHalf = mat.density * (mat.thickness/2);  % kg/m^2, per half-thickness node

for i = 1:Nx
    qw_i = referenceTempHeating(Me, Te, pe, Ue, Touter(i), xStations(i), opts);

    if opts.includeRadiation
        qrad_i = geom.emissivity * sigmaSB * (Touter(i)^4 - T_inf^4);
    else
        qrad_i = 0;
    end

    qcond_i = 2*mat.k*(Touter(i) - Tinner(i)) / mat.thickness;
    qinnerLoss_i = mat.h_internal * (Tinner(i) - mat.T_payload);

    dTouter(i) = (qw_i - qrad_i - qcond_i) / (massAreaHalf * mat.cp);
    dTinner(i) = (qcond_i - qinnerLoss_i)  / (massAreaHalf * mat.cp);
end

dydt = [dTouter; dTinner];

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


% =========================================================================
% LOCAL FUNCTION: 1976 U.S. Standard Atmosphere (valid to 86 km geometric)
%   Returns geopotential height, temperature, pressure, density.
%   Above 86 km, falls back to an exponential extrapolation as a rough
%   approximation only.
% =========================================================================
function [Hgp, T, p, rho] = standardAtmosphere(h_geometric)

Re = 6356766; % m, effective Earth radius used for geopotential conversion
g0 = 9.80665; % m/s^2
Rgas = 287.0528; % J/kg-K, air

Hgp = Re * h_geometric / (Re + h_geometric); % geopotential height

% Layer base data: [base geopotential height (m), base T (K), lapse rate (K/m), base p (Pa)]
% Standard 1976 atmosphere layers up to 86 km
baseH   = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];
baseT   = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946];
lapse   = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, 0];
baseP   = zeros(1, numel(baseH));
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
    % Rough exponential extrapolation above 86 km -- approximate only.
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
if ~isfield(s, field) || isempty(s.(field))
    if isfield(s, field) && isempty(s.(field)) && ...
            (strcmp(field,'halfAngle')) % allow explicit [] to mean "skip"
        return;
    end
    s.(field) = value;
end
end
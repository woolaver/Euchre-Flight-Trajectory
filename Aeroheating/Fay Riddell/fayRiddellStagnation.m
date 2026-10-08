function results = fayRiddellStagnation(fs, geom, wall, opts)
%FAYRIDDELLSTAGNATION  Stagnation-point convective heat flux on a blunted
%nose (sphere) behind the bow shock, using the Fay-Riddell correlation.
%
%   results = FAYRIDDELLSTAGNATION(fs, geom, wall, opts)
%
%   METHOD
%   ---------------------------------------------------------------
%   1. Freestream state from stdatmfull(altitude) (or direct T, p).
%   2. Normal shock ahead of the stagnation point: Rankine-Hugoniot jump
%      with the full enthalpy of the selected gas model.
%   3. Subsonic flow behind the shock is brought to rest isentropically
%      (entropy of the gas model):
%         h0 = h_inf + V^2/2 (conserved),  T0 = T(h0,p0),  s(T0,p0) = s2.
%      T0 and p0 are the boundary-layer-edge (stagnation) conditions.
%   4. Newtonian stagnation-point velocity gradient:
%         due/dx = (1/Rn) * sqrt( 2*(p0 - p_inf) / rho_e )
%   5. Fay-Riddell (axisymmetric stagnation point, isothermal wall):
%         q_w = 0.763 Pr^-0.6 (rho_e mu_e)^0.4 (rho_w mu_w)^0.1
%               * sqrt(due/dx) * (h0 - h_w) * [1 + (Le^a - 1)*hD/h0]
%      a = 0.52 (equilibrium boundary layer) or 0.63 (frozen BL, fully
%      catalytic wall). hD is the dissociation enthalpy carried at the
%      edge (see opts.hD).
%   6. As a sanity check, the Sutton-Graves estimate
%         q = 1.7415e-4 * sqrt(rho_inf/Rn) * V^3
%      is also returned (Earth air, SI units). It should agree with
%      Fay-Riddell to within roughly 10-30 percent.
%
%   INPUT: fs (struct) -- freestream. Scalars or equal-length vectors.
%   ---------------------------------------------------------------
%     fs.altitude  [m]    altitude; T, p, rho from stdatmfull. Default: 30000
%     fs.T, fs.p   [K,Pa] direct freestream state, used ONLY if fs.altitude
%                         is absent or empty.
%     fs.velocity  [m/s]  freestream speed relative to the vehicle.
%                         Default: 3000
%     fs.time      [s]    (optional) time vector, used only as the x-axis
%                         when plotting a vector of conditions.
%
%   INPUT: geom (struct)
%   ---------------------------------------------------------------
%     geom.Rn      [m]    nose radius. Default: 0.01
%
%   INPUT: wall (struct)
%   ---------------------------------------------------------------
%     wall.Tw      [K]    wall temperature (isothermal). Default: 300
%
%   INPUT: opts (struct)
%   ---------------------------------------------------------------
%     opts.gasModel       'equilibrium' (default): 5-species equilibrium air
%                         (N2, O2, NO, N, O) from statistical-mechanics
%                         partition functions (rigid rotor, harmonic
%                         oscillator, electronic levels). Includes O2 and
%                         N2 dissociation and NO formation; no ionization,
%                         so valid to roughly 8000-9000 K.
%                         'vibEq': N2/O2 vibration only, no dissociation.
%                         'perfect': constant gamma = opts.gamma.
%     opts.gamma          Default: 1.4 (gasModel 'perfect' only)
%     opts.R              [J/kg-K] Default: 287
%     opts.Pr             Prandtl number. Default: 0.71
%     opts.Le             Lewis number. Default: 1.4
%     opts.boundaryLayer  'equilibrium' (default) or 'frozen'
%     opts.hD             [J/kg] dissociation enthalpy at the BL edge,
%                         hD = cO*hfO + cN*hfN (atomic O, N mass fractions
%                         times their heats of formation). Default: [] =
%                         computed from the equilibrium edge composition
%                         for gasModel 'equilibrium', 0 for the other
%                         models (which removes the Lewis-number term).
%                         A numeric value overrides in every model.
%     opts.makePlots      plots when a vector of conditions is given.
%                         Default: true
%
%   OUTPUT: results (struct) -- [N x 1] column vectors, SI units
%   ---------------------------------------------------------------
%     altitude, velocity, Tinf, pinf, rhoinf, Minf   freestream
%     T2, p2             static state just behind the normal shock
%     Te, pe, rho_e      edge (stagnation) T, p, rho   (Te = T0, pe = p0)
%     hD                 dissociation enthalpy used in the Lewis term [J/kg]
%     mu_e, h0           edge viscosity, total enthalpy
%     Tw, hw, rho_w, mu_w   wall state
%     dudx               stagnation velocity gradient [1/s]
%     q_w                Fay-Riddell heat flux [W/m^2]
%     q_SuttonGraves     Sutton-Graves check value [W/m^2]
%
%   NOTES / SIMPLIFICATIONS
%   ---------------------------------------------------------------
%   - Gas models 'vibEq' and 'perfect' do not include dissociation, so
%     their T0 is strongly overpredicted (a warning is issued above
%     ~2500 K). Heat flux is far less sensitive: it is driven by the
%     enthalpy difference (h0 - hw), and h0 is exact in every model.
%   - Ionization is not modeled; a warning is issued if T0 exceeds 9000 K.
%   - Viscosity is Sutherland's law for all models, extrapolated above
%     ~2000 K (no mixture transport model for dissociated air).
%
%   Example:
%       fs.altitude = 30000; fs.velocity = 3000;
%       geom.Rn = 0.01; wall.Tw = 300;
%       r = fayRiddellStagnation(fs, geom, wall);
%       fprintf('q_w = %.1f W/cm^2\n', r.q_w/1e4);

% =========================================================================
% 1. INPUT HANDLING / DEFAULTS
% =========================================================================
if nargin < 1 || isempty(fs),   fs   = struct(); end
if nargin < 2 || isempty(geom), geom = struct(); end
if nargin < 3 || isempty(wall), wall = struct(); end
if nargin < 4 || isempty(opts), opts = struct(); end

useAltitude = isfield(fs, 'altitude') && ~isempty(fs.altitude);
if ~useAltitude && ~(isfield(fs, 'T') && isfield(fs, 'p'))
    fs = setDefault(fs, 'altitude', 30000);   % placeholder
    useAltitude = true;
end
fs   = setDefault(fs, 'velocity', 3000);      % placeholder
geom = setDefault(geom, 'Rn', 0.01);          % placeholder
wall = setDefault(wall, 'Tw', 300);           % placeholder

opts = setDefault(opts, 'gasModel',      'equilibrium');
opts = setDefault(opts, 'gamma',         1.4);
opts = setDefault(opts, 'R',             287);
opts = setDefault(opts, 'Pr',            0.71);
opts = setDefault(opts, 'Le',            1.4);
opts = setDefault(opts, 'boundaryLayer', 'equilibrium');
opts = setDefault(opts, 'hD',            []);
opts = setDefault(opts, 'makePlots',     true);

switch lower(opts.boundaryLayer)
    case 'equilibrium', alphaLe = 0.52;
    case 'frozen',      alphaLe = 0.63;
    otherwise
        error('fayRiddellStagnation:badBoundaryLayer', ...
            'opts.boundaryLayer must be ''equilibrium'' or ''frozen'' (got ''%s'').', ...
            opts.boundaryLayer);
end

R   = opts.R;
gas = buildGasModel(opts);
useEq = strcmpi(opts.gasModel, 'equilibrium');
if useEq
    sp = airSpecies();
end

% Expand scalars so every input is an N x 1 vector
if useAltitude
    sizes = [numel(fs.altitude), numel(fs.velocity), numel(geom.Rn), numel(wall.Tw)];
else
    sizes = [numel(fs.T), numel(fs.p), numel(fs.velocity), numel(geom.Rn), numel(wall.Tw)];
end
N = max(sizes);
V  = expandTo(fs.velocity, N, 'fs.velocity');
Rn = expandTo(geom.Rn,     N, 'geom.Rn');
Tw = expandTo(wall.Tw,     N, 'wall.Tw');
if useAltitude
    alt = expandTo(fs.altitude, N, 'fs.altitude');
else
    Tin = expandTo(fs.T, N, 'fs.T');
    pin = expandTo(fs.p, N, 'fs.p');
    alt = nan(N,1);
end

% Preallocate outputs
Tinf = zeros(N,1); pinf = zeros(N,1); rhoinf = zeros(N,1); Minf = zeros(N,1);
T2 = zeros(N,1); p2 = zeros(N,1);
Te = zeros(N,1); pe = zeros(N,1); rho_e = zeros(N,1); mu_e = zeros(N,1); h0 = zeros(N,1);
hw = zeros(N,1); rho_w = zeros(N,1); mu_w = zeros(N,1);
dudx = zeros(N,1); q_w = zeros(N,1); q_SG = zeros(N,1); hDout = zeros(N,1);

% =========================================================================
% 2. STAGNATION-POINT HEATING AT EACH CONDITION
% =========================================================================
for k = 1:N
    % --- Freestream ---
    if useAltitude
        [rhoinf(k), pinf(k), Tinf(k)] = stdatmfull(alt(k));
        if pinf(k) <= 0
            error('fayRiddellStagnation:altitudeTooHigh', ...
                'stdatmfull returns zero pressure at %.4g m (above its 500 km table).', alt(k));
        end
    else
        Tinf(k) = Tin(k); pinf(k) = pin(k);
        rhoinf(k) = pinf(k) / (R * Tinf(k));
    end
    Minf(k) = V(k) / sqrt(gas.gamma(Tinf(k)) * R * Tinf(k));

    if useEq
        % --- Equilibrium air: normal shock + isentropic stop to rest ---
        st = eqStagnation(V(k), Tinf(k), pinf(k), sp);
        T2(k) = st.T2;  p2(k) = st.p2;
        h0(k) = st.h0;  Te(k) = st.T0;
        pe(k) = st.p0;  rho_e(k) = st.rho_e;
        hD_k  = st.hD;
        [hw(k), rho_w(k)] = eqThermo(Tw(k), pe(k), sp);
    else
        % --- Normal shock (normal-shock velocity = freestream velocity) ---
        [~, T2(k), p2(k)] = normalShockJump(V(k), Tinf(k), pinf(k), gas, R);

        % --- Isentropic deceleration to rest: edge (stagnation) state ---
        h0(k) = gas.h(Tinf(k)) + 0.5*V(k)^2;      % total enthalpy, conserved
        Te(k) = gas.T(h0(k));
        pe(k) = p2(k) * exp((gas.s(Te(k)) - gas.s(T2(k))) / R);
        rho_e(k) = pe(k) / (R * Te(k));
        hD_k  = 0;

        % --- Wall state ---
        hw(k)    = gas.h(Tw(k));
        rho_w(k) = pe(k) / (R * Tw(k));
    end
    mu_e(k) = sutherlandViscosity(Te(k));
    mu_w(k) = sutherlandViscosity(Tw(k));
    if ~isempty(opts.hD)
        hD_k = opts.hD;                           % user override
    end
    hDout(k) = hD_k;

    % --- Newtonian stagnation-point velocity gradient ---
    dudx(k) = sqrt(2 * (pe(k) - pinf(k)) / rho_e(k)) / Rn(k);

    % --- Fay-Riddell ---
    lewisFactor = 1 + (opts.Le^alphaLe - 1) * hD_k / h0(k);
    q_w(k) = 0.763 * opts.Pr^(-0.6) ...
           * (rho_e(k) * mu_e(k))^0.4 * (rho_w(k) * mu_w(k))^0.1 ...
           * sqrt(dudx(k)) * (h0(k) - hw(k)) * lewisFactor;

    % --- Sutton-Graves cross-check ---
    q_SG(k) = 1.7415e-4 * sqrt(rhoinf(k) / Rn(k)) * V(k)^3;
end

% =========================================================================
% 3. PACKAGE RESULTS, WARN, REPORT, PLOT
% =========================================================================
results.altitude = alt;      results.velocity = V;
results.Tinf = Tinf;         results.pinf = pinf;
results.rhoinf = rhoinf;     results.Minf = Minf;
results.T2 = T2;             results.p2 = p2;
results.Te = Te;             results.pe = pe;
results.rho_e = rho_e;       results.mu_e = mu_e;
results.h0 = h0;
results.Tw = Tw;             results.hw = hw;
results.rho_w = rho_w;       results.mu_w = mu_w;
results.dudx = dudx;
results.hD = hDout;
results.q_w = q_w;
results.q_SuttonGraves = q_SG;

if ~useEq && max(Te) > 2500
    warning('fayRiddellStagnation:dissociation', ...
        ['Stagnation temperature reaches %.0f K but gasModel ''%s'' has no ' ...
         'dissociation, so T0, rho_e and p0 are off (use gasModel ''equilibrium''). ' ...
         'The heat flux is much less sensitive than the temperatures.'], ...
        max(Te), opts.gasModel);
end
if useEq && max(Te) > 9000
    warning('fayRiddellStagnation:ionization', ...
        ['Stagnation temperature reaches %.0f K. Ionization is not modeled ' ...
         'and becomes significant above ~9000 K.'], max(Te));
end

if N == 1
    fprintf('Fay-Riddell stagnation-point heating\n');
    fprintf('  Freestream: M = %.2f, V = %.0f m/s, T = %.1f K, p = %.4g Pa\n', ...
        Minf, V, Tinf, pinf);
    fprintf('  Stagnation: T0 = %.0f K, p0 = %.4g Pa, due/dx = %.4g 1/s\n', Te, pe, dudx);
    fprintf('  q_w (Fay-Riddell)  = %.2f W/cm^2\n', q_w/1e4);
    fprintf('  q   (Sutton-Graves) = %.2f W/cm^2\n', q_SG/1e4);
elseif opts.makePlots
    if isfield(fs, 'time') && numel(fs.time) == N
        xAxis = fs.time(:); xLabel = 'Time [s]';
    else
        xAxis = (1:N).'; xLabel = 'Sample index';
    end

    figure('Name', 'Stagnation-Point Heat Flux');
    plot(xAxis, q_w/1e4, 'LineWidth', 1.8); hold on;
    %plot(xAxis, q_SG/1e4, '--', 'LineWidth', 1.6);
    xlabel(xLabel); ylabel('Stagnation-Point Heat Flux [W/cm^2]');
    title('Stagnation-Point Heat Flux');
    %legend('Fay-Riddell', 'Sutton-Graves check', 'Location', 'best');
    grid on; hold off;

    figure('Name', 'Stagnation Temperature');
    plot(xAxis, Te, 'LineWidth', 1.8); hold on;
    plot(xAxis, Tw, '--', 'LineWidth', 1.6);
    xlabel(xLabel); ylabel('Temperature [K]');
    title('Edge (Stagnation) and Wall Temperature');
    legend('Edge, T_0', 'Wall, T_w', 'Location', 'best');
    grid on; hold off;
end

end % ================== END MAIN FUNCTION ==================


% =========================================================================
% LOCAL FUNCTION: gas model (handles for h, T(h), gamma, entropy)
%   h [J/kg] measured from h(0) = 0.  s(T) [J/kg-K] is the temperature-only
%   part of the entropy, s_total = s(T) - R*ln(p).  Scalars only.
% =========================================================================
function gas = buildGasModel(opts)

R = opts.R;

switch lower(opts.gasModel)
    case 'perfect'
        cp0 = opts.gamma * R / (opts.gamma - 1);
        gas.h     = @(T) cp0 * T;
        gas.T     = @(h) h / cp0;
        gas.gamma = @(T) opts.gamma;
        gas.s     = @(T) cp0 * log(T);

    case {'vibeq', 'equilibrium'}   % 'equilibrium' uses these only for M_inf
        thetaV = [3395, 2239];   % N2, O2 characteristic vibrational temps [K]
        xMol   = [0.79, 0.21];   % N2, O2 mole fractions
        % Translation + rotation: cp/R = 7/2. Vibration: harmonic oscillator.
        cpOf = @(T) R * (3.5 + sum(xMol .* vibCpTerm(thetaV ./ T)));
        hOf  = @(T) R * (3.5*T + sum(xMol .* thetaV ./ expm1(thetaV ./ T)));
        sOf  = @(T) R * (3.5*log(T) + sum(xMol .* ...
                   ((thetaV./T) ./ expm1(thetaV./T) - log1p(-exp(-thetaV./T)))));
        gas.h     = hOf;
        gas.T     = @(h) temperatureFromEnthalpy(h, hOf, cpOf, R);
        gas.gamma = @(T) cpOf(T) / (cpOf(T) - R);
        gas.s     = sOf;

    otherwise
        error('fayRiddellStagnation:badGasModel', ...
            ['opts.gasModel must be ''equilibrium'', ''vibEq'' or ''perfect'' ' ...
             '(got ''%s'').'], opts.gasModel);
end

end


% Vibrational contribution to cp/R for one mode, u = theta_v/T.
% Written with exp(-u) so it does not overflow at low temperature.
function f = vibCpTerm(u)
f = u.^2 .* exp(-u) ./ (1 - exp(-u)).^2;
end


% Invert h(T) by Newton's method. h is convex and increasing in T and
% h(T) >= 3.5*R*T, so starting at T = h/(3.5 R) (an upper bound on the
% root) gives monotone convergence.
function T = temperatureFromEnthalpy(hT, hOf, cpOf, R)
if hT <= 0
    T = 1;
    return;
end
T = hT / (3.5*R);
for it = 1:100
    dT = (hOf(T) - hT) / cpOf(T);
    T  = T - dT;
    if abs(dT) < 1e-10 * T
        break;
    end
end
end


% =========================================================================
% LOCAL FUNCTION: normal-shock jump for a gas with h(T)
%   Conservation of mass, momentum and energy across a normal shock with
%   upstream normal velocity u1n. Unknown: eps = rho1/rho2, solved by
%   Newton's method starting from the perfect-gas value.
% =========================================================================
function [eps2, T2, p2] = normalShockJump(u1n, T1, p1, gas, R)

rho1 = p1 / (R*T1);
g1   = gas.gamma(T1);
a1   = sqrt(g1 * R * T1);

if u1n <= a1                     % no shock (Mach <= 1)
    eps2 = 1; T2 = T1; p2 = p1;
    return;
end

h1  = gas.h(T1);
phi = @(e) rho1 * R * gas.T(h1 + 0.5*u1n^2*(1 - e^2)) ...
           / (p1 + rho1*u1n^2*(1 - e));

Mn2  = (u1n/a1)^2;
eps2 = ((g1-1)*Mn2 + 2) / ((g1+1)*Mn2);

for it = 1:50
    G  = phi(eps2) - eps2;
    d  = 1e-7;
    dG = (phi(eps2 + d) - (eps2 + d) - G) / d;
    if dG == 0
        break;
    end
    epsNew = min(max(eps2 - G/dG, 0.02), 1);
    if abs(epsNew - eps2) < 1e-12
        eps2 = epsNew;
        break;
    end
    eps2 = epsNew;
end

T2 = gas.T(h1 + 0.5*u1n^2*(1 - eps2^2));
p2 = p1 + rho1*u1n^2*(1 - eps2);

end


% =========================================================================
% EQUILIBRIUM AIR (N2, O2, NO, N, O) FROM PARTITION FUNCTIONS
%   Species order everywhere: 1 N2, 2 O2, 3 NO, 4 N, 5 O.
%   Energy zero: ground state of N2 and O2 molecules (so cold air has
%   h ~ 3.5*R*T like the other gas models). Units: eps0, theta in K.
% =========================================================================
function sp = airSpecies()
sp.kB = 1.380649e-23;       % J/K
sp.h  = 6.62607015e-34;     % J*s
sp.NA = 6.02214076e23;      % 1/mol
sp.Ru = 8.314462618;        % J/mol-K
sp.M  = [28.0134, 31.9988, 30.0061, 14.0067, 15.9994] * 1e-3;  % kg/mol
D0    = [113250, 59370, 75396];  % D0/k for N2, O2, NO [K]
% Energy of each species above N2/O2 ground state (heats of formation at 0 K)
sp.eps0 = [0, 0, 0.5*(D0(1) + D0(2)) - D0(3), 0.5*D0(1), 0.5*D0(2)];
sp.thv  = [3395, 2239, 2817];    % vibrational temperatures [K]
sp.thr  = [2.88, 2.07, 2.45];    % rotational temperatures [K]
sp.sig  = [2, 2, 1];             % rotational symmetry numbers
% Electronic levels: degeneracies and energies (K)
sp.elecG  = {1, [3 2 1], [2 2], [4 10 6], [5 3 1 5 1]};
sp.elecTh = {0, [0 11390 18985], [0 174.2], [0 27660 41490], [0 227.7 326.6 22830 48620]};
% Undissociated air composition (mole fractions, argon neglected)
sp.xN2 = 0.79;
sp.xO2 = 0.21;
end


% Electronic partition function (log) and mean electronic energy / k [K]
function [lnqe, Eel] = electronic(i, T, sp)
g  = sp.elecG{i};
th = sp.elecTh{i};
w  = g .* exp(-th ./ T);
q  = sum(w);
lnqe = log(q);
Eel  = sum(w .* th) / q;
end


% ln of the partition function per unit volume, including the energy zero
function v = lnQspecies(i, T, sp)
m    = sp.M(i) / sp.NA;
lnqt = 1.5 * log(2*pi*m*sp.kB*T / sp.h^2);
lnqe = electronic(i, T, sp);
v = lnqt + lnqe - sp.eps0(i)/T;
if i <= 3   % diatomic: rigid rotor + harmonic oscillator
    v = v + log(T / (sp.sig(i)*sp.thr(i))) - log1p(-exp(-sp.thv(i)/T));
end
end


% Equilibrium mole fractions x = [N2 O2 NO N O] at temperature T [K] and
% pressure p [Pa]. Law of mass action for O2<->2O, N2<->2N, NO<->N+O plus
% the air N:O atom ratio and total pressure; damped Newton in ln(pO), ln(pN).
function x = eqComposition(T, p, sp)

lQ = zeros(1, 5);
for i = 1:5
    lQ(i) = lnQspecies(i, T, sp);
end
lkT  = log(sp.kB * T);
lKO2 = lkT + 2*lQ(5) - lQ(2);      % pO^2/pO2   [Pa]
lKN2 = lkT + 2*lQ(4) - lQ(1);      % pN^2/pN2   [Pa]
lKNO = lkT + lQ(4) + lQ(5) - lQ(3);% pN*pO/pNO  [Pa]
ratio = sp.xN2 / sp.xO2;           % N-atom : O-atom ratio in air

% Starting guess: undissociated air, capped so atoms cannot exceed p
a = min(0.5*(lKO2 + log(sp.xO2*p)), log(0.2*p));   % ln pO
b = min(0.5*(lKN2 + log(sp.xN2*p)), log(0.7*p));   % ln pN

for it = 1:100
    [f, J] = eqResidual(a, b, lKO2, lKN2, lKNO, p, ratio);
    if max(abs(f)) < 1e-12
        break;
    end
    d  = -(J \ f);
    mx = max(abs(d));
    if mx > 2
        d = d * 2/mx;               % limit the step in log space
    end
    n0 = norm(f);
    t  = 1;
    while t > 1e-4                  % backtracking line search
        fn = eqResidual(a + t*d(1), b + t*d(2), lKO2, lKN2, lKNO, p, ratio);
        if norm(fn) < n0 || max(abs(fn)) < 1e-12
            break;
        end
        t = t/2;
    end
    a = a + t*d(1);
    b = b + t*d(2);
end

[~, ~, pp] = eqResidual(a, b, lKO2, lKN2, lKNO, p, ratio);
x = pp / sum(pp);

end


% Residuals f = [ln(N-atoms/O-atoms/ratio); ln(Ptot/p)], Jacobian wrt
% [ln pO, ln pN], and partial pressures pp = [pN2 pO2 pNO pN pO].
function [f, J, pp] = eqResidual(a, b, lKO2, lKN2, lKNO, p, ratio)
pO  = exp(a);
pN  = exp(b);
pO2 = exp(2*a - lKO2);
pN2 = exp(2*b - lKN2);
pNO = exp(a + b - lKNO);
Nn  = pN + 2*pN2 + pNO;
Oo  = pO + 2*pO2 + pNO;
P   = pO + pN + pO2 + pN2 + pNO;
f   = [log(Nn / (Oo*ratio)); log(P / p)];
J   = [ pNO/Nn - (pO + 4*pO2 + pNO)/Oo,  (pN + 4*pN2 + pNO)/Nn - pNO/Oo; ...
        (pO + 2*pO2 + pNO)/P,             (pN + 2*pN2 + pNO)/P ];
pp  = [pN2, pO2, pNO, pN, pO];
end


% Equilibrium-air enthalpy h [J/kg], density rho [kg/m^3] and (only if
% requested) entropy s [J/kg-K], composition x and mean molar mass Mb.
function [h, rho, s, x, Mb] = eqThermo(T, p, sp)

x  = eqComposition(T, p, sp);
Mb = sum(x .* sp.M);

H = zeros(1, 5);                   % molar enthalpy / R_u [K]
for i = 1:5
    [~, Eel] = electronic(i, T, sp);
    if i <= 3
        u = sp.thv(i) / T;
        H(i) = 3.5*T + u*T/expm1(u) + sp.eps0(i) + Eel;
    else
        H(i) = 2.5*T + sp.eps0(i) + Eel;
    end
end
h   = sp.Ru * sum(x .* H) / Mb;
rho = p * Mb / (sp.Ru * T);

s = [];
if nargout >= 3
    S = 0;                         % entropy per mole of mixture / R_u
    for i = 1:5
        if x(i) > 0
            [lnqe, Eel] = electronic(i, T, sp);
            m   = sp.M(i) / sp.NA;
            s_i = 1.5*log(2*pi*m*sp.kB*T/sp.h^2) + log(sp.kB*T/(x(i)*p)) ...
                  + 2.5 + lnqe + Eel/T;
            if i <= 3
                u   = sp.thv(i) / T;
                s_i = s_i + log(T/(sp.sig(i)*sp.thr(i))) + 1 ...
                      - log1p(-exp(-u)) + u/expm1(u);
            end
            S = S + x(i) * s_i;
        end
    end
    s = sp.Ru * S / Mb;
end

end


% Invert h(T, p) for T at fixed pressure
function T = eqTemperatureFromEnthalpy(hT, p, sp)
T = fzero(@(Tt) eqThermo(Tt, p, sp) - hT, [150, 25000], optimset('TolX', 1e-9));
end


% Equilibrium normal shock: upstream speed u1 (normal to the shock),
% unknown eps = rho1/rho2 from mass + momentum + energy conservation.
function [eps2, T2, p2] = eqNormalShock(u1, T1, p1, sp)

[h1, rho1] = eqThermo(T1, p1, sp);
G = @(e) eqShockResidual(e, u1, p1, h1, rho1, sp);

if G(0.03) * G(0.9) > 0
    error('fayRiddellStagnation:noShock', ...
        'No normal-shock solution found: is the freestream speed (%.4g m/s) supersonic?', u1);
end
eps2 = fzero(G, [0.03, 0.9], optimset('TolX', 1e-12));

p2 = p1 + rho1*u1^2*(1 - eps2);
h2 = h1 + 0.5*u1^2*(1 - eps2^2);
T2 = eqTemperatureFromEnthalpy(h2, p2, sp);

end


function r = eqShockResidual(e, u1, p1, h1, rho1, sp)
p2 = p1 + rho1*u1^2*(1 - e);
h2 = h1 + 0.5*u1^2*(1 - e^2);
T2 = eqTemperatureFromEnthalpy(h2, p2, sp);
[~, rho2] = eqThermo(T2, p2, sp);
r = rho1/rho2 - e;
end


% Equilibrium stagnation state behind a normal shock, plus the
% dissociation enthalpy hD = cO*hfO + cN*hfN at the edge.
function st = eqStagnation(V, T1, p1, sp)

[eps2, T2, p2] = eqNormalShock(V, T1, p1, sp); %#ok<ASGLU>
h1 = eqThermo(T1, p1, sp);
h0 = h1 + 0.5*V^2;                          % total enthalpy, conserved
[~, ~, s2] = eqThermo(T2, p2, sp);

% Isentropic stop: find p0 such that s(T(h0,p0), p0) = s2
Rs   = @(lnp0) eqEntropyResidual(lnp0, h0, s2, sp);
lnp0 = fzero(Rs, [log(p2), log(p2) + 0.5], optimset('TolX', 1e-12));
p0   = exp(lnp0);
T0   = eqTemperatureFromEnthalpy(h0, p0, sp);
[~, rho_e, ~, x, Mb] = eqThermo(T0, p0, sp);

cO  = x(5) * sp.M(5) / Mb;                  % atomic O mass fraction
cN  = x(4) * sp.M(4) / Mb;                  % atomic N mass fraction
hfO = sp.eps0(5) * sp.Ru / sp.M(5);         % heat of formation of O [J/kg]
hfN = sp.eps0(4) * sp.Ru / sp.M(4);         % heat of formation of N [J/kg]

st.T2 = T2;   st.p2 = p2;
st.T0 = T0;   st.p0 = p0;
st.h0 = h0;   st.rho_e = rho_e;
st.hD = cO*hfO + cN*hfN;

end


function r = eqEntropyResidual(lnp0, h0, s2, sp)
p0 = exp(lnp0);
T0 = eqTemperatureFromEnthalpy(h0, p0, sp);
[~, ~, s] = eqThermo(T0, p0, sp);
r = s - s2;
end


% =========================================================================
% LOCAL FUNCTION: Sutherland viscosity
% =========================================================================
function mu = sutherlandViscosity(T)
mu0 = 1.716e-5;   % kg/(m-s)
T0  = 273.15;     % K
S   = 110.4;      % K
mu  = mu0 * (T/T0)^1.5 * (T0 + S) / (T + S);
end


% =========================================================================
% LOCAL FUNCTION: standard atmosphere (user-supplied)
% =========================================================================
function [rho,p,T] = stdatmfull(h)
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
% LOCAL FUNCTION: scalar -> N x 1 expansion with a size check
% =========================================================================
function v = expandTo(x, N, name)
x = x(:);
if numel(x) == 1
    v = repmat(x, N, 1);
elseif numel(x) == N
    v = x;
else
    error('fayRiddellStagnation:sizeMismatch', ...
        '%s has %d elements but the other inputs imply %d.', name, numel(x), N);
end
end


% =========================================================================
% LOCAL FUNCTION: fill in default struct fields
% =========================================================================
function s = setDefault(s, field, value)
if ~isfield(s, field)
    s.(field) = value;
end
end

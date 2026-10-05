function [CD, CL, CN, CA, Mach, results] = hypersonicDragCoeff(altitude, velocity, alpha, geom)
%HYPERSONICDRAGCOEFF Estimate hypersonic drag/lift coefficients using
%   Modified Newtonian Theory (MNT) for a slender conical body.
%
%   [CD, CL, CN, CA, Mach, results] = hypersonicDragCoeff(altitude, velocity, alpha, geom)
%
%   METHOD
%   ------
%   Modified Newtonian theory approximates local surface pressure as:
%       Cp(theta_local) = Cp_max * sin(theta_local)^2
%   where theta_local is the local surface inclination to the freestream,
%   and Cp_max is the stagnation-point pressure coefficient computed from
%   normal-shock / isentropic stagnation relations (valid for calorically
%   perfect gas, gamma = 1.4, at the freestream Mach number).
%
%   This is a standard, widely published approximation (see J.D. Anderson,
%   "Hypersonic and High-Temperature Gas Dynamics") used for preliminary
%   estimates of pressure-dominated (inviscid) hypersonic aerodynamics on
%   simple shapes. It excludes skin friction, real-gas/chemistry effects,
%   and viscous-interaction/rarefaction corrections, so results are only
%   an order-of-magnitude / early-design estimate, not a substitute for
%   CFD or wind-tunnel/flight data.
%
%   INPUTS
%   ------
%   altitude  : geometric altitude above sea level [m]
%   velocity  : freestream velocity magnitude [m/s]
%   alpha     : angle of attack [deg], angle between body axis and
%               freestream velocity vector
%   geom      : struct describing the body (all fields optional; defaults
%               shown). Body is modeled as a sharp right-circular cone
%               with a flat base (base drag neglected/optional).
%       .halfAngle_deg   cone half-angle [deg]                (default 10)
%       .length          axial length of cone [m]              (default 1)
%       .baseRadius      base radius [m] (overrides halfAngle  (default [])
%                         if both length and baseRadius given)
%       .refArea         reference area for coefficients [m^2] (default:
%                         base circular cross-section, pi*R_base^2)
%       .includeBaseDrag logical, add crude base-drag term      (default true)
%       .nPanels         number of azimuthal panels for the     (default 360)
%                         numerical surface integration
%
%   OUTPUTS
%   -------
%   CD       : drag coefficient (freestream-wind-axis), referenced to refArea
%   CL       : lift coefficient (freestream-wind-axis), referenced to refArea
%   CN       : normal-force coefficient (body axes)
%   CA       : axial-force coefficient (body axes)
%   Mach     : freestream Mach number computed from altitude/velocity
%   results  : struct with intermediate quantities (rho, T, a, Cp_max,
%              q_inf, geometry used, etc.) for inspection/debugging
%
%   EXAMPLE
%   -------
%   geom.halfAngle_deg = 8;
%   geom.length = 2.5;
%   [CD, CL] = hypersonicDragCoeff(30000, 2200, 3, geom);
%
%   NOTE: This function estimates aerodynamic *coefficients* for
%   generic educational/engineering analysis (e.g., glide range,
%   trajectory-shaping studies). It does not compute or output any
%   guidance, control, or weaponization-specific quantities.

    %% ---- Handle inputs / defaults ----
    if nargin < 4 || isempty(geom)
        geom = struct();
    end
    geom = setDefault(geom, 'halfAngle_deg', 10);
    geom = setDefault(geom, 'length', 1);
    geom = setDefault(geom, 'baseRadius', []);
    geom = setDefault(geom, 'includeBaseDrag', true);
    geom = setDefault(geom, 'nPanels', 360);

    if isempty(geom.baseRadius)
        geom.baseRadius = geom.length * tand(geom.halfAngle_deg);
    end
    geom = setDefault(geom, 'refArea', pi * geom.baseRadius^2);

    alpha_rad = deg2rad(alpha);
    theta_c   = deg2rad(geom.halfAngle_deg); % cone half-angle

    %% ---- Atmosphere model (1976 US Standard Atmosphere, simplified) ----
    [T, a, ~, rho] = atmoscoesa(altitude);
    Mach = velocity / a;
    q_inf = 0.5 * rho * velocity^2; % dynamic pressure

    if Mach < 5
        warning('hypersonicDragCoeff:lowMach', ...
            ['Mach = %.2f. Modified Newtonian theory is intended for ' ...
             'hypersonic flow (Mach > ~5); accuracy degrades below that.'], Mach);
    end

    %% ---- Stagnation pressure coefficient (Cp_max) ----
    % Rayleigh Pitot-tube formula behind a normal shock (gamma = 1.4),
    % valid for calorically perfect gas. This gives Cp_max = (p0'-p_inf)/q_inf
    gamma = 1.4;
    Cp_max = newtonianCpMax(Mach, gamma);

    %% ---- Integrate pressure over cone surface (Modified Newtonian) ----
    % Geometry: cone apex at origin, axis along body +x, base at x = L.
    % A surface point at axial station x and azimuth phi (phi=0 along body +y)
    % has position (x, r*cos(phi), r*sin(phi)) with r(x) = x*tan(theta_c).
    % The outward unit surface normal (derived from the cross product of the
    % two surface tangent directions, dP/dx and dP/dphi) is:
    %   n_hat = [ -sin(theta_c) ; cos(theta_c)*cos(phi) ; cos(theta_c)*sin(phi) ]
    % which is constant with x (a cone's normal only depends on azimuth).
    %
    % Angle of attack alpha is taken as a rotation of the freestream direction
    % in the body x-y plane (i.e., phi=0 is the windward ray for alpha>0):
    %   V_hat = [ cos(alpha) ; sin(alpha) ; 0 ]
    %
    % Local impingement/deflection factor c(phi) = -(n_hat . V_hat):
    %   c(phi) = sin(theta_c)*cos(alpha) - cos(theta_c)*sin(alpha)*cos(phi)
    % Flow impinges on the surface where c(phi) > 0 (windward side); on the
    % leeward side (c<=0) Newtonian theory sets Cp = 0 (no aerodynamic shadow
    % contribution, consistent with standard Modified Newtonian practice).
    % Local pressure coefficient: Cp(phi) = Cp_max * max(c(phi),0)^2

    phi = linspace(0, 2*pi, geom.nPanels+1);
    phi = phi(1:end-1); % avoid duplicate 0/2pi point
    dphi = 2*pi / geom.nPanels;

    sinC = sin(theta_c);
    cosC = cos(theta_c);

    c = sinC.*cos(alpha_rad) - cosC.*sin(alpha_rad).*cos(phi);
    impinge = max(c, 0);
    Cp_local = Cp_max .* (impinge.^2);

    % --- Force integration over the conical surface ---
    % Surface-area element: dA = r(x) * (dx/cos(theta_c)) * dphi, r(x)=x*tan(theta_c)
    % Force on the body from pressure: dF = -Cp_local * q_inf * n_hat * dA
    % Axial-force coefficient (body +x, rearward = drag-producing):
    %   CA = (1/refArea) * integral[ Cp_local * sin(theta_c) dA ]   (n_hat_x = -sin(theta_c))
    % Normal-force coefficient (body +y, in the angle-of-attack plane):
    %   CN = -(1/refArea) * integral[ Cp_local * cos(theta_c)*cos(phi) dA ]
    % Sign conventions above are chosen so CA>0 opposes forward motion and
    % CN>0 produces positive lift for positive angle of attack.

    L = geom.length;
    tanC = tan(theta_c);

    % x-integral of r(x)/cos(theta_c) dx over x=[0,L]: = tan(theta_c)*L^2/(2*cos(theta_c))
    xIntegralFactor = tanC * (L^2/2) / cosC;

    CA =  sum(Cp_local .* sinC)            * dphi * xIntegralFactor / geom.refArea;
    CN = -sum(Cp_local .* cosC .* cos(phi)) * dphi * xIntegralFactor / geom.refArea;

    %% ---- Optional crude base drag (blunt base at freestream static pressure deficit) ----
    if geom.includeBaseDrag
        % Simple empirical hypersonic base-pressure estimate: Cp_base ~ -1/Mach^2
        % (crude correlation; replace with empirical data for real designs)
        Cp_base = -1 / Mach^2;
        Abase = pi * geom.baseRadius^2;
        CA_base = -Cp_base * (Abase / geom.refArea); % base pressure acts rearward-facing surface
        CA = CA + CA_base;
    else
        CA_base = 0;
    end

    %% ---- Rotate body-axis (CA, CN) into wind-axis (CD, CL) ----
    CD = CA .* cos(alpha_rad) + CN .* sin(alpha_rad);
    CL = CN .* cos(alpha_rad) - CA .* sin(alpha_rad);

    %% ---- Package intermediate results ----
    results = struct( ...
        'rho', rho, 'temperature_K', T, 'speedOfSound', a, ...
        'q_inf', q_inf, 'Cp_max', Cp_max, 'gamma', gamma, ...
        'geometry', geom, 'CA_base', CA_base, 'alpha_deg', alpha, ...
        'note', 'Modified Newtonian Theory estimate; inviscid, no real-gas effects.');

end

% ==================== helper functions ====================

function s = setDefault(s, field, value)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = value;
    end
end

function Cp_max = newtonianCpMax(M, gamma)
%NEWTONIANCPMAX Stagnation-point pressure coefficient behind a normal
%   shock (Rayleigh Pitot formula), calorically perfect gas.
    if M <= 1
        % Subsonic/transonic isentropic stagnation Cp (not the intended regime,
        % included only to avoid complex results if called at low Mach)
        Cp_max = (2/(gamma*M^2)) * ((1 + (gamma-1)/2*M^2)^(gamma/(gamma-1)) - 1);
        return
    end
    % Pressure ratio across normal shock
    p2_p1 = (2*gamma*M^2 - (gamma-1)) / (gamma+1);
    % Stagnation pressure ratio behind shock relative to freestream static (Rayleigh Pitot)
    p02_p1 = p2_p1 * ( ((gamma+1)^2 * M^2) / (4*gamma*M^2 - 2*(gamma-1)) )^(gamma/(gamma-1));
    Cp_max = (2/(gamma*M^2)) * (p02_p1 - 1);
end

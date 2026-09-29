function [CL, CD, Cp] = modifiedNewtonian(segs, normals, alpha, Mach, refLength)
% MNTCLCD Modified Newtonian Theory sweep over a 2-D sliced outline.
%
% [CL, CD, CA, CN, Cp] = mntCLCD(segs, normals, alphaDeg, CpMax, refLength)
%
% Inputs:
%   segs      - Ns x 4 matrix from sliceSTLPlane: [c1a c2a c1b c2b] per
%               panel, with c1 taken as the body's flow-aligned
%               (axial-ish) direction and c2 the transverse direction in
%               the cut plane
%   normals   - Ns x 2 matrix from sliceSTLPlane: outward unit normal
%               per panel in the (c1,c2) plane. Rows of NaN (ambiguous
%               normal) are skipped.
%   alpha - angle of attack, degrees. Freestream direction is taken
%               as [cos(alpha), sin(alpha)] in (c1,c2).
%   refLength - reference length for nondimensionalizing the panel sum
%               (e.g. body length along c1). Keep consistent with
%               whatever convention the rest of your framework uses.
%
% Outputs:
%   CL, CD - lift and drag coefficients (per unit span -- see caveat
%            below), in wind axes
%   Cp     - Ns x 1 vector of per-panel pressure coefficients, useful
%            for plotting the pressure distribution or debugging
%
% Method: for each panel, sin(theta) = dot(freestream_dir, outward_normal).
% Windward panels (sin(theta) > 0) get Cp = CpMax * sin(theta)^2;
% shadowed panels (sin(theta) <= 0) get Cp = 0 (the basic Newtonian
% shadowing rule -- note this is where a Prandtl-Meyer leeside
% correction would improve accuracy instead of just zeroing Cp out, per
% the same caveat flagged for the panel-method work). Force on each
% panel acts along -outward_normal (pressure pushes into the body); sum
% and rotate CA/CN by alpha to get CL/CD.\

alpha = deg2rad(alpha);
Vinf  = [cos(alpha), sin(alpha)];

gamma = 1.4;
CpMax = (2/(gamma*Mach^2)) * (((((gamma+1)^2*Mach^2)/(4*gamma*Mach^2-2*(gamma-1)))^(gamma/(gamma-1))) * ((1-gamma+2*gamma*Mach^2)/(gamma+1)) - 1);

Ns = size(segs,1);
Cp = nan(Ns,1);
CA = 0;
CN = 0;

for i = 1:Ns
    n = normals(i,:);
    if any(isnan(n))
        continue % ambiguous normal for this panel -- skip
    end

    p1 = segs(i,1:2);
    p2 = segs(i,3:4);
    dL = hypot(p2(1)-p1(1), p2(2)-p1(2));

    sinTheta = dot(Vinf, n);
    if sinTheta > 0
        Cp(i) = CpMax * sinTheta^2;
    else
        Cp(i) = 0;
    end

    dF = -Cp(i) * n * dL;   % force contribution, [c1, c2] components
    CA = CA + dF(1);
    CN = CN + dF(2);
end

CA = CA / refLength;
CN = CN / refLength;

CL = CN*cos(alpha) - CA*sin(alpha);
CD = CN*sin(alpha) + CA*cos(alpha);

end

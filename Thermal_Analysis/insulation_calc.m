function R =  insulation_calc(kval,thick)
%takes a matrix of k values and matrix of thicknesses

if ~isequal(size(kval), size(thick))
    error('kval and thick must have the same dimensions.');
end
R = sum(thick ./ kval);
end
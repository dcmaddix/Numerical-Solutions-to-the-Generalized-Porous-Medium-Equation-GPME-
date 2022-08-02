function [f_plus, f_min] = computeEOsplitting(uL, uR, f, j)
%assumes we have flux function for Burgers 
%f = @(u) 1/2 *u^2;
f_plus = 0;
f_min = 0;
if (uL > 0)
    f_plus = f(uL, j);
end
if (uR < 0)
    f_min = f(uR, j);
end
end
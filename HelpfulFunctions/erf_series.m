function res = erf_series(z)
    invsqr = 1 ./ (z.^2);
    res = 1 - 0.5*invsqr + 0.75 * invsqr.^2 - (15 / 8) *invsqr.^3 ...
            + (105 / 16) *invsqr.^4; %last term on order 1e-6 check for conv analysis
end
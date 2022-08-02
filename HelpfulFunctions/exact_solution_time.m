%in limitiing case where eps -> 0
function u_t = exact_solution_time(t, x, alpha_star, A, B, eps, lim)
    %plot at final time
    psi = alpha_star * sqrt(t);
    if (x <= psi)
        u_t = 1 - A * erf(x(x <= psi) / (2 * sqrt(t)));
    elseif(~lim)
        z = x / (2*sqrt(eps*t));
        u_t = B * (exp(-z.^2) ./ (z*sqrt(pi))) .* erf_series(z);
    else
        u_t = 0;
    end
end
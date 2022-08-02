%in limitiing case where eps -> 0
function [u_t u] = exact_solution(t, x, ind, alpha_star, A, B, eps, lim)
    u = zeros(length(x),1);
    %plot at final time
    psi = alpha_star * sqrt(t);
    u(x <= psi) = 1 - A * erf(x(x <= psi) / (2 * sqrt(t)));
    %using expressions for 1-erf(z) = exp(-z^2) / (sqrt(pi)*z) * series
    %u2  = B * (1-erf(z))
    if (~lim) %otherwise in limit stay 0
        z = x(x > psi) / (2*sqrt(eps*t));
        u(x > psi) = B * (exp(-z.^2) ./ (z*sqrt(pi))) .* erf_series(z);
    end
    u_t = u(ind); %return scalr
end
    
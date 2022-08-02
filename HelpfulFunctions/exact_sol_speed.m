function [alpha_star A B] = exact_sol_speed(lim, theta, eps)
    %lim eps->0
    if (lim)
        f = @(alpha) (alpha / 2) .* erf(alpha / 2) .* exp((alpha / 2).^2) ...
            - (1 - theta)/(theta * sqrt(pi));
    else
        z = @(alpha) alpha / (2*sqrt(eps));
        f = @(alpha) theta * (alpha / 2) .* erf(alpha / 2) .* exp((alpha/2).^2) ...
                -( (1 - theta) / sqrt(pi) * erf_series(z(alpha)));
    end
    [alpha_star fval] = fzero(f, 1.2);
    if (fval > 1e-8)
        fprint('Opt Failed\n')
        exit(1)
    end
    A = (1 - theta) / erf(alpha_star / 2); %change to erf(alpha/(2*sqrt(k_max))
    B = 0;
    if (~lim)
        z = z(alpha_star); %evaluate at specific alpha star
        %B = theta / (1 - erf(z(alpha_star)))
        B = theta ./ (exp(-z.^2) ./ (z*sqrt(pi))) * erf_series(z);   
    end
end
function p_exact = computeExactIC(x_cent, lim, k_lower, t)
    cell = 1;
    %lim = true; %smooth initial condition
    %k_lower = 0.0; %0.01 for initial condition this is for plot
    p_star = 0.5;
    [alpha_star A_exact B_exact] = exact_sol_speed(lim, p_star, k_lower);
    [~, p_exact] = exact_solution(t, x_cent, cell, alpha_star,...
                                  A_exact, B_exact, k_lower, lim);
end
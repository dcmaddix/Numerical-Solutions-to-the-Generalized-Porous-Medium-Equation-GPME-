function computeIC(N)
%we take the exact solution at an evolved time hard coded as 0.0479 with
%k_lower = 0.01 rather than 0 so smooth and no shock
t = 0.0479;
lim = false;
p_star = 0.5;
k_lower = 0.01;
[alpha_star A_exact B_exact] = exact_sol_speed(lim, p_star, k_lower);
x_cent = 0:1/N:1;
[~, p_exact] = exact_solution(t, x_cent, 20, alpha_star,...
                         A_exact, B_exact, k_lower, lim);
%save the solution
save(sprintf('IC_0.0479/discont_%i.mat',N), 'p_exact')
end
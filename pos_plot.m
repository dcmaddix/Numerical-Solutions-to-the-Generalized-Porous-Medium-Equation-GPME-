gcf = 2;
x_exact = 0:1/3200:1;
if (k_pow == 1)
    p_exact = load('data/arith/p_1/arithmetic_pressure_refsol_p1');
elseif (k_pow == 2)
    p_exact = load('data/arith/p_2/arithmetic_pressure_refsol_p2');
elseif (k_pow == 3)
    p_exact = load('data/arith/p_3/arithmetic_pressure_refsol_p3');
else
    p_exact = load('discont');
end
if (k_pow == -1)
    x_exact = 0:1/100:1;
    p2(count+1) = plot(x_exact, p_exact.p_exact, '-o', 'Linewidth', 2);
    legend('show', 'N = 25', 'N = 50', 'N = 100', 'Exact. Sol')
else
    if (k_pow ~= 3)
        p2(count+1) = plot(x_exact, p_exact.p, 'Linewidth', 2);
    else
        p2(count+1) = plot(x_exact, p_exact.p_exact, 'Linewidth', 2);
    end
    legend('show', 'N = 25', 'N = 50', 'N = 100', 'Ref. Sol')
end
saveas(gcf,sprintf('~/Desktop/Results/p_%d/upwindFlux/%s_%s_pos.jpg', k_pow, JSolver, limiter))
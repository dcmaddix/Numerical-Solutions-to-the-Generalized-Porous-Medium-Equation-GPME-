    N = 100;
    dx = 1 / N;
    dt = dx^2/32;
    lim = true;
    t0 = 0.0479;
    nt = 0.1;
    t = t0 : dt : t0 + nt;
    k_lower = 0;
    p_star = 0.5;
    alpha_star = exact_sol_speed(lim, p_star, k_lower);
    x_star = alpha_star*t.^(0.5); %interface position at time t = 0
    v = 0.5 * alpha_star * t.^(-0.5);
    %calculate initial distance to the initial shock position on uniform
    %cartesian grid
    x = 0:dx:1;
    %initial signed dsitance function
    phi = x - x_star(1);
    for n = 1:length(t)-1
       phi(2:end) = phi(2:end) - dt*v(n)*(phi(2:end) - phi(1:end-1)) / dx;
    end
    %shock position is given by {x|phi(x,t) = 0}
    x_approx = interp1(phi, x, 0);
    fprintf('Approx shock loc = %f\n', x_approx);
    fprintf('Exact shock loc = %f\n', x_star(end));

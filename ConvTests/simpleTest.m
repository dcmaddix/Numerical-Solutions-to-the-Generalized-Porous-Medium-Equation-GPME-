function p= simpleTest(k_i, N, F0, alpha,nt, p0, delta_t, h, x_coord)    

    p = p0;
    F = F0; %initial source term
    
    %define domain bounds
    lb = 0.0; %for bdry heat eqtn test case
    ub = 1.0;%for bdry heat eqtn test case
    
    %define spatial vector at cell faces to plot against velocity
    x = linspace(lb + h/2, ub-h/2, N);
    
    %define x at cell-centers at cell centers to plot against pressure
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub];
    
    %arithmetic averaging
    k_A = @(k1,k2) 0.5 * (k1 + k2);
    
    %delta_t = nt /999;
    t = 0:delta_t:nt;
    %length(t)
    p_t = zeros(1,length(t)); %time vector where initial condition is given at t=0
    %compute constant lambda
    lambda = delta_t / h^2;
    k_H = @(k1,k2) (2*k1.*k2)./(k1+k2);
    %test case where k is constant in time
    k = zeros(size(p));
    k(x_cent <= 0.5) = 1000;
    k(x_cent > 0.5) = 1;
    %Does not depend on time can store at the beginning
   % k_plus = (k_A(k(2:end-1), k(3:end)));
    %k_minus = (k_A(k(1:end-2), k(2:end-1)));
    k_plus = (k_i(k(2:end-1), k(3:end)).^2 ./ k_A(k(2:end-1), k(3:end)));
    k_minus = (k_i(k(1:end-2), k(2:end-1)).^2 ./ k_A(k(2:end-1), k(1:end-2)));
    size(k_minus)
    size(p(2:end-1))
    %reorder for matrix off-diagonal
    k_minus_new = k_minus;
    k_minus_new(1:end-1) = k_minus(2:end);
    k_plus_new = k_plus;
    k_plus_new(2:end) = k_plus(1:end-1); %to get spdiag to work
%     D_k = spdiags([-k_minus_new, k_minus + k_plus, -k_plus_new], [-1 0 1], N-1, N-1);
%     M = alpha*speye(N-1,N-1) + lambda*D_k;
%     F(1) = F(1) + k_minus(1)*p(1) / h^2;
%     F(end) = F(end) + k_plus(end)*p(end) / h^2;
%     E = M \ speye(N-1,N-1);
    p_t = zeros(1,length(t));
     ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6);
    if (length(ind) > 1)
        ind = ind(1);
    elseif (isempty(ind)) %make sure xcoord exists
        ind = 1;
        fprintf('Ind warning\n')
    end
    %boundaries also constant in time
    for i = 1:length(t)-1
        %write simple finite volume solver and separate BC    
        %this makes it explicit change to Newton Raphson for fully-implicit
        
        p(2:end-1) = p(2:end-1) + lambda * (k_minus.*p(1:end-2) - (k_minus + k_plus).*p(2:end-1)...
                    + k_plus.*p(3:end));
        
        %%%%%%POWERS OD P%%%%
        %k = (p.^k_pow)';%update at each timestep depends on new p depends on which timestepping scheme
        %do taylor linearization fo matrix-can't change for BE/CN
        %important how update the nonlinear k in time
        p_t(i) = p(ind);
        %p(2:end-1) = E* (p(2:end-1) + delta_t * F); %BAckwardEuler
    end
    p_t(end) = p(ind);
    figure(4);
    plot(t,p_t)
    hold on;
end

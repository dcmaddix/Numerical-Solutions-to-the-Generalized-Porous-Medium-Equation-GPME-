%linear convergence study
n_grids = 10; %use 10 for 0.0001 7:end/10 for .1/.01 AND 12  for 0.001 with 3:end in symptotic conve rate
N = zeros(n_grids,1);
N(end) = 100*2^(n_grids-1);
h = zeros(n_grids,1);
h(end) = 1 / N(end);
bdry = zeros(2,1);
k_lb = 0.01;
lb = 0;
ub = 1;
x = linspace(h(end)/2,1-h(end)/2, N(end)); %to plot against velocity on finest grid
x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures
p_exact =( -(1 / (log(k_lb) * (1-k_lb))) * log((1-k_lb)*x_cent + k_lb) -...
         x_cent / (1-k_lb)  + 1 / (1-k_lb))';
%never actually comparing on the finest     
u_exact = (-(-x -k_lb/(1-k_lb) - 1/log(k_lb)) ./ ((1-k_lb)*x+k_lb))';
error_p = zeros(3,n_grids); %rows is methods
error_u = zeros(3,n_grids); %rows is methods
for i = 0:n_grids - 1
    N(i+1) = 100*2^i;
    h(i+1) = 1 / N(i+1);
    k = linspace(k_lb,1,N(i+1)+1);
    F = ones(N(i+1)-1,1);
    [p_harm p_arith p_FV, ~, u_harm u_arith u_FV] = testSchemeI(N(i+1),k, bdry, F);
    %computed exact solution on finest grid and indexing it to match at x
    %points on coarser grid. Could also compute exact solution on all the
    %grids as done for u
    error_p(1,i+1) = max(abs(p_harm - p_exact(1:2^(n_grids-1-i):end))) /max(abs(p_harm));
    error_p(2,i+1) = max(abs(p_arith - p_exact(1:2^(n_grids -1 -i):end))) / max(abs(p_arith));
    error_p(3,i+1) = max(abs(p_FV - p_exact(1:2^(n_grids-1-i):end))) / max(abs(p_FV));
%      error_p(1,i+1) = h(i+1)^(0.5)*norm(p_harm - p_exact(1:2^(n_grids-1-i):end));
%      error_p(2,i+1) = h(i+1)^(0.5)*norm(p_arith - p_exact(1:2^(n_grids -1 -i):end));
%      error_p(3,i+1) = h(i+1)^(0.5)*norm(p_FV - p_exact(1:2^(n_grids-1-i):end));
%     max(abs(p_FV));
    %relative error
    %calculate exact u on the grid rather than the spacing by 2
    x = linspace(h(i+1)/2,1-h(i+1)/2, N(i+1)); %to plot against velocity on finest grid
    u_exact = (-(-x -k_lb/(1-k_lb) - 1/log(k_lb)) ./ ((1-k_lb)*x+k_lb))';
    error_u(1,i+1) = max(abs(u_harm - u_exact)) / max(abs(u_harm));
    error_u(2,i+1) = max(abs(u_arith - u_exact)) / max(abs(u_arith));
    error_u(3,i+1) = max(abs(u_FV - u_exact))/ max(abs(u_FV)) ;
%     error_u(1,i+1) = h(i+1)^(0.5)*norm(u_harm - u_exact);
%     error_u(2,i+1) = h(i+1)^(0.5)*norm(u_arith - u_exact);
%     error_u(3,i+1) = h(i+1)^(0.5)*norm(u_FV - u_exact);
    %second order convergence in the maximum norm and l2 norm! perfect
    %linear loglog plot in asymtotic convergence region
end
error_p
error_u
convRate_p = zeros(1,3);
convRate_u = zeros(1,3);
for j = 1:3
    p = polyfit(log10(h(1:end)), log10(error_p(j,1:end))',1);
    u = polyfit(log10(h(1:end)), log10(error_u(j,1:end))',1);
    convRate_p(j) = p(1);
    convRate_u(j) = u(1);
end
convRate_p
convRate_u
figure(1);
loglog(h,error_p(1,:),h, error_p(2,:),h, error_p(3,:))
xlabel('Spatial step h')
ylabel('Log of Error')
title('Loglog plot of pressure error as measured in the maximum norm for the various schemes')
legend('MFD harm','arith','FV harm')
figure(2);
loglog(h,error_u(1,:), h, error_u(2,:), h, error_u(3,:))
xlabel('Spatial step h')
ylabel('Log of Error')
title('Loglog plot of velocity error as measured in the maximum norm for the various schemes')
legend('MFD harm','arith','FV harm')

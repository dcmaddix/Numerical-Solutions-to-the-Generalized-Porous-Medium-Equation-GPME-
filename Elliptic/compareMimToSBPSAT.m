function [error_p error_u] = compareMimToSBPSAT(N, seed) %input fine grid size k defined via coarse grid size
%this script compares the mimetic schemes to SBP-SAT where k is randomly
%generated and so represents heterogrenous coefficients
h = 1 / N;
%to define k assuming that N >= 100
%compute k so same with grid refinement
N_coarse = 100; %coarest k grid is 100 refine from there
h_coarse = 1 ./ N_coarse;
rand('seed',seed)%2 has big jump from 0 to 1, 19, 20 13
k = zeros(N-1,1)';% more interesting discontinuites and heterogenities than
%analytic piecewise constant sol
x_coarse = linspace(h_coarse/2,1-h_coarse/2, N_coarse);
x_cent_coarse = 0.5 * (x_coarse(1:end-1) + x_coarse(2:end));
x = linspace(h/2,1-h/2, N); %to plot against velocity
%define x at cell-centers
x_cent = 0.5 * (x(1:end-1) + x(2:end)); %to plot against pressures
%loop on interior
for i = 1:N_coarse-2 %want k constant on cells of fine grid sicne length of x_cent_coarse is N_coarse-1
    k(x_cent >= x_cent_coarse(i) & x_cent <= x_cent_coarse(i+1)) = rand(1); %k and x_cent have the same indices
end
%do left and right boundary
k(x_cent <= x_cent_coarse(1)) = rand(1);
k(x_cent >= x_cent_coarse(end)) = rand(1);
%add in k exactly zero as test
k = [k(1) k k(end)]; % extend variable coefficients to boundary grid
%points - simply by extrapolating
%test non-homogenous Dir BC
F = zeros(N-1,1); %0 RHS add non-homogenous Dir condition to obtain nontrivial solution
bdry = [1 2]; %Dirichlet bdry values fixed g function
lb = 0; %physical lower bound 
ub = 1; %physical upper bound
x_cent = [lb x_cent ub];
[p_harm p_arith p_FV p_SAT u_harm u_arith u_FV] = testSchemeI(N,k, bdry, F);
%compute p_exact and u_exact = -dp/dx
p_exact = zeros(N+1,1);
p_exact(1) = bdry(1);
%constant of integration holds for the Dirichlet boundary case to enforce
%right bdry
%unfold recursion itnergating twice get p_{N} = C / K_{N-1}(XN-XN-1) +
%p_N-1 = bdry(2) sub in p_N-1 unlike get to p(1) = bry(1)
C = (bdry(2)-bdry(1)) / sum((x_cent(2:end) - x_cent(1:end-1)) ./ k(1:end-1));
%define constant of integration to satisfy right bc
for i = 1:length(x_cent)-1
  p_exact(i+1) = C*(x_cent(i+1)-x_cent(i))/k(i) + p_exact(i);
  %u is derivative
end
%exact u solution
u_exact = (-C ./ k(1:end-1))';

%compute errors
error_p = zeros(4,1);
error_u = zeros(3,1); %need to define u for SBP-SAT

error_p(1) = h^0.5*norm(p_exact - p_harm);
error_p(2) = h^0.5*norm(p_exact - p_arith);
error_p(3) = h^0.5*norm(p_exact - p_FV);
error_p(4) = h^0.5*norm(p_exact - p_SAT(1:end-1));

error_u(1) = h^0.5*norm(u_exact - u_harm);
error_u(2) = h^0.5*norm(u_exact- u_arith);
error_u(3) = h^0.5*norm(u_exact - u_FV);


%max norm error:velocity convergence goes flat in max norm.
%location occurs because u is piecewise constant on cell and at
%discontinuities have sharp transition in exact but in numerical at cell
%enpts have a weighted average value. Numerical differs at the jumps
%showing lack of convergence in u
%On final grid for n_grids = 10, exact solution is constant for 4 grid
%cells. Numerical solution is matching at 3 of those 4 locations and has a
%weighted average value at the jumps!
% error_p(1) = max(abs(p_exact - p_harm)) / max(abs(p_harm));
% error_p(2) = max(abs(p_exact - p_arith)) / max(abs(p_arith));
% error_p(3) = max(abs(p_exact - p_FV)) / max(abs(p_FV));
% error_p(4) = max(abs(p_exact - p_SAT(1:end-1))) / max(abs(p_SAT(1:end-1)));

% [max1 ind1] = max(abs(u_exact - u_harm));
% [max2 ind2] = max(abs(u_exact- u_arith));
% [max3 ind3] = max(abs(u_exact - u_FV));
% 
% error_u(1) = max1 / max(abs(u_harm));
% error_u(2) = max2 / max(abs(u_arith));
% error_u(3) = max3 / max(abs(u_FV));

%plot results and corresponding k
figure(2);
plot(x_cent,p_harm, x_cent, p_arith, x_cent, p_FV, x_cent, p_SAT(1:end-1), x_cent, p_exact)
title('Plot of pressure from Mimetic Scheme I with various k_i versus the exact solution for random heterogenous k')
legend('MFD Harmonic', 'Arithmetic', 'FV Harmonic', 'SBP-SAT', 'Exact Solution')
xlabel('x')
ylabel('pressure p')
axis([0 1 1 2])

figure(3);
plot(x_cent,k)
xlabel('x')
ylabel('permeability k')
title('Plot of permeabilities')

figure(4)
plot(x, u_harm, x, u_arith, x, u_FV, x, u_exact)
legend('MFD Harmonic', 'Arithmetic', 'FV Harmonic', 'Exact Solution')
xlabel('x')
ylabel('velocity u')
title('Plot of velocity from Mimetic Scheme I with various k_i versus the exact solution for random heterogenous k')

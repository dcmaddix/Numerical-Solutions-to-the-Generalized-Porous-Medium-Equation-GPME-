% script to compute pressure for 1D Poisson equation
n = 202;
h = 1;
% permability field - change to realistic data (cell values)
c = rand(n-2,1);

% extend variable coefficients to boundary grid points - simply by extrapolating
c = [c(1); c; c(end)];

% pressure boundary conditions
p1 = 1;
p2 = 2;
g = [p1 p2];

% Boundary condition type: D for Dirichlet, N for Neumann
g_type = 'D';
%g_type = 'N';

% Flux  type - 'H' leads to TPFA scheme in interior (and better suited than 
% 'A' for heterogeneous coefficients but 'A' is common for SBP-SAT methods)
%type = 'A';
type = 'H';
%harmonic vs arithmetic averaging?

% Imposition of boundary conditions - 'weak' is SAT conditions, 'strong' is
% exact imposition (via SAT terms)
bc = 'weak';
%bc = 'strong';
%depends on monotonicity guarantee


% compute the variable coefficient 2nd derivative operator D with BCs and RHS F
[D, F, H, HI] = set_up_1D_matrix_SBPSAT(n,h,c,bc,g_type,g,type);

% solve for pressure
P = D\F;

%use same c and see if get the same solution
figure
plot(P);









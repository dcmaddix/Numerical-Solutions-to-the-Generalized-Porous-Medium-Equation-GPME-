%This function computes a high order wide stencil first derivative to be
%used in truncation error cancellation
function dp_dx = central_4thOrder2D(p, ghostNodes, dx,xderiv)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h

N_int  = length(p(2:end-1,1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
dp_dx = zeros(N_int,N_int);
%vector of interior gridpoints for wide stencil
i = 3 : N - 2;
j = 2:N-1;
%Fourth Order First Derivative
if (xderiv)
    dp_dx(1,:) = (-p(4,j) + 8 * p(3,j) - 8 * p(1,j) + ghostNodes{1}) / (12 * dx);
    dp_dx(end,:) = (-ghostNodes{2} + 8 * p(N,j) - 8 * p(N - 2,j) + p(N - 3,j) ) / (12 * dx);
    dp_dx(2:end-1,:) = (-p(i + 2,j) + 8 * p(i + 1,j) - 8 * p(i - 1,j) + p(i - 2,j) ) / (12 * dx);
else
    dp_dx(:,1) = (-p(j,4) + 8 * p(j,3) - 8 * p(j,1) + ghostNodes{1} ) / (12 * dx);
    dp_dx(:,end) = (-ghostNodes{2} + 8 * p(j,N) - 8 * p(j,N - 2) + p(j,N - 3) ) / (12 * dx);
    dp_dx(:,2:end-1) = (-p(j,i + 2) + 8 * p(j,i + 1) - 8 * p(j,i - 1) + p(j,i - 2) ) / (12 * dx);
end
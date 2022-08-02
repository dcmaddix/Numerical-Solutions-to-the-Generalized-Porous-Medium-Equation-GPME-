%This function computes a high order wide stencil first derivative to be
%used in truncation error cancellation
function dp_dx = central_4thOrder(p, ghostNodes, dx)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h

N_int  = length(p(2:end-1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
dp_dx = zeros(N_int,1);
%vector of interior gridpoints for wide stencil
i = 3 : N - 2;

%Fourth Order First Derivative
dp_dx(1) = (-p(4) + 8 * p(3) - 8 * p(1) + ghostNodes(1) ) / (12 * dx);
dp_dx(end) = (-ghostNodes(2) + 8 * p(N) - 8 * p(N - 2) + p(N - 3) ) / (12 * dx);
dp_dx(2:end-1) = (-p(i + 2) + 8 * p(i + 1) - 8 * p(i - 1) + p(i - 2) ) / (12 * dx);
end
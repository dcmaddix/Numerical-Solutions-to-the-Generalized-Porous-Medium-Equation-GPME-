%This function computes a high order wide stencil second derivative to be
%used in truncation error cancellation
function d3p_dx3 = central_3rdDeriv_2ndOrder(p, ghostNodes, dx)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h
N_int  = length(p(2:end-1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
d3p_dx3 = zeros(N_int,1);
%vector of interior gridpoints for wide stencil
i = 3 : N - 2;

%Fourth Order Second Derivative  
%left interior point i = 2
d3p_dx3(1) = (p(4) - 2 * p(3) + 2 * p(1) - ghostNodes(1) ) / (2 * dx^3);
%right interior point i = N - 1
d3p_dx3(end) = (ghostNodes(2) - 2 * p(N) + 2 * p(N - 2) - p(N - 3) ) / (2 * dx^3);
%i = 3:N-2
d3p_dx3(2:end-1) = (p(i + 2) - 2 * p(i + 1) + 2 * p(i - 1) - p(i - 2) ) / (2 * dx^3);
end
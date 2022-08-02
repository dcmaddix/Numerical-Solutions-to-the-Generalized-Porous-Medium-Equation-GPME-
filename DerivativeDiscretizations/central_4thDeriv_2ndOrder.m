%This function computes a high order wide stencil second derivative to be
%used in truncation error cancellation
function d4p_dx4 = central_4thDeriv_2ndOrder(p, ghostNodes, dx)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h
N_int  = length(p(2:end-1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
d4p_dx4 = zeros(N_int,1);
%vector of interior gridpoints for wide stencil
i = 3 : N - 2;

%Fourth Order Second Derivative  
%left interior point i = 2
d4p_dx4(1) = (p(4) - 4 * p(3) + 6*p(2) - 4 * p(1) + ghostNodes(1) ) / (dx^4);
%right interior point i = N - 1
d4p_dx4(end) = (ghostNodes(2) - 4 * p(N) + 6 * p(N-1) - 4* p(N - 2) + p(N - 3) ) / (dx^4);
%i = 3:N-2
d4p_dx4(2:end-1) = (p(i + 2) - 4 * p(i + 1) + 6*p(i) -4 * p(i - 1) + p(i - 2) ) / (dx^4);
end
%This function computes a high order wide stencil second derivative to be
%used in truncation error cancellation
function d3p_dx3 = central_3rdDeriv_4thOrder(p, ghostNodes, dx)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h
%ghostNodes(3) is left ghost node at lb - 2h 
%ghostNodes(4) is right ghost node at ub + 2h
%increased length of array
N_int  = length(p(2:end-1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
d3p_dx3 = zeros(N_int,1);
%vector of interior gridpoints for wide stencil
i = 4 : N - 3;
%Fourth Order Second Derivative  
%left interior point i = 2
d3p_dx3(1) = (-p(5) + 8 * p(4) - 13 * p(3) + 13 * p(1) - 8 * ghostNodes(1)...
            +ghostNodes(3)) / (8 * dx^3);
%i = 3
d3p_dx3(2) = (-p(6) + 8 * p(5) - 13 * p(4) + 13 * p(2) - 8 * p(1) ...
            +ghostNodes(1)) / (8 * dx^3);  
%right interior point i = N - 1
d3p_dx3(end) = (-ghostNodes(4) + 8 * ghostNodes(2) - 13 * p(N) + 13 * p(N - 2) ...
               - 8 * p(N - 3) +  p(N - 4) ) / (8 * dx^3);
%right intertior point i = N-2
d3p_dx3(end-1) = (-ghostNodes(2) + 8 * p(N) - 13 * p(N-1) + 13 * p(N - 3) ...
               - 8 * p(N - 4) +  p(N - 5) ) / (8 * dx^3);
%i = 3:N-2
d3p_dx3(3:end-2) = (-p(i + 3) + 8 * p(i + 2) - 13 * p(i + 1) + 13 * p(i - 1) ...
                    - 8 * p(i - 2)  + p(i - 3)) / (8 * dx^3);
end
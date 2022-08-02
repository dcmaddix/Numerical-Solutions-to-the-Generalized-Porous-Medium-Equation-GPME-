%This function computes a high order wide stencil second derivative to be
%used in truncation error cancellation
function d4p_dx4 = central_4thDeriv_4thOrder(p, ghostNodes, dx)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h
%ghostNodes(3) is left ghost node at lb - 2h 
%ghostNodes(4) is right ghost node at ub + 2h
%increased length of array
N_int  = length(p(2:end-1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
d4p_dx4 = zeros(N_int,1);
%vector of interior gridpoints for wide stencil
i = 4 : N - 3;

%Fourth Order Second Derivative  
%left interior point i = 2
%Usign Dir BC and assumign ghostnodes equal to right bdry this scheme
%requires 2 ghostnodes
d4p_dx4(1) = (-p(5) + 12 * p(4) - 39 * p(3) +56*p(2) - 39 * p(1) + 12*ghostNodes(1)...
            -ghostNodes(3)) / (6*dx^4);
%i = 3
d4p_dx4(2) = (-p(6) + 12 * p(5) - 39 * p(4) +56*p(3) - 39 * p(2) + 12*p(1)...
            -ghostNodes(1)) / (6*dx^4);
%right interior point i = N - 1
d4p_dx4(end) = (-ghostNodes(4) + 12*ghostNodes(2) - 39 * p(N) + 56 * p(N-1) ...
                - 39* p(N - 2) + 12 * p(N - 3) - p(N - 4)) / (6 * dx^4);
%right intertior point i = N-2
d4p_dx4(end-1) = (-ghostNodes(2) + 12*p(N) - 39 * p(N - 1) + 56 * p(N - 2) ...
                - 39* p(N - 3) + 12 * p(N - 4) - p(N - 5)) / (6 * dx^4);
%i = 4:N-3
d4p_dx4(3:end-2) = (-p(i + 3) + 12 * p(i + 2) - 39 * p(i + 1) + 56*p(i) ...
                   -39 * p(i - 1) + 12 * p(i - 2)  - p(i - 3)) / (6 * dx^4);
end
%This function computes a high order wide stencil second derivative to be
%used in truncation error cancellation
function d2p_dx2 = central_2ndDeriv_4thOrder2D(p, ghostNodes, dx, xderiv)
%ghostNodes(1) is left ghost node at lb - h 
%ghostNodes(2) is right ghost node at ub + h
N_int  = length(p(2:end-1,1));
%Total gridpoints is interior plus two boundaries
N = N_int + 2;
%create derivative vectors
d2p_dx2 = zeros(N_int,N_int);
%vector of interior gridpoints for wide stencil
i = 3 : N - 2;
j = 2:N-1;
%Fourth Order Second Derivative  
if (xderiv)
    d2p_dx2(1,:) = (-p(4,j) + 16 * p(3,j) - 30 * p(2,j) + 16 * p(1,j) - ghostNodes{1} ) ...
                / (12 * dx^2);
    d2p_dx2(end,:) = (-ghostNodes{2} + 16 * p(N,j) - 30 * p(N - 1,j) + 16 * p(N - 2,j)...
                    - p(N - 3,j) ) / (12 * dx^2);
    d2p_dx2(2:end-1,:) = (-p(i + 2,j) + 16 * p(i + 1,j) - 30 * p(i,j) + 16 * p(i - 1,j) ...
                        - p(i - 2,j) ) / (12 * dx^2);
else
    d2p_dx2(:,1) = (-p(j,4) + 16 * p(j,3) - 30 * p(j,2) + 16 * p(j,1) - ghostNodes{1} ) ...
                / (12 * dx^2);
    d2p_dx2(:,end) = (-ghostNodes{2} + 16 * p(j,N) - 30 * p(j,N - 1) + 16 * p(j,N - 2)...
                    - p(j,N - 3) ) / (12 * dx^2);
    d2p_dx2(:,2:end-1) = (-p(j,i + 2) + 16 * p(j,i + 1) - 30 * p(j,i) + 16 * p(j,i - 1) ...
                        - p(j,i - 2) ) / (12 * dx^2); 
end
end
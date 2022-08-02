%finite difference pressure solve p(0) = p(end) = 0
%d^@p/dx^2 =1 poisson eqtn
N = 10; %number of gridpoints including the two boundary
h = 1 / (N-1);
p = zeros(N,1);
x = zeros(N,1);
for i = 1:N
    x(i) = h*(i-1);
end
e = ones(N-1,1);
D = spdiags([-e, 2*e, -e], [-1 0 1], N-2, N-2);
F = h^2*ones(N-2,1);
p(2:end-1) = D\F;
p_exact = 0.5*(-x.^2 + x); %only in constant coeff case
err = max(abs(p - p_exact))
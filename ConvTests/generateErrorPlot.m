%this script creates a loglog plot of the error versus timestep size for
%the simple test case d^2p/dx^2 =1
n_runs = 3;
N = [100 200 400];
h = 1 ./ N; %uniform grid spacing
err_u = zeros(n_runs,1);
err_p = zeros(n_runs,1);
for i = 1:n_runs
    err_u(i) = testConstantK(N(i));
    err_p(i) = testConstantK(N(i));
end
loglog(h, err_p, h, err_u);
%compute slope of line
press_slope = polyfit(log10(h),log10(err_p),1);
press_slope = press_slope(1)
vel_slope = polyfit(log10(h),log10(err_u),1);
vel_slope = vel_slope(1)

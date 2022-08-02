%Convergence study which calls compareMimToSBPSAT script
seed = 1; %fix randomly generated seed for k
%compute exact solution FV with harmonic on finest grid 3200
n_grids = 10; %12
N = zeros(n_grids,1);
h = zeros(n_grids,1);
error_p = zeros(4,n_grids); % each row corresponds to different averaging
error_u = zeros(3,n_grids);
for i = 0:n_grids - 1
    N(i+1) = 100 *2^i; 
    h(i+1) = 1 / N(i+1);
    [error_pi error_ui] = compareMimToSBPSAT(N(i+1), seed);
    %need to project p_exact from fine grid to coarse grid see assignment
    %ps1c_conv
    error_p(:,i+1) = error_pi;
    error_u(:,i+1) = error_ui;    
end
error_p
error_u
convRate_p = zeros(1,4);
convRate_u = zeros(1,3);
%rate of convergence for each method and loglog error plots as measured in
%the maximum norm, MFD with harmonic, arithmetic, FV with harmonic, SBp
%with harmonic
for i = 1:4 %h(4) for asymptotic region
    p = polyfit(log10(h(5:end)), log10(error_p(i,5:end))',1); %took image of 3:end in asymptotic convergence region
    convRate_p(i) = p(1);%4:end for convergence plot seed 0
    if (i <= 3)
        p = polyfit(log10(h(4:end)), log10(error_u(i,4:end))',1);
        convRate_u(i) = p(1);
    end
end
convRate_p
convRate_u
figure(3);
loglog(h,error_p(1,:), h, error_p(2,:), h, error_p(3,:), h, error_p(4,:))
xlabel('Spatial step h')
ylabel('Log of Error')
title('Loglog plot of pressure error as measured in l_2-norm for the various schemes')
legend('MFD harm','arith','FV harm', 'SBP-SAT')
figure(10);
loglog(h,error_u(1,:), h, error_u(2,:), h, error_u(3,:))
xlabel('Spatial step h')
ylabel('Log of Error')
title('Loglog plot of velocity error as measured in l_2 norm for the various schemes')
legend('MFD harm','arith','FV harm')
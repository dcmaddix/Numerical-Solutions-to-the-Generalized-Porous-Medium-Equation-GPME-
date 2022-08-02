
%convergence study for parabolic problem for k = p^3 and arithmetic
%averaging
N = [25 50 100 200];
n_grids = length(N);
h = 1 ./ N;
l2_error = zeros(3, length(N));
dir_int = 'data/foam/int_nx';
dir_arith = 'data/foam/arith_nx';
dir_sam = 'data/foam/sam_nx';
dir_exact = 'data/foam/exact_nx';
lb = 0.0; %for bdry heat eqtn test case
ub = 1.0; %for bdry heat eqtn test case
normtype = 2;
for i = 1:n_grids
    %store harmonic error
    num = load(strcat(dir_arith, sprintf('%i', N(i))));
    exact = load(strcat(dir_exact, sprintf('%i', N(i))));
    l2_error(1,i) = norm(num.p - exact.p_exact, normtype);
    %store arith error
    num = load(strcat(dir_int, sprintf('%i', N(i))));
    l2_error(2,i) = norm(num.p - exact.p_exact, normtype);
    num = load(strcat(dir_sam, sprintf('%i', N(i))));
    l2_error(3,i) = norm(num.p - exact.p_exact,normtype);
    if (normtype ~= Inf)
        l2_error(:,i) = l2_error(:,i) *h(i)^(0.5);
    end
end
%legend(plots, 'N = 25', 'N = 50', 'N = 100', 'N = 200', 'N = 400', 'N = 800')
conv_arith_p3 = polyfit(log10(h), log10(l2_error(1,:)),1);
conv_arith_p3(1)
conv_int_p3 = polyfit(log10(h), log10(l2_error(2,:)),1);
conv_int_p3(1)
conv_sam_p3 = polyfit(log10(h), log10(l2_error(3,:)),1);
conv_sam_p3(1)
figure(2)
p = plot(log10(h),log10(l2_error(1,:)), '-rs',...,
      log10(h), log10(l2_error(2,:)),'-gd', ...
      log10(h), log10(l2_error(3,:)), '-o');
set(p(2), 'Color',[0 0.498039215803146 0]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
legend('Arithmetic', 'Integral', 'SAM', 'Location', 'southeast')
xlabel('log_{10}(\Deltax)')
ylabel('log_{10}(Error)')
title('l_{2} norm of error')
print(2,'-depsc2','test');
fixPSlinestyle('test.eps','~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/l2_error.eps');

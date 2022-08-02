%plot arithmetic average in between timesteps
%close all
N = 50;
halfdt = true;
one = false;
save = true;
addpath('fixPSlinestyle')
if (halfdt)
    %load(sprintf('kt_arith%idthalf.mat',N)
    if (one)
        p = load(sprintf('data/foam/arith_nx%idthalf.mat',N));
        load(sprintf('kt_arith%idthalf.mat',N));
        load('ktplus_arith50plusdthalf.mat');
    else
        p = load('ktplus_arith50plusdthalf.mat');
    end
else
    if (one)
        load(sprintf('kt_arith%i.mat',N))
        load(sprintf('ktplus_arith%iplus.mat',N))
    else
        load(sprintf('ktplus_arith%iplus.mat',N))
    end
    if (one)
        p = load(sprintf('data/foam/arith_nx%i.mat',N));
    else
        p = load('p_tplus_arith.mat');
    end
end
k = ones(size(p.p_t));
k(p.p_t < 0.5) = 1/2;
close all
t = p.t - 0.0479;
p_ub = 0.48;
if (N == 50)
    p_ub = 0.47;
end
[ax h1 h2] = plotyy(t, p.p_t,t, p.k_t);
ub = 0.0698-0.0479;
if (N == 100)
    ub = 0.06935-0.0479;
end
%if (one)
    set(ax,'XLim',[0.0214 ub]);
    
%else
%    set(ax,'XLim',[0.015 0.0155]);
%end
%set(ax,'XLim',[min(t_shock-0.0479), max(t_shock-0.0479)])
if (one)
    set(ax(1),'YLim',[0.49 0.51]);
else
   %set(ax(1),'YLim',[0.07 0.12]);
   set(ax(1),'YLim',[0.0 0.03]);
end
if (one)
    set(ax(2),'YLim',[0 0.5]);
else
    set(ax(2),'YLim',[0 0.5]);
end
set(ax, 'FontSize',16,'FontWeight','bold');
set(ax,{'ycolor'},{'r';[0.0392156876623631 0.141176477074623 0.415686279535294]})
set(h1, 'Marker','o');
set(h1, 'Color','r');
set(h2, 'Color',[0.0392156876623631 0.141176477074623 0.415686279535294]);
set(h2,'Marker','o');
hold on
set(findall(ax, 'Type', 'Line'),'LineWidth',2);
plot(t, 0.5*ones(size(p.t)), '--k', 'LineWidth',4)
xlabel('t')
if (one)
    set(ax(1), 'YTick', [0.49 0.5 0.51])
    set(ax(2), 'YTick', [0 0.5])
else
    %set(ax(1), 'YTick', [0.07 0.1 0.12])
    set(ax(1), 'YTick', [0.0 0.01 0.02 0.03])
    set(ax(2), 'YTick', [0 0.5])
end
if (N == 100)
    set(ax,'XTick',[0.0213 0.02144])
    if (N == 100)
        set(ax(1), 'YTick', [0.45 0.5 0.515])
    end
end
if (one)
    ylabel(ax(1), 'p')
    ylabel(ax(2), 'k^A_{i+3/2}')
else
    ylabel(ax(1), 'p')
    ylabel(ax(2), 'k^A_{i+3/2}')
end
title('Solution profile and arithmetic average as a function of time')
print(1,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/Discontinous/Results/';
if (save)
    if (one && ~halfdt)
        filename = sprintf('arith_oscill%i',N);
    elseif(one && halfdt)
        filename = sprintf('arith_oscill%idthalf',N);
    elseif(~one && ~halfdt)
        filename = sprintf('arith_oscill%ileft.eps',N);
    else
        filename = sprintf('arith_oscill%ileftdthalf.eps',N);
    end
    fixPSlinestyle('test.eps',[fileLoc filename]);
end
%figure(2)
% close all
% ka = load('kt_arith100.mat');
% ki = load('kt_int100.mat');
% kh = load('kt_harm100.mat');
% ks_plus = load('kt_sam100plus.mat');
% ks_min = load('kt_sam100minus.mat');
% plot(ka.k_t(1:end-1), 'o')%,t(1:end-1), ki.k_t(1:end-1), '*',...
%     %t(1:end-1), kh.k_t(1:end-1),'x', t(1:end-1), ks_plus.k_t(1:end-1),'-x', ...
%     %t(1:end-1), ks_min.k_t(1:end-1),'-x');
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% %set(gca,'XLim',[0.0186 0.02]);
% legend('Arithmetic', 'Integral','Harmonic', 'SAM F_i^+', 'SAM F_{i+1}^-')
% %close all
% p = load(sprintf('data/foam/int_nx%i.mat',N));
% %[ax h1 h2] = plotyy(t, p.p_t,t, ki.k_t);
%set(h1, 'Marker','o');
%set(h2, 'Marker','o');
%set(findall(ax, 'Type', 'Line'),'LineWidth',2);
%legend('p_t', 'SAM F_i^+', 'SAM F_{i+1}^-')
%close all
% plot(t(1:end-1), ki.k_t(1:end-1))
% set(gca,'XLim',[min(ki.t_shock)-0.0479 max(ki.t_shock)-0.0479]);
% figure;
% plot(t(1:end-1), ka.k_t(1:end-1),'o', t(1:end-1), kh.k_t(1:end-1),'x')
% set(gca,'XLim',[0.01867 max(ka.t_shock)-0.0479]);
% plot(t(1:end-1), ki.k_t(1:end-1),'-o', t(1:end-1), ks_plus.k_t(1:end-1),'x')
% set(gca,'XLim',[min(ks_plus.t_shock)-0.0479 max(ki.t_shock)-0.0479]);
%set(gca,'YLim',[0.45 0.6]);



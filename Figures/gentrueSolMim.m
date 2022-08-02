%This function plots the true solution to the test problem
function [x_cent p_exact] = gentrueSolMim(N) %Run for various times
    %Solution as given in Section 4.5 p. 124 of the MFD with staggered
    %discretization of diffusion coeff
    h = 1/ N;
    eps = 1e-9; %approach 0 to get true analytical solution
    c = 1.0; %constant in analytical formula
    %define domain bounds
    lb = 0.0;
    ub = 1.0;
    x = linspace(h/2,1-h/2, N); %to plot against velocity
    %define x at cell-centers
    x_cent = [lb 0.5 * (x(1:end-1) + x(2:end)) ub]; %to plot against pressures
    %search for indices less than ct
    t = [0.3 0.5 0.7]; %define times to plot
    k = zeros(N+1,1);
    p_exact = zeros(3,N+1);
    %compute at x_cent = 0.12
    for i = 1:length(t)
        k(x_cent < c*t(i)) = 3*c*(c*t(i)-x_cent(x_cent < c*t(i)));
        k(x_cent >= c*t(i)) = eps;
        p_exact(i,:) = k.^(1/3);
    end
    %exact solution plot
    figure(11);
    plot(x_cent,p_exact(1,:), '*', x_cent,p_exact(2,:),'*', x_cent, p_exact(3,:), '*')
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'FontSize',16,'FontWeight','bold');
    title('Plot of exact solution at various times')
    legend('t = 0.3', 't = 0.5', 't = 0.7')
    xlabel('x')
    ylabel('p') 
    fileLoc = 'data/p^m/mimeticp3/';
    ind = 1;
    for i = 3:2:7
        pa = load([fileLoc sprintf('p50_arith%i.mat',i)]); 
        pmh = load([fileLoc sprintf('p50_mh%i.mat',i)]);
        plot1 = plot(x_cent, p_exact(ind,:),'-k', pa.x_cent, pa.p, ...
         pmh.x_cent, pmh.p, 'MarkerSize',8,'LineWidth',2);
        set(plot1(2),'Marker','o','LineStyle','-.','Color',[1 0 0],...
        'DisplayName','Arithmetic');
        set(plot1(3),'Marker','square','LineStyle',':',...
        'Color',[0 0.498039215803146 0],...
        'DisplayName','Mod. Harmonic');
        legend('Exact', 'Arithmetic', 'Mod. Harm')
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca, 'FontSize',16,'FontWeight','bold');
        hold on
        ind = ind + 1;
    end
    xlabel('x')
    ylabel('p') 
    title('Pressure vs. position at various times')
    dt = 2.5e-5;
    t = 0:dt:0.7;
    x_coord = 0.12;
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6, 1, 'first');
    pt_exact = zeros(1, length(t));
    for i = 1:length(t)
        k(x_cent < c*t(i)) = 3*c*(c*t(i)-x_cent(x_cent < c*t(i)));
        k(x_cent >= c*t(i)) = eps;
        pt_exact(i) = k(ind).^(1/3);
        pt_exact(i);
    end
    figure(3)
    pat = load([fileLoc 'p50_arith7t.mat']); 
    pmht = load([fileLoc 'p50_mh7t.mat']);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    plot1 = plot(t, pt_exact, '-k', pat.t, pat.p_t, pmht.t, pmht.p_t, 'LineWidth',4);
    set(plot1(2),'LineStyle','-.','Color',[1 0 0],'DisplayName','Arithmetic');
    set(plot1(3),'LineStyle',':','Color',[0 0.498039215803146 0],'DisplayName','Mod. Harmonic');
    legend('Exact', 'Arithmetic', 'Mod. Harm', 'Location', 'Southeast')
    set(gca, 'FontSize',16,'FontWeight','bold');
    xlabel('t')
    ylabel('p') 
    title('Solution profile as a function of t')
    print(3,'-depsc2','test');
    fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/mimetic/';
    fixPSlinestyle('test.eps',[fileLoc 'mh_t502.eps']);
end
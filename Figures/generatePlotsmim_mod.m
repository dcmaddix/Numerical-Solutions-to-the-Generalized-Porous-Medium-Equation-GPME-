%generate plots
figure(1);
N = 25;
fileLoc = 'data/p^m/mimeticp3/';
pa3 = load([fileLoc 'p25_arith3.mat']);
pa5 = load([fileLoc 'p25_arith5.mat']);
pa7 = load([fileLoc 'p25_arith7.mat']);%* in mimetic
pmh3 = load([fileLoc 'p25_mh3.mat']);
pmh5 = load([fileLoc 'p25_mh5.mat']);
pmh7 = load([fileLoc 'p25_mh7.mat']);
[x_cent p_exact] = gentrueSolMim(N);
if (N == 50)
    createfigure_mimpos1(pa3.x_cent, [p_exact(1,:)' pa3.p pmh3.p ...
            p_exact(2,:)' pa5.p pmh5.p p_exact(3,:)' pa7.p pmh7.p])
else
    createfigure_mimpos1(pa3.x, [p_exact(1,:)' pa3.p pmh3.p ...
            p_exact(2,:)' pa5.p pmh5.p p_exact(3,:)' pa7.p pmh7.p])
end
legend('Exact', 'Arithmetic', 'Mod. Harm')
print(3,'-depsc2','test');
fileLoc = '~/Documents/Stanford/Research/ResearchPapers/PorousMediumEquation/Results/p_m/m=3/mimetic/';
fixPSlinestyle('test.eps',[fileLoc 'mh_t25.eps']);
%figure(1)
%createfiguretimeMIM(pa7.t, [p_exact(3,:)' pa7.p pmh7.p])
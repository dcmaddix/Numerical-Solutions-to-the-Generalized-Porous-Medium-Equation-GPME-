close all
N = 25;
kpow = 1;
i = 0;
if (i == 1)
    deg = 30;
elseif(i == 2)
    deg = 45;
elseif (i == 3)
    deg = 75;
end
dir = sprintf('data/2D/m%i/',kpow);
if (kpow ~= 1)
    pa = load(strcat(dir, sprintf('arith%ideg%i.mat', N,deg)));
    ph = load(strcat(dir, sprintf('harm%ideg%i.mat', N,deg)));
    pmhm = load(strcat(dir, sprintf('mhm%i.mat', N)));
else
    pa = load(strcat(dir, sprintf('arith%ideg.mat', N)));
    ph = load(strcat(dir, sprintf('harm%ideg.mat', N)));
    pmhm = load(strcat(dir, sprintf('mhm%ideg.mat', N)));
end
%ref = load(strcat(dir, sprintf('arith200deg%i.mat',deg)));
ref = load(strcat(dir, 'arith200deg.mat'));
if (kpow ~= 1)
    createfiguretime(ref.t(1:end-1), ref.p_t(i,1:end), pa.t(1:end-1), [pa.p_t' ph.p_t' pmhm.p_t(i,1:end)'])
    %createfiguretime(ref.t(1:end-1), ref.p_t(i,1:end), pa.t(1:end-1), [pa.p_t' ph.p_t' pmhm.p_t'])
else
    createfiguretime(ref.t(1:end-1), ref.p_t(i,1:end), pa.t(1:end-1), [pa.p_t(i,1:end)' ph.p_t(i,1:end)' pmhm.p_t(i,1:end)'])
end
legend('Ref. Sol.', 'Arithmetic', 'Harmonic', 'Mod. Harmonic', 'Location', 'southeast')
print(1,'-depsc2','test');
fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/p_t%ideg%i.eps', kpow,N,deg));
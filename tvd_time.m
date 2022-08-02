%compute total variation time-constant at 1 in space
%sum(abs(p(2:end)-p(1:end-1))
%compute at x probe point sum(p^n+1_i - p^n_i))
N = 50;
pa = load('p_t_arith.mat');
ph = load('p_t_harm.mat');
pi = load('p_t_int.mat');
psam = load('p_t_sam.mat');
pe = load('p_texact');
tv_arith = sum(abs(pa.p_t(:,2:end-1)-pa.p_t(:,1:end-2)),2);
tv_harm = sum(abs(ph.p_t(:,2:end-1)-ph.p_t(:,1:end-2)),2);
tv_int = sum(abs(pi.p_t(:,2:end-1)-pi.p_t(:,1:end-2)),2);
tv_sam = sum(abs(psam.p_t(:,2:end-1)-psam.p_t(:,1:end-2)),2);
tv_exact = sum(abs(pe.p_t_exact(:,2:end-1)-pe.p_t_exact(:,1:end-2)),2);
plot(x_cent, tv_sam,x_cent, tv_int,x_cent, tv_arith,x_cent, tv_exact,x_cent, tv_harm)
legend('SAM', 'Int', 'Arith','exact', 'harm')
set(gca, 'FontSize',16,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
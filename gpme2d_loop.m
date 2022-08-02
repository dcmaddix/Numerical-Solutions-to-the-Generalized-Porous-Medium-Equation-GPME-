addpath('DerivativeDiscretizations/')
close all
unitcircle = true;
planar = false;
oneD = false;
avg = 'arith';
%avg = 'harm';
%avg = 'mhm';
crossterms2d = true;
coeff = 1;
nolabel = false;
mhm = false;
if (strcmp(avg, 'mhm'))
    mhm = true;
end
matrix = true;
ninepoint = false;
kpow = 3;
timeScheme = 'ForwardEuler';
%timeScheme = 'BackwardEuler';
%timeScheme = 'RK2_TVD';
eps = 1e-9; %parameter for minimum value
c = 1; %0.4 in mimetic paper parameter
degenerate = true;
sam = false;
saveRes = false;
N = 25;
N1 = N;
deg = 30;
dx = 1/ N;  
dy = dx;
%For planar case make unit normal vector for direction of travel
nw = ones(2,1);
nw = nw / norm(nw);
if (unitcircle)
    lb = -1;
    ub = 1;
else
    lb = 0;
    ub = 1; %y upper bound is 3 in mimetic case
end
x = lb:dx:ub; %lb -1 for whole circle
y = x;
dt = dx^2/32; %32
%dt = dx/8;
nt = 0.08; %5
if (planar)
    nt = 0.3; %0.5 0.7
end
%Ste problem parameters
pstar = 0.5;
kmax = 1.0;
kmin = 0.0;
[x, y] = meshgrid(x,y);
%in this example for each k(p(x,y,t)) below the line is k_min and each
%k(p(x,y,t)) above the line is k_max check corresponding pressure value
%linear initial pressure distribution
phi = (x + y - pstar) / sqrt(2);
p = x + y; %linear ic
if (unitcircle)
    p = (x.^2 + y.^2);
    if (degenerate)
        %p((x.^2+y.^2) >= 0.5) = 0.1;
        %p((x.^2 + y.^2) < 0.5) = 1; luiz test case similar to 1d
        p((x.^2 + y.^2) < 1) = 0.1;
    end
    phi = 1 - sqrt(x.^2+y.^2); %distance to unit circle
end
if (planar) %going to try case where traveling wave not alligned with grid
    p = eps^(1/kpow)*ones(size(x,1), size(x,2));
end
%contour x+y = 0.5 contour(x,y,p, 0.5)
figure(1);
h1 = surf(x,y,p,'LineStyle','none');
view([-0.5 90]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca, 'FontSize',20,'FontWeight','bold');
if (~nolabel)
    xlabel('x')
    ylabel('y')
    zlabel('p')
    %title('Solution profile as a function of (x,y)')
    title('Cross-terms')
    colorbar
end
set(h1,'edgecolor', 'none')
colormap(paruly)
if (nolabel)
    set(gca, 'Xtick',[])
    set(gca, 'Ytick',[])
end
%grid off
%title('Absolute error as a function of (x,y)')
if (saveRes)
    if (planar)
%       saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/planar/%s%i/p_initt%i.jpg', avg,N,nt));
    else
        print(1,'-depsc2','test');
        if (nolabel)
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/IC/p_init%ino.eps', N));
        else
            fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/IC/p_init%i.eps', N));
        end
        %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/pme/IC/p_init%i', N));
    end
end
if (kpow == -1)
    figure(2)
    contour(x,y,p,0.5)
    colorbar
end
k = zeros(size(p));
k(p >= pstar) = kmax;
k(p < pstar) = kmin;
figure(3)
if (kpow == -1)
    surf(x,y,k)
    %view([-0.5 90]);
    xlabel('x')
    ylabel('y')
    title('k')
    figure(6)
    contour(x, y, phi)
    title('Initial phi')
    colorbar
end
%initialize distnace from every point to line
t = 0:dt:nt;
%Store pressure value at (x,y) = (.5,.35)
[coordx, coordy] = find(x == 0.76 & y == -0.6); %(.8,.64)
[coordx, coordy] = find(x == 0.9 & y == 0.28);
%30 deg angle unit circle
[coordx30, coordy30] = find(x == 0.84 & y == 0.52); %approx (sqrt(3)/2,0.5)
[coordx45, coordy45] = find(abs(x - 0.68) < 1e-15 & abs(y - 0.68) < 1e-15); %approx (sqrt(2)/2,sqrt(2)/2)
[coordx75, coordy75] = find(x == 0.96 & y == 0.24); %approx ((sqrt(6)-sqrt(2))/4,(sqrt(6)+sqrt(2))/4)
if (planar)
    [coordx, coordy] = find(x == 0.12 & y == 0.08);
end
if (~planar)
    p_t = zeros(3,length(t)-1); %still 1d plot cause p is scalar
else
    p_t = zeros(1,length(t)-1);
end
k_A = @(x,y) (x+y)/2;
k_H = @(x,y) (2*x.*y) ./ (x + y); %oscillations with harm avg
i = 2:size(p,1)-1;
j = 2:size(p,2)-1;
%make FD matrix
if (unitcircle)
    N = 2*N; %double the points on domain [-1,1] of length L = 2
end
I = speye(N-1,N-1);
I2d = speye((N-1)^2,(N-1)^2);
e = ones(N-1,1);
%remove boundaries
left_bound = find(x == lb);
top_bound = find(y == lb);
right_bound = find(x == ub);
bot_bound = find(y == ub);
ind = setdiff(1:length(p(:)), [left_bound;top_bound;right_bound;bot_bound])';
lambda = dt / dx^2;
for n = 1:length(t)-1
    if (planar) %growing in tiem left boundary condition
        if (oneD)
            p(:,1) = (c^2*kpow*t(n))^(1/kpow);
            p(:,end) = eps^(1/kpow);
            %homogenous Neiman on remaining sides;
            p(1,:) = p(2,:);
            p(end,:) = p(end-1,:);
        else %(x,y) dot n
            %x = 0
            p(:,1) = (kpow*c * pos(c*t(n)-nw(2)*y(:,1))).^(1/kpow);
            %x = 1
            p(:,end) = (kpow*c * pos(c*t(n)-(nw(1)+nw(2)*y(:,1)))).^(1/kpow);
            %y = 0
            p(1,:) = (kpow*c * pos(c*t(n)-nw(1)*x(1,:))).^(1/kpow);
            %y = 1
            p(end,:) = (kpow*c * pos(c*t(n)-(nw(1)*x(1,:) + nw(2)))).^(1/kpow);
            p(p == 0) = eps^(1/kpow);
        end
    end
    %pold = p;
    %phiold = phi;
    %temporal plot still in one dimension
    if (planar)
        p_t(n) = p(coordx,coordy);
    else
        p_t(1,n) = p(coordx30,coordy30);
        p_t(2,n) = p(coordx45,coordy45);
        p_t(3,n) = p(coordx75,coordy75);
    end
    %Update interior leaving boundaries fixed
    %find speed near interface
    %[indx,indy] = find(p(2:end,:) >= pstar & p(1:end-1,:) < pstar);
    %indx = indx(1);
    %indy = indy(1);
    %Vx = (pold(indx+1,indy+1) - pold(indx,indy+1)) / (dx * pold(indx,indy+1));
    %Vy = Vx;
    if (~sam)
        if (strcmp(avg,'arith'))
            k_i = @(x,y)k_A(x,y);
        else
            k_i = @(x,y)k_H(x,y);
        end
            if (kpow == -1)
                k(p >= pstar) = kmax;
                k(p < pstar) = kmin;
            else
                k = p.^kpow;
            end
            k_plusi = k_i(k(i+1,j), k(i,j)); %j vector so this is matrix
            k_minusi = k_i(k(i-1,j), k(i,j));
            k_plusj = k_i(k(i,j+1), k(i,j));
            k_minusj = k_i(k(i,j-1), k(i,j));
            if (matrix)
                k_plusjnew = k_plusj(:,1:end-1); %last column in F
                k_plusinew = k_plusi;
                k_plusinew(end,:) = 0; %remove boundary values
                k_minusinew = k_minusi;
                k_minusinew(1,:) = 0;
                k_minusjnew = k_minusj(:,2:end);
                %pad because only tkaes bottom of column if too long
                D_k = spdiags([[zeros(N-1,1);-k_plusjnew(:)], ...
                        -[0; k_plusinew(1:end-1)'],...
                    k_plusi(:)+k_plusj(:)+k_minusi(:)+k_minusj(:), ...
                    -[k_minusinew(2:end)'; 0], -[k_plusjnew(:); zeros(N-1,1)]],...
                    [N-1 1 0 -1 -(N-1)], (N-1)^2, (N-1)^2);
                %From PoissonMatrix2D make 2D finite diff matrix
                if (strcmp(timeScheme, 'ForwardEuler') || ...
                    strcmp(timeScheme, 'RK2_TVD')) 
                    M = I2d - lambda*D_k;
                elseif (strcmp(timeScheme, 'BackwardEuler')) 
                    M = I2d + lambda*D_k;
                elseif (strcmp(timeScheme, 'CrankNicolson')) %Crank Nicolson
                    M  = I2d + (1-theta)*lambda*D_k;
                end
                F = zeros(size(ind));
                F(1:N-1) = k_minusj(1:N-1) .* p(2:N);
                F(end-(N-2):end) = k_plusj(end-(N-2):end) .* p(end-1-(N-2):end-1);
                for m = 1:N-1
                    %row N+1 boundary
                    F((m-1)*(N-1)+N-1) = F((m-1)*(N-1)+N-1) + ...
                                     k_plusi(m*(N-1))*p(m*(N+1)+N+1);
                    %row 1 boundary
                    F((m-1)*(N-1)+1) = F((m-1)*(N-1)+1) + k_minusi((m-1)*(N-1)+1)*p(m*(N+1)+1);
                end
                if (strcmp(timeScheme, 'ForwardEuler'))
                    p(ind) = M*p(ind) + lambda * F; %need to remove boundaries or pad M
                elseif (strcmp(timeScheme, 'BackwardEuler'))
                    p(ind) = M \ (p(ind) + lambda * F); 
                elseif (strcmp(timeScheme, 'RK2_TVD'))
                    p1 = p;
                    p1(ind) = M*p(ind) + lambda * F; %p^(1) = p^n+ dt * D_k (FE step) m = 1 first stage
                    %do second stage
                    p(ind) = 0.5*(p(ind) + p1(ind)) + 0.5*lambda*(F - D_k*p1(ind));
                end
            else
                FR = k_i(k(i+1,j), k(i,j)) .* (p(i+1,j) - p(i,j)) / dx; 
                FL = k_i(k(i-1,j), k(i,j)) .* (p(i-1,j) - p(i,j)) / dx;
                FU = k_i(k(i,j+1), k(i,j)) .* (p(i,j+1) - p(i,j)) / dy;
                FD = k_i(k(i,j-1), k(i,j)) .* (p(i,j-1) - p(i,j)) / dy;
                p(i,j) = p(i,j) + dt / dx * (FR + FL + FD + FU);
            end
            %fill in boundary conditions in forcing vector
            if (mhm || crossterms2d) %MHM: Add terms to fic harmonic
                %try higher order derivatives to see if helps!
                dp_dx = (p(i+1,j) - p(i-1,j)) / (2*dx); %2nd order central
                dp_dy = (p(i,j+1) - p(i,j-1)) / (2*dy);
                %experimenting with new derivative because bad convergence
                %behavior on new solution
                if (kpow > 1 && planar)
                    ghostNodes_t = {p(1,j) p(end,j)};
                    dp_dx = central_4thOrder2D(p, ghostNodes_t, dx,true);
                    %dp2_dx2 =  central_2ndDeriv_4thOrder2D(p, ghostNodes_t, dx, true);
                    ghostNodes_t = {p(i,1) p(i,end)};
                    dp_dy = central_4thOrder2D(p, ghostNodes_t, dy,false);
                end
                %try nine point scheme!
                dp2_dx2 = (p(i+1,j) - 2*p(i,j) + p(i-1,j)) / dx^2;
                dp2_dy2 = (p(i,j+1) - 2*p(i,j) + p(i,j-1)) / dy^2;
                if (ninepoint)
                    %assumign dx = dy
                      %dp2_dx2 = (5/6 * (p(i+1,j) + p(i-1,j)) - 1/6 * (p(i,j+1) + p(i,j-1)) +...
                      %          1/12 * (p(i+1,j+1) + p(i+1,j-1) + p(i-1,j+1) + p(i-1,j-1)) - 5/3 *p(i,j)) / dx^2;
                     % dp2_dy2 = (-1/6 * (p(i+1,j) + p(i-1,j)) + 5/6 * (p(i,j+1) + p(i,j-1)) + ...
                     %          1/12 * (p(i+1,j+1) + p(i+1,j-1) + p(i-1,j+1) + p(i-1,j-1)) - 5/3* p(i,j)) / dy^2;
                     cross_terms = 1/(12*dx^2) * (p(i+1,j+1) - 2*p(i+1,j) + p(i+1,j-1) -...
                          2*(p(i,j+1) - 2*p(i,j) + p(i,j-1)) + p(i-1,j+1) - 2*p(i-1,j) + p(i-1,j-1));
                     dp2_dx2 = dp2_dx2 + cross_terms;
                     dp2_dy2 = dp2_dy2 + cross_terms;
                end
                %dp2_dy2 =  central_2ndDeriv_4thOrder2D(p, ghostNodes_t,
                %dy, false);
                %Compute pme derivatives
                k_p = kpow*p(i,j).^(kpow-1);
                k_pp = kpow*(kpow-1)*p(i,j).^(kpow-2);
                k_ppp = kpow*(kpow-1)*(kpow-2)*p(i,j).^(kpow-3);
                B = 3/4*k_pp;
                deltaB = -3/4*k_p.^2.*k(i,j).^(-1);
                F = 1/6*k_ppp;
                deltaF = -1/4 * (2*k_p.*k_pp.*k(i,j).^(-1) - k_p.^3.*k(i,j).^(-2));
                %Cancel anti-diffusive terms causing temporal oscillations
                %slightly different discretization check
                if (mhm)
                    p(i,j) = p(i,j) - dt * (B + deltaB).*(dx^2.*dp_dx.^2.*dp2_dx2...
                                                      + dy^2.*dp_dy.^2.*dp2_dy2); %y derivatives
                    %cancel advective term
                    p(i,j) = p(i,j) - dt * (F + deltaF).*(dx^2 .* dp_dx.^4 ...
                                                      + dy^2 .* dp_dy.^4);
                end
                if (crossterms2d)
                    cross_term = 3*(k_p.^2 + k(i,j).*k_pp);
                    p(i,j) = p(i,j) - coeff*dt * (-dt/2 * cross_term) ....
                              .*(dp_dy.^2.*dp2_dx2 + dp_dx.^2.*dp2_dy2);
                end
                %cancel cross term (can switch on or off-see effect without
                %other terms)
            end
    else
        for i = 2:size(p,1)-1
            for j = 2:size(p,2)-1
                %FIND POINT CLOSEST TO THE INTERFACE USING THE MINIMUM
                %DISTANCE AND USE RH CONDITION to calculate difference
                Vx = (pold(i+1,j) - pold(i,j)) / (dx * pold(i,j)); %need flux definition when next to interface not zero
                Vy = (pold(i,j+1) - pold(i,j)) / (dy * pold(i,j));
                %fR = -(kmin .* (p(i+1,j) - p(i,j)) / dx);
                %fL = -(kmax .* (p(i,j) - p(i - 1, j)) / dx);
                %Vx = (fR - fL) / (p(i+1,j)-p(i,j));
                %fR = -(kmin .* (p(i,j+1) - p(i,j)) / dx);
                %fL = -(kmax .* (p(i,j) - p(i, j - 1)) / dx);
                %Vy = (fR - fL) / (p(i,j+1)-p(i,j));
                phi_x = (phiold(i+1,j) - phiold(i,j)) / dx;
                phi_y = (phiold(i,j+1) - phiold(i,j)) / dy;
                %BUG in speed of level set equation
                phi(i,j) = phiold(i,j) + dt*[Vx Vy]*[phi_x; phi_y]; %tested correct sign
                dx_star = abs(phi(i,j)); %distance to the interface
                %each gridpoint check all four neighbors
                %1) Right neighbor
                if (pold(i,j) < pstar && pold(i+1,j) >= pstar) %can generalize to have both cases
                    FR = (pstar - pold(i,j)) / dx_star;
                else
                    FR = (pold(i+1,j) - pold(i,j)) / dx; 
                end
                %Left neighbor
                if (pold(i,j) >= pstar && pold(i-1,j) < pstar)
                    FL = (pstar - pold(i,j)) / dx_star;
                else
                    FL = (pold(i-1,j) - pold(i,j)) / dx;
                end
                %Upper neighbor
                if (pold(i,j) < pstar && pold(i,j+1) >= pstar)
                    FU = (pstar - pold(i,j)) / dx_star;
                else
                    FU = (pold(i,j+1) - pold(i,j)) / dy;
                end
                %Down neighbor
                if (pold(i,j) >= pstar && pold(i,j-1) < pstar)
                    FD = (pstar - p(i,j)) / dx_star; %goes above dx 
                else
                    FD = (pold(i,j-1) - pold(i,j)) / dy;
                end
                %assuming monotonicty
                if (pold(i,j) >= pstar)
                    k_SAM = kmax;
                else
                    k_SAM = kmin;
                end
                p(i,j) = pold(i,j) + k_SAM * dt / dx * (FR + FL + FD + FU);
            end
        end
    end
    %Neumann boundary conditions on signed distance function
    %phi(1,:) = phi(2,:);
    %phi(:,1) = phi(:,2);
    %phi(:,end) = phi(:,end-1);
    %phi(end,:) = phi(end-1,:);
end %make a function for each to reduce errors
f = figure(4);
%surf(x,y,p, 'LineStyle','none')
surf(x(2:end-1,2:end-1),y(2:end-1,2:end-1),(-dt/2 * 2*(3*p(i,j).^2 + p.^3.*6) ....
                              .*(dp_dy.^2.*dp2_dx2 + dp_dx.^2.*dp2_dy2))
colormap(paruly)
view([-0.5 90]);
if (planar)
    view([14.5 16]);
end
if (~planar)
    %shading interp
end
set(gca, 'FontSize',20,'FontWeight','bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
if (~nolabel)
    xlabel('x')
    ylabel('y')
    title('Solution profile as a function of (x,y)')
    title(sprintf('Cross terms for %s', avg))
    colorbar
    set(gca, 'FontSize',20,'FontWeight','bold');
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    saveas(f,sprintf('~/Desktop/cross_terms/%s.jpg',avg));
   % print(4,'-depsc2','test');
    %fixPSlinestyle('test.eps',sprintf('~/Desktop/cross_terms/%s.jpg',avg));
end
if (nolabel)
    set(gca, 'Xtick',[])
    set(gca, 'Ytick',[])
end
if (saveRes)
    if (~planar)
        print(4,'-depsc2','test');
        fixPSlinestyle('test.eps',sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%i.eps', kpow,avg,N1));
        %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_final%i.jpg', kpow,avg,N));
    else
        %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/planar/%s%i/p_finalt%i.jpg', avg,N,nt));
    end
end
if (kpow == -1)
    figure(5)
    contour(x, y, phi)
    colorbar
end
f = figure(7);
plot(t(1:end-1),p_t)
ylabel('p_t')
xlabel('t')
title('Pressure at (x,y) vs time')
if (saveRes)
    if (~planar)
        %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/pme/m%i/%s/p_t%i.jpg', kpow,avg,N));
    else
        %saveas(f, sprintf('~/Documents/Stanford/Research/2Dres/planar/%s%i/p_tt%i.jpg', avg,N,nt));
    end
end
if (planar)
   save(sprintf('data/2D/m%i/planar/%s%i.mat',kpow,avg,N),'x','y','t', 'p_t','p')
else
    if (strcmp(timeScheme, 'BackwardEuler')||strcmp(timeScheme, 'RK2_TVD'))
        save(sprintf('data/2D/m%i/%s/%s%i.mat',kpow,timeScheme,avg,N1),'x','y','t', 'p_t','p') %deg
    elseif (crossterms2d)
        save(sprintf('data/2D/m%i/%s%ict.mat',kpow,avg,N1),'x','y','t', 'p_t','p') %deg
    else
        save(sprintf('data/2D/m%i/%s%i.mat',kpow,avg,N1),'x','y','t', 'p_t','p') %deg
    end
end
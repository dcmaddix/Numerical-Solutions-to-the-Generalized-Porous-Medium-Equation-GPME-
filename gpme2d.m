N = 10;
dx = 1/ N;
dy = dx;
x = 0:dx:1;
y = x;
dt = dx^2/32;
nt = 0.05;
%Ste problem parameters
p_star = 0.5;
kmax = 1.0;
kmin = 0.0;
[x, y] = meshgrid(x,y);
%in this example for each k(p(x,y,t)) below the line is k_min and each
%k(p(x,y,t)) above the line is k_max check corresponding pressure value
%linear initial pressure distribution
p = x + y;
%contour x+y = 0.5 contour(x,y,p, 0.5)
figure(1)
surf(x,y,p)
view([-0.5 90]);
xlabel('x')
ylabel('y')
title('p')
colorbar
figure(2)
contour(x,y,p,0.5)
k = zeros(size(p));
k(p >= p_star) = kmax;
k(p < p_star) = kmin;
figure(3)
surf(x,y,k)
%view([-0.5 90]);
xlabel('x')
ylabel('y')
title('k')
%ind_shock = find(p(i - 1) >= p_star & p_star > p(i));
i = 2:N; %check four cases where cross boundary
j = i; %equal sized domain
count = 1;
pairs = [];
phi = zeros(size(x));
for k= 1:length(i)
    for l = 1:length(j)
        pairs(count,1:2) = [i(k),j(l)];
        count = count + 1;
    end
end
for k = 1:size(x,1)
    for l = 1:size(y,1)
        %distance from ax+by +c = 0 to (x0,y0) = |ax0+by0+c|/sqrt(a^2=b^2)
         phi(k,l) = (x(k, l) + y(k, l)- 0.5) / sqrt(2); %signed distnace function
    end
end
figure(6)
contour(x, y, phi)
title('Initial phi')
%initialize distnace from every point to line
t = 0:dt:nt;
for n = 1:1000
    pold = p;
    %need to check four corner cases otherwise can do regular update
    [iR, jR] = find(p(i - 1,j) < p_star & p(i,j) >= p_star); %careful with iR and iR+1
    %how to access (i,j) p(iR, jR+1) and for right neighbor p(iR+1, jR + 1)
    [iU,jU] = find(p(i,j - 1) < p_star & p(i,j) >= p_star);
    %how to access (i,j) p(iU+1, jU) and for right neighbor p(iU+1, jU + 1)
    %special stencils
    %Find all not near interface check all four cases
    ind_shock = [[iR, jR+1];[iR+1, jR+1];[iU+1,jU]; [iU+1,jU+1]]; %check if boundary indices don't change
    reg = setdiff(pairs, ind_shock, 'rows');
    %standard laplacian away from itnerface
    [kmax_i, kmax_j] = find(p >= 0.5);
    X = intersect(reg, [kmax_i, kmax_j], 'rows');
    %level set evolution
    %all same vlaue regardless of index constant along lines for position
    %first
    ind = 2;
    %redefine in terms of fluxes: change back to numerical speed
    V_x = (p(iR(ind),jR(ind)+1) - p(iR(ind)+1,jR(ind)+1))...
        ./ (dx * p(iR(ind),jR(ind)+1))*ones(size(phi(2:end,2:end))); %need flux definition when next to interface not zero
    V_y = (p(iU(ind)+1,jU(ind)) - p(iU(ind)+1,jU(ind)+1)) ./ (dy * p(iU(ind)+1,jU(ind)))*ones(size(phi(2:end,2:end)));
    %ls equation
    phi_x = (phi(2:end,2:end) - phi(1:end-1,2:end)) / dx;
    phi_y = (phi(2:end,2:end) - phi(2:end,1:end-1)) / dy;
    for l = 1:size(phi,1)-1
        for m = 1:size(phi,1)-1
            phi(l,m) = phi(l,m) - dt*[V_x(l,m) V_y(l,m)]*[phi_x(l,m); phi_y(l,m)];
        end
    end
    phi(:,end) = phi(:,end-1);
    phi(end,:) = phi(end-1,:);
    for m = 1:length(X)
        p = updateStencil(pold,kmax, dt, dx, dx, false, ...
                          false, false, false, X(m,1), X(m,2));
    end 
    [kmin_i, kmin_j] = find(p < 0.5);
    X = intersect(reg, [kmin_i, kmin_j], 'rows');
    for m = 1:length(X)
       p = updateStencil(pold,kmin, dt, dx, dx, false, ...
                          false, false, false, X(m,1), X(m,2));
    end %look up distance in phi function: write as function
    if (iR > 1) %check if double overlap indices
        p = updateStencil(pold,kmin, dt, dx, abs(phi(iR, jR+1)), true, ...
                          false, false, false, iR, jR+1);
    end%finish stencil for 
    p = updateStencil(pold,kmax, dt, dx, abs(phi(iR+1,jR+1)), false, ...
                     true, false, false, iR+1,jR+1);
    %need to detect overlap c-ases!!! AND Correct speed function
    if (jU > 1)
        p = updateStencil(pold,kmin, dt, dx, abs(phi(iU+1, jU)), false, ...
                     false, false, true, iU+1,jU);
    end
    %check input arguments for this stnecil~
    p = updateStencil(pold,kmax, dt, dx, abs(phi(iU+1, jU+1)), false, ...
                     false, true, false, iU+1,jU+1);
end %make a function for each to reduce errors
figure(4)
surf(x,y,p)
view([-0.5 90]);
xlabel('x')
ylabel('y')
title('p')
colorbar
figure(5)
contour(x, y, phi)
% p(iR, jR+1) = pold(iR, jR+1) + kmin * dt / dx * ((pold(iR+1, jR+1) - pold(iR, jR+1))...
        %        ./ abs(phi(iR, jR+1)) - (pold(iR, jR+1) - pold(iR-1, jR+1)) / dx ...
         %       - (pold(iR, jR+1) - pold(iR, jR+2)) / dx - (pold(iR, jR+1)
         %       - pold(iR, jR)) / dx);
  %p(iU+1, jU+1) = pold(iU+1, jU+1) + kmax * dt / dx * (-(pold(iU+1, jU+1) - pold(iU+1, jU))...
    %            ./ abs(phi(iU+1, jU+1)) - (pold(iU+1, jU+1) - pold(iU, jU+1)) / dx ...
    %           +(pold(iU+2, jU+1) -pold(iU+1, jU+1)) / dx - (pold(iU+1,
    %           jU+1) - pold(iU+1, jU)) / dx);
     %p(iU+1, jU) = pold(iU+1, jU) + kmin * dt / dx * ((pold(iU+1, jU+1) - pold(iU+1, jU))...
         %       ./ abs(phi(iU+1, jU)) - (pold(iU+1, jU) - pold(iU, jU)) / dx ...
          %      +(pold(iU+2, jU) -pold(iU+1, jU)) / dx - (pold(iU+1, jU) -
          %      pold(iU+1, jU-1)) / dx);
          %p(iR+1, jR+1) = pold(iR+1, jR+1) + kmax * dt / dx * (-(pold(iR+1, jR+1) - pold(iR, jR+1))...
     %           ./ abs(phi(iR+1, jR+1)) - (pold(iR+1, jR+1) - pold(iR, jR+1)) / dx ...
      %          -(pold(iR+1, jR+1) - pold(iR+1, jR)) / dx + (pold(iR+1,
      %          jR+2) - pold(iR, jR+1)) / dx);

clear all
%close all

% Required modules
mrstModule add ad-core % ad-fi
tic

L = 1;
nVc = 200;

pLeft = 1.0;
pRight = 0.0;

dtOverDx2 = 1/4;
dx = L/nVc;
dt = dtOverDx2 * dx^2;
alpha = 1.0; 
eps = 0.05; %(left "half" band ramp)
nPow = 3;

probeLocation = 0.13;

nNewton = 1;

%  |   *   |   *   | ... |   *   |
% f1   c1  f2  c2  f3        c_N  f_{N+1} 
faces = 1:nVc+1;
cells = 1:nVc;
facesR = 2:nVc+1;
facesL = 1:nVc;
cellsR = 2:nVc;
cellsL = 1:nVc-1;
xC = [dx/2; dx/2 + dx*cellsL'];

% Initial condition
p0 = zeros(numel(cells),1);

p = initVariablesADI(p0);
pL = ADI(pLeft,zeros(1,nVc));
pR = ADI(pRight,zeros(1,nVc));

% k definition
kMin = 0.01;
kMax = 1;
% k = @(P) (kMin + (kMax-kMin)*(P>=0.5));
% k = @(P) ( max(min(0.5*(kMin+kMax) + 10*(P-0.5), kMax),kMin) );
beta = (kMax-kMin)/((2*eps)^nPow);
% k = @(P) ( max(min(0.5*(kMin+kMax) + beta*(P.^2-0.25), kMax),kMin) );
% k = @(P) ( max(min((kMin + beta*((P-0.5+eps).^nPow)), kMax),kMin) );
beta = (kMax-kMin)/((eps)^nPow);
k = @(P) ( max(min((kMin + beta*((P-0.5+eps).^nPow)), kMax),kMin) );

figure;
pp = 0:0.01:1;
plot(pp, k(pp));

pOld = p.val;
pfOld = [pLeft; 0.5*( p.val(cellsL)+p.val(cellsR) ); pRight];
kfOld = k(pfOld);

pf = [pL; 0.5*( p(cellsL)+p(cellsR) ); pR];   
px = [p; p(end)+2*(pR-p(end))];
px = px - pf;
px = 2*px; 

time = 0:dt:0.05;
nTime = numel(time);
probe = zeros(nTime,1);
[val,iProbe] = (min(abs(xC-probeLocation)));
figure;
plot(xC,p0)
% time Loop:
for iTime = 1:nTime
        
    for iNewton = 1:nNewton    
        pf.val = [pLeft; 0.5*( p.val(cellsL)+p.val(cellsR) ); pRight];   
        px.val = [p.val; p.val(end)+2*(pRight-p.val(end))];
        px.val = px.val - pf.val;
        px.val = 2*px.val;
%         kf = ((iNewton-1)*kfOld + (nNewton-iNewton+1)*k(pf.val))/nNewton;
        kf = ((alpha-1)*kfOld + (alpha)*k(pf.val));
        kfOld = kf;        
        eqP = (p-pOld) ...
                   - (dtOverDx2)*(...
                                    kf(facesR).*px(facesR)...
                                  - kf(facesL).*px(facesL)...
                                 ) ; 
        dp = -eqP.jac{1}\eqP.val;
        p.val = p.val + dp;
    end


    pOld = p.val;
    kfOld = kf;
    

    probe(iTime,1) = p.val(iProbe);

    % plot(xC, p.val);
    % hold on
    % pause(0.5);
end
toc

figure(2);
plot(xC, p.val);
hold on

figure(3);
plot(time, probe); 
hold on
%why is h always input as 1?
function [D, F, H, HI] = set_up_1D_matrix_SBPSAT(n,h,c,bc,g_type,g,type);
% output: D - full 2nd derivative discretization with boundary terms (SAT terms)
% F - right hand side
% H - norm operator
% HI - inverse of norm operator 

% c - vector with coefficients (including bdy values) 
% n - number of grid points (including bdy)
% h - grid spacing
% bc - bdy condition implementation type: 'weak' or 'strong'
% g - boundary conditions (values)
% g_type - boundary conditions type: Dirichlet (D) of Neumann (N)

% aa, bb are alpha, beta in new penalty system
% type - averaging types: Harmonic (H), arithmetic (A) 

% construct norm operator: HD1 = Q used as quad rule for discrete
% integration
H=spdiags(ones(n,1),0, n, n); %same as withou t0 specifiying main diagonal
%where did coefficients in upper block come from tau penalty param?
H(1:2,1:2)=[0.1e1 / 0.4e1 0; 0 0.3e1 / 0.4e1;]; %want smallest size cell on first and last due to bdr discret
H(n-1:n,n-1:n)=fliplr(flipud(H(1:2,1:2))); %flipud: flid matrix up/down dir columns preserved and rows flipped in up/down dir
H=H*h;%fliplr flips left and right [0.25 0; 0 .75] goes to flipud [0 0.75; 0.25 0] flip lr [0.75 0; 0 0.25]
HI=inv(H); %inverting diagonal matrix just reciprol entries so will be fast

% construct first derivative on boundary
S_U=[-2, 2]/h;
S_1=zeros(1,n);
S_1(1:2)=S_U;
S_m=zeros(1,n);
S_m(n-1:n)=fliplr(-S_U);%flip left to right so 0's first then -2 and 2
%only need to store first row and last row of S since identity on interior
%forward difference on bdry


% set up matrix D1 approximating first derivative
Q = 0.5*(spdiags(ones(n-1,1),1, n-1, n-1)-spdiags(ones(n-1,1),-1, n-1, n-1)); %tridiagonal 0 on main diag -1/2 on lower diag
Q(1,1) = -0.5;
Q(n,n) = 0.5;
%1/2 on upper diag. Q+Q^T = E = diag(-1,..,1) (SBP prop)
%want -1/2 in(1,1) and 1/2 in (n,n) diagonal so must be skew symmetric

%compute first derivative operator D1 = H_inv*Q
D1=HI*Q;

% start setting up the discretization of variable coefficient 2nd
% derivative D2 from Ken's paper
% Matrix without boundary conditions will be on form D2 = HI*(-M + BS);
M=sparse(n,n);

% RHS
F = sparse(n,1); %sparse RHS vector
g_vec = sparse(n,1); %bdry vector with Dirichlet bc
g_vec(1) = g(1);
g_vec(end) = g(2);


btilde = zeros(size(c)); %matrix consisting of harm or arith averages of perms
nn = length(btilde);

cbar = zeros(nn,1);
cbar(1) = 2*c(1)*c(2)./(c(1) + c(2)); %harmonic averaging on bdries
cbar(end) = 2*c(end-1)*c(end)./(c(end-1) + c(end));

cbar(2:nn-1) = 16*c(1:nn-2).*c(2:nn-1).*c(3:nn)./((c(1:nn-2) + 2*c(2:nn-1) + c(3:nn)).^2);

if strcmp(type,'H') %do harmonic averaging of the permeabilties just as TPFA in FV on the interior
    
    btilde(1:nn-1) = 0.2e1*c(1:nn-1).*c(2:nn) ./(c(1:nn-1) + c(2:nn));
elseif strcmp(type,'A')
    btilde(1:nn-1) = (c(1:nn-1) + c(2:nn)) / 0.2e1;
        
%elseif strcmp(type,'new') cbar to be used in new test
        
%    btilde(j) = (cbar(j) + cbar(j+1))/0.2e1;
else
    fprintf('ERR: Invalid Averaging Type!\n')
    return;
end

%B0_tilde = diag(B0tilde0,0,...,0 (B0)tildeN) choose coeffs depending on
%weak imposition of bc
if strcmp(bc,'weak')
% we keep c(2) and c(end-1) for the bdy
elseif strcmp(bc,'strong')
   % this is if strong imposition:
    c(2) = btilde(1);
    c(end-1) = btilde(end-1);
end


% add boundary conditions: define standard basis vectors
e0 = zeros(n,1);
e0(1) = 1;
E0 = sparse(n,n); %need sparse storage!
E0(1,1) = 1;

eN = zeros(n,1);
eN(n) = 1;
EN = sparse(n,n);
EN(n,n) = 1;

% set second penalty parameter in symmetric discretization, cc_0/N < 1
% more flexibility (less strong imposition) in boundary conditions the closer cc_0/N are to 1
% if cc_0/N get closer to 0 it becomes closer to enforcing the boundary conditions
% exactly
cc_0 = 0.5;
cc_N = 0.5;
%penalty param tau has 1/c coeff

if strcmp(bc,'weak')
    
% alpha is for weak implementation that leads to a monotone scheme
% c0 < btilde(1) is a condition for monotonicity and extending c(2)->c(1)
% directly (cell permeability value bis set to bdy value) also leads to 
% btilde(1) = c(2). So monotonicity condition is c0 < c(2).
% Should preferrably be handled as a condition in upscaling instead of set here.
alpha = 0.9;

c0 = alpha*c(2);
cN = alpha*c(end-1);

% these conditions leads to monotonicity for the second and second last row
% but really bad behavior - like putting in barriers just by the injection
% boundary
%c0 = min(c(2),c(3));
%cN = min(c(end-1),c(end-2));

% set penalty parameters
tau_0 = -1;
tau_N = -1;
%from Annas' paper notation btilde(1) = b0 and c0 = (Btilde0)0
sigma_0 = (2*(btilde(1) - c0).^2)/(btilde(1)) - 2*btilde(1) + 4*c0;
sigma_N = (2*(btilde(end-1) - cN).^2)/(btilde(end-1)) - 2*btilde(end-1) + 4*cN;

% sigma needs to be scaled with h and also there is a condition sigma 
sigma_0 = sigma_0/h/cc_0; %(1/c) param in papers
sigma_N = sigma_N/h/cc_N;

elseif strcmp(bc,'strong')
c0 = btilde(1);
cN = btilde(end-1);
    
% set penalty parameters
tau_0 = -1;
tau_N = -1;

sigma_0 = (2*(btilde(1) - c0).^2)/(btilde(1)) - 2*btilde(1) + 4*c0;
sigma_N = (2*(btilde(end-1) - cN).^2)/(btilde(end-1)) - 2*btilde(end-1) + 4*cN;

% sigma needs to be scaled with h
sigma_0 = sigma_0/h/cc_0;
sigma_N = sigma_N/h/cc_N;

end

for i=2:n-1 %applying 3 columns per row
    M(i,i-1:i+1)=[-btilde(i-1) btilde(i-1)+btilde(i) -btilde(i);]; %averages same as D_k on interior
end

%need to specify separate bdry stencils
M(1:2,1:2)=[2*btilde(1) -2*btilde(1); -2*btilde(1) 2*btilde(1) + btilde(2);];

M(n-1:n,n-1:n)=[btilde(n-2) + 2*btilde(n-1) -2*btilde(n-1); -2*btilde(n-1) 2*btilde(n-1);];
M=M/h;


%bdry term and deriv
BS = sparse(n,n);
BS(1,:) = -c0*S_1;
BS(end,:) = cN*S_m;
%BS=-c0*e0*S_1+cN*eN*S_m; %outer product cannot form explicitly


if strcmp(g_type,'D')
% Dirichlet BC

SAT_1 = tau_0*HI*BS'*E0 + tau_N*HI*BS'*EN;
SAT_2 = sigma_0*HI*E0 + sigma_N*HI*EN;

F = (SAT_1 + SAT_2) * g_vec;

elseif strcmp(g_type,'N')
% Neumann BC, use tau=1 for boundary terms to cancel out with operator

%SAT_1 = -tau_0*HI*BS*E0 - tau_N*HI*BS*EN;
SAT_1 = -tau_0*HI*BS;
SAT_2 = sparse(n,n);


end


% full matrix without bc's
D2 = HI*(-M + BS);
%D2=HI*(-M-c(2)*e_1*S_1+c(n-1)*e_m*S_m);

% with bc's
D = -D2 + SAT_1 + SAT_2;


end



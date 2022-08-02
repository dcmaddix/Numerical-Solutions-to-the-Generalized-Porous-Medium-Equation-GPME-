%function to test scheme I from the 2016 mimetic paper
%given cell-centered values k solving spd system for cell centers p,
%where u = -kgrad(p) can compute flux from gradient of p
function [p_harm p_arith p_FV p_SAT u_harm u_arith u_FV] = testSchemeI(N,k, bdry, F)
   %define harmonic average of cell-centered constants
    k_H = @(k1,k2) (2*k1.*k2) ./ (k1+k2);
    %define arithmetic average of cell-centered constants for comparison with
    %2016 mimetic paper
    k_A = @(k1,k2) (k1 + k2) / 2; %arithmetic averaging causes smoothing
    k_FV = @(k1,k2) (k_H(k1,k2).*k_A(k1,k2)).^(0.5);%input to get simple FV scheme
    [p_harm u_harm] = mimeticSteadyTransportCoeff(k_H, N, k, bdry, F);
    [p_arith u_arith] = mimeticSteadyTransportCoeff(k_A, N,k, bdry, F);
    [p_FV u_FV] = mimeticSteadyTransportCoeff(k_FV, N,k, bdry, F);
    
    %call Anna's SBP-SAT solver
    g_type = 'D'; %Dirichlet conditions
    h = 1;
    % permability field - change to realistic data (cell values

    % Flux  type - 'H' leads to TPFA scheme in interior (and better suited than 
    % 'A' for heterogeneous coefficients but 'A' is common for SBP-SAT
    % methods)
    type = 'H';

    % Imposition of boundary conditions - 'weak' is SAT conditions, 'strong' is
    % exact imposition (via SAT terms)
    bc = 'weak';

    g = bdry; %values of Dir BC
    % compute the variable coefficient 2nd derivative operator D with BCs and RHS F
    [D, F, ~,~] = set_up_1D_matrix_SBPSAT(N+2,h,k,bc,g_type,g,type);
    p_SAT = D\F;
end
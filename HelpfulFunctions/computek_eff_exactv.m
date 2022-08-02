%This function computes k_eff by finding the location of dx* ta each
%timestep. Pass in x_cent to determine which cell shock in
function [k_face k_plot ratio] = computek_eff_exactv(p, i, p_star, k_lower, k_upper, v, dx)
    %smaller index first - Case 1: (i-1,i) and 
    %Case 2: (i,i+1) need to add to ind_shock
    %flags if not set
    k_plot = -1;
    ratio = -1;
    k_face = zeros(length(i),1);
    ind_shock = find(p(i) > p_star & p_star >= p(i + 1)); %determine which cell in
    ind_shock_p = ind_shock;
    %Shock location
    %ind_shock indexed on interior need to add 1 to the index for (i-1,i)
    %and 2 in (i,i+1) case
    if (min(i) > 1) %(i,i+1) case add 1 more to ind_shock called on final cell
        ind_shock_p = ind_shock_p + 1;
    end %change and only interpolate if find shock location do don't have to for RHS
    if (~isempty(ind_shock))
%       %y = store dx*/dx assume in cell x_i in cell i+1
%       Flux from x_star to x_i+1 assuming that the star is in cell i
%         F12 = @(y) (y / k_upper) * (8*dx * (p_star - p(ind_shock_p + 1))*k_lower ...
%                 + v*(p(ind_shock_p) - p(ind_shock_p + 1)) * (1-y)) ...
%                 / (8*dx*(p(ind_shock_p) - p(ind_shock_p + 1))*(1-y)) ...
%                 - (p_star - p(ind_shock_p + 1)) / (p(ind_shock_p) - p(ind_shock_p + 1));
        %Flux from x_i to x_star assuming that stra is in cell i+1
%         Fi1 = @(y) ((1-y) / k_lower) * (k_upper / y - dx*v/8*...
%                     (p(ind_shock_p) - p(ind_shock_p + 1)) / (p(ind_shock_p) - p_star)) ...
%                     - (p_star - p(ind_shock_p + 1)) / (p(ind_shock_p) - p_star);
        %change for quaf formula
        a = -dx * v / 8;
        b = k_lower + dx * v / 8 - ((p(ind_shock_p) - p_star)/ (p(ind_shock_p) - p(ind_shock_p + 1))...
                    * (k_lower - k_upper) / k_upper);
        c = -(p(ind_shock_p) - p_star) / (p(ind_shock_p) - p(ind_shock_p + 1));
        discrim = b^2 - 4*a*c;
        ratio = (-b + sqrt(discrim)) / (2*a);
        if (0 > ratio || ratio > 1)
            ratio = (-b - sqrt(discrim)) / (2*a);
        end
%         Fi1 = @(y) y*((8*k_lower + dx*(1-y)*v) / 8 - (p(ind_shock_p) - p_star) ....
%                     / (p(ind_shock_p) - p(ind_shock_p + 1)) * (k_lower - k_upper) / k_upper) ...
%                     -(p(ind_shock_p) - p_star) / (p(ind_shock_p) - p(ind_shock_p + 1));
% %         F12 = @(y) (1-y) * (8*k_upper - dx*y*v) / 8 - ((p_star - p(ind_shock_p + 1)) ...
% %                 / (p(ind_shock_p) - p(ind_shock_p + 1)) * y * (k_lower - k_upper) / k_lower) ...
% %                 -(p_star - p(ind_shock_p + 1)) / (p(ind_shock_p) - p(ind_shock_p + 1)) * (k_upper / k_lower);
%         %[ratio_12 fval exitflag] = fzero(F12, 0.25);%y in range 0 to 0.5
%         [ratio fval exitflag] = fzero(Fi1, 0.5);%y in range 0.5 to 1
        %ratio  = 1 / ratio;
        %when ratio = 0.5 get harm avg and when ratio = 0 get k_lower ratio = 1
        %get k_upper
        %in cell i
        %if (0 <= ratio_i1 && ratio_i1 <= 1.0)
         psi_1 = ratio * dx;
         psi_2 = dx - psi_1;
         harm_kpsi = 1/((psi_1 / k_upper) + (psi_2/k_lower));
         %compare k_Eff values instead
         k_face(ind_shock) = dx * harm_kpsi * (1 + (psi_2 / k_lower) * (v / 8));
         %k_face_iplus(ind_shock) = dx * harm_kpsi * (1 - (psi_1 / k_upper) * (v / 8));
         k_plot = k_face(ind_shock);
    end
    %Case b: p* >= p_i >= p_{i+1}
    k_face(p_star >= p(i)) = k_lower;
    %Case c: p_i >= p_{i+1} >= p*
    k_face(p(i + 1) >= p_star) = k_upper;
end
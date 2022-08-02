%This function computes k_eff by finding the location of dx* ta each
%timestep
function [k_face k_plot ratio] = computek_eff(p, i, p_star, k_lower, k_upper)
    %smaller index first - Case 1: (i-1,i) and 
    %Case 2: (i,i+1) need to add to ind_shock
    k_plot = -1; %initialize flags
    ratio = -1;
    k_face = zeros(length(i),1);
    ind_shock = find(p(i) > p_star & p_star >= p(i + 1));
    ind_shock_p = ind_shock;
    %Shock location
    %ind_shock indexed on interior need to add 1 to the index for (i-1,i)
    %and 2 in (i,i+1) case
    if (min(i) > 1) %(i,i+1) case add 1 more to ind_shock called on final cell
        ind_shock_p = ind_shock_p + 1;
    end %change and only interpolate if find shock location do don't have to for RHS
    if (~isempty(ind_shock))
        %Use flux equlaity condition at p* to solve for dx/dx*
        ratio = 1 + (k_lower / k_upper) * (p(ind_shock_p + 1) - p_star) ...
                                      / (p_star - p(ind_shock_p));
        %store dx*/dx
        %when ratio = 0.5 get harm avg and when ratio = 0 get k_lower ratio = 1
        %get k_upper
        %0 for k_lower = 0 since ratio = 1 
        %can't use this form for k_lower = 0 ratio = 1 and 0/0
        %ratio  = 1 / ratio;
        %k_face(ind_shock) = (k_upper * k_lower) / ...
                                (k_upper * (1 - ratio) + k_lower * ratio);
        %same as integral average
        k_face(ind_shock) = k_upper * (p_star - p(ind_shock_p)) * ratio ...
                                    / (p(ind_shock_p+1) - p(ind_shock_p));
        k_plot = k_face(ind_shock);
    end
    %Case b: p* >= p_i >= p_{i+1}
    k_face(p_star >= p(i)) = k_lower;
    %Case c: p_i >= p_{i+1} >= p*
    k_face(p(i+1) >= p_star) = k_upper;
end
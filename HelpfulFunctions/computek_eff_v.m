%This function computes k_eff by finding the location of dx* ta each
%timestep
function [k_face k_plot v] = computek_eff_v(p, i, p_star, k_lower, k_upper, vold, dx)
    %smaller index first - Case 1: (i-1,i) and 
    %Case 2: (i,i+1) need to add to ind_shock
    %identify shock location
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
        fR = k_lower .* (p(ind_shock_p + 2) - p(ind_shock_p + 1)) / dx;
        %fR = 0;
        if (ind_shock_p ~= 1)
            fL = k_upper .* (p(ind_shock_p) - p(ind_shock_p - 1)) / dx;
        else
            fL = 0.0;
        end
        v = (fL - fR) / (p(ind_shock_p + 1) - p(ind_shock_p));
        v = min(v, vold); %vold fom last timestep;
        k_face(ind_shock) = v * dx;
        k_plot = k_face(ind_shock);
    end
    %Case b: p* >= p_i >= p_{i+1}
    k_face(p_star >= p(i)) = k_lower;
    %Case c: p_i >= p_{i+1} >= p*
    k_face(p(i+1) >= p_star) = k_upper;
end
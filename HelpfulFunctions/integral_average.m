%This function computes the integral average of k_lower and k_upper at the
%flux face assuming p is monotonically non-increasing
function k_face = integral_average(p, i, p_star, k_lower, k_upper)
    %smaller index first - Case 1: (i-1,i) and 
    %Case 2: (i,i+1) need to add to ind_shock
    k_face = zeros(length(i),1);
    j = i + 1;
    %Shock location
    ind_shock = find(p(i) >= p_star & p_star >= p(j));
    %ind_shock indexed on interior need to add 1 to the index for (i-1,i)
    %and 2 in (i,i+1) case
    ind_shock_p = ind_shock;
    if (min(i) > 1) %(i,i+1) case add 1 more to ind_shock called on final cell
        ind_shock_p = ind_shock_p + 1;
    end
    %Case a: p_i >= p* >= p_{i+1}
    k_face(ind_shock) = 1 ./ ( p(ind_shock_p + 1) - p(ind_shock_p) ) .* ...
                         (k_upper * (p_star - p(ind_shock_p) ) + ...
                          k_lower * (p(ind_shock_p + 1) - p_star) );
    %Case b: p* >= p_i >= p_{i+1}
    k_face(p_star >= p(i)) = k_lower;
    %Case c: p_i >= p_{i+1} >= p*
    k_face(p(j) >= p_star) = k_upper;
end
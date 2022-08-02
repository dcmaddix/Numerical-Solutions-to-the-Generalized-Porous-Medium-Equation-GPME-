%This function computes the conservative form for the discontinuous k case
%with no averaging p_t - A_xx = 0, where A(p) = int_0^p k(p~)dp~
function p_new = compute_newP_noAvg(i, p, F, dt, alpha, lambda, k_lower, k_upper, p_star)
    p_new = zeros(length(i), 1); %defined on itnerior only
    %CASE 1: p_{i-1} >= p* >= p_i >= p_{i+1}
    j = find(p(i - 1) >= p_star & p_star >= p(i));
    %To access p need to add 1 to the index since p(i) defined on itnerior
    %returns indices in range 1:end-2
    p_new(j) = p(j + 1) + (1 / alpha) * (dt * F(j) + lambda * ...
                              (k_lower *(p(j + 2) - p(j + 1)) ...
                             + k_lower * (p_star - p(j + 1)) ...
                             + k_upper * (p(j) - p_star)));
   %CASE 2: p_{i-1} >= p_i >= p* >= p_{i+1}
   %change dpart of p!
   j = find(p(i) >= p_star & p_star >= p(i + 1));
   p_new(j) = p(j + 1) + (1 / alpha) * (dt * F(j) + lambda * ...
                             (k_upper *(p_star - p(j + 1)) ...
                            + k_lower * (p(j + 2) - p_star) ...
                            - k_upper * ( p(j+1) - p(j) )));
   %CASE 3: p_{i-1} >= p_i >= p_i+1 >= p*      
   %All k_upper since larger than p* and use standard laplacian stencil
   j = find(p(i + 1) >= p_star);
   p_new(j) = p(j + 1) + (1 / alpha) * (dt * F(j) + lambda * ...
                              k_upper *(p(j) - 2*p(j + 1) + p(j + 2)));
   %CASE 4: p* >= p_{i-1} >= p_i >= p_{i+1}
   %All k_lower since less than p* and use standard laplacian stencil
   j = find(p_star >= p(i - 1));
   p_new(j) = p(j + 1) + (1 / alpha) * (dt * F(j) + lambda * ...
                              k_lower *(p(j) - 2*p(j + 1) + p(j + 2)));
end
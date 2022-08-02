function F = calcFluxes(p1,p2, dx, isShock, p_star, k_lower, k_upper)
    F = 0;
   if (isShock)
       F = k_upper * (p1 - p_star) / dx + k_lower * (p_star - p2) / dx;
   elseif (p2 >= p_star)
       F = k_upper * (p1 - p2) / dx;
   elseif (p1 <= p_star)
       F = k_lower * (p1 - p2) / dx;
   end
end
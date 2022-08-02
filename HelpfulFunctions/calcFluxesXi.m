function F = calcFluxesXi(p1,p2, dx, isShock, p_star, k_lower, k_upper, xi, isLeft) %may be dx - xi
   %very sensitive to this tolerance@ for 2e-6 too high iwth 4 points 2e-7
   %low 5e-7 best for 4 points
   %test 10 points now!
   tol = 2e-6; %2e-6 for exact solution and dx_shock = dx / 4
   skip_left = xi <= tol;
   skip_right = (dx - xi) <= tol;
   grad_p = (p1 - p2) / dx;
   if (isShock)
       if (isLeft)
           F = k_upper * (p1 - p_star) / xi; %dx will be xi for the shock case
           %skipping small segment with pi > p* use k_lower for cell
            if(skip_left) %doesnt enter here in nonexact case
                F = k_lower * grad_p;
            end
       else %should never enter here!!! problem! ALERT!
           F = k_lower * (p_star - p2) / (dx - xi);
           if (skip_right)
               F = k_upper * grad_p
           end
       end
   %p1 > p2 > p_star left of shock
   elseif (p2 >= p_star)
       F = k_upper * grad_p;
   % p_star > p1 > p2 to left of shock
   elseif (p1 <= p_star)
       F = k_lower * grad_p;
   end
end

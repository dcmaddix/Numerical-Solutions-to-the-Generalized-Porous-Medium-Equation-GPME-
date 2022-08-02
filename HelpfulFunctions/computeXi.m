function [xi_new check] = computeXi(xi, dt, dx_min, dx_shock, dx_plus, p_min, p, p_plus, ...
                            p_plus2, k_upper, k_lower, p_star) %x coor of the pressure
    %Need communications between functions as to whether or not it will be
    %skipped! and then reset when crosses next cell use ind_ref index to
    %change outside
    check = false;  
    tol = dx_shock / 10;
    F_i_minus = calcFluxes(p_min, p, dx_min, ...
                          false, p_star, k_lower, k_upper);
    %if (dx_shock - xi < tol)
        F_i_plus = k_upper * (p - p_plus) / dx_shock;
    %else
        F_i_plus = calcFluxesXi(p, p_plus, dx_shock, ...
                              true, p_star, k_lower, k_upper, xi, true);
    %end
    F_iplus1_minus = calcFluxes(p, p_plus, dx_shock, ...
                                true, p_star, k_lower, k_upper);
    F_iplus1_plus = calcFluxes(p_plus, p_plus2, dx_plus, ...
                                false, p_star, k_lower, k_upper);
    h1 = (dx_min + dx_shock) / 2;
    dpi_dt = (F_i_minus - F_i_plus) / h1;
    h2 = (dx_shock + dx_plus) / 2;
    dpiplus_dt = (F_iplus1_minus - F_iplus1_plus) / h2;
    %add tolerance here:dx_min wrong
    xi_new = xi + dt * (8*(F_i_plus - F_iplus1_minus) - xi * dpi_dt - ...
                        (dx_shock - xi) * dpiplus_dt) / (p - p_plus);
                   
    if (xi_new < xi)
        %check = true
        %use linear extrapolation for the value
        %xi_new = (p_star - p) / (p - p_min) * dx_min
    end
end
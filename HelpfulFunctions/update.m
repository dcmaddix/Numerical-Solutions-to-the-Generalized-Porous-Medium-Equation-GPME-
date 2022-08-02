function p = update(pold, dt, area, F_min, F_plus)
    p = pold + dt / area * (F_min - F_plus);
end


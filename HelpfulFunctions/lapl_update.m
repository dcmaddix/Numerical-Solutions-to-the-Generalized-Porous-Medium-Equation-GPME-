function p = lapl_update(pold, dt, dx, i, k)
    p = pold(i) + k * dt / dx^2 * (pold(i-1) - 2* pold(i) + pold(i+1));
end


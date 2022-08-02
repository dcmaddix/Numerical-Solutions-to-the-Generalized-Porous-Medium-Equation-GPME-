function p = updateStencil(p,k, dt, dx, dx_star, shift_right, shift_left, shift_down, shift_up, i,j)
%false on shaft_right is shift_left
FR = p(i+1, j) - p(i, j);
FL = p(i-1, j) - p(i, j);
FD = p(i, j-1) - p(i, j);
FU = p(i, j+1) - p(i, j);
if (shift_right)
    FR = FR ./ dx_star;
else
    FR = FR / dx;
end
if (shift_left)
    FL = FL ./ dx_star;
else
    FL = FL / dx;
end
if (shift_down)
    FD = FD ./ dx_star;
else
    FD = FD / dx;
end
if (shift_up)
    FU = FU ./ dx_star;
else
    FU = FU / dx;
end
p(i,j) = p(i,j) + k * dt / dx * (FR + FL + FU + FD);
end
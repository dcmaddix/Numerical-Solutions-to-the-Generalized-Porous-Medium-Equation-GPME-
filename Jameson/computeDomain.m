function x = computeDomain(dx,IC)
if (IC == 2) %need two ghost cells for periodic BC
    x=-2-2*dx:dx:2 +2*dx;%spatial vector
else %Dirichlet BC will fix the endpoints
    x = -1:dx:1;
end
end
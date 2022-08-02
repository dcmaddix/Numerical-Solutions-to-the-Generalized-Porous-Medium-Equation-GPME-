%verifying that the errors in 2D are symmetric and that there is an error
%in the plotting
tol = 1e-15;
dx = 1 / N;
xcoord = unique(x);
ycoord = unique(y);
for i = 1:length(xcoord)
    for j = 1:length(ycoord)
        coordx = xcoord(i);
        coordy = ycoord(j);
        [i1,j1] = find(x == coordx & y == coordy); %1st quadrant
        [i2,j2] = find(x == -coordx & y == coordy); %2nd quadrant
        [i3,j3] = find(x == -coordx & y == -coordy); %3rd quadrant
        [i4,j4] = find(x == coordx & y == -coordy); %4th quadrant
        if (max(abs(p(i1,j1) - p(i2,j2))) > tol)
            fprintf('anti-sym: 1st and 2nd quad\n')
        end
        if (max(abs(p(i2,j2) - p(i3,j3))) > tol)
            fprintf('anti-sym: 2nd and 3rd quad\n')
        end
        if (max(abs(p(i3,j3) - p(i4,j4))) > tol)
            fprintf('anti-sym: 3rd and 4th quad\n')
        end
        if (max(abs(p(i4,j4) - p(i1,j1))) > tol)
            fprintf('anti-sym: 1st and 4th quad\n')
        end
        if (max(abs(p(i1,j1) - p(j1,i1))) > tol)
            fprintf('anti-sym: x + y quad\n')
        end
    end
end
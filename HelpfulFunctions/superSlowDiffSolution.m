function u = superSlowDiffSolution(x, t)
    %Solve nonlinear system for w
    c = 2;
    w = x;
    for i = 1:length(x)
        w(i) = fsolve(@(w)(2+log(2*t)) * w + (c-w) * log(c-w) - ...
                    (c+w) * log(c+w) - abs(x(i)), 0.01);
    end
    v = pos(c^2 - w.^2) / (2*t)
    u = -1 ./ log(v)
    plot(u)
end
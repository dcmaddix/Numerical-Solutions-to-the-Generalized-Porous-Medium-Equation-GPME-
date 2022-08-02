function L = computeLimiter(limiter)
tol = 1e-6;
S = @(u,v) 0.5 * (sign(u) + sign(v)); % 0 if they have opposite signs
if (strcmp(limiter,'minmod'))
    L = @(u,v) S(u,v).*min(abs(u),abs(v));
elseif (strcmp(limiter,'vanleer'))
    L = @(u,v) 2*S(u,v).*abs(u).*abs(v) ./ (abs(u) + abs(v) + tol);
elseif (strcmp(limiter,'superbee'))
    L = @(u,v) S(u,v) .* max(min(2*abs(u),abs(v)), min(abs(u), 2*abs(v)));
else
    error('Invalid limiter');
end
end
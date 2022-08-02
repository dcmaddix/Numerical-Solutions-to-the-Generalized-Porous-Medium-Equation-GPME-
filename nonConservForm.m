%Solving nonconservative expanded form to high order: p_t + k_xp_x + kp_xx
%= 0 with advective and diffusive term for nonlinear k = p^3 with Dirichlet
%BC
function [p p_t t] = nonConservForm(bdry, ghostNodes, nt, p, dt, dx, x_coord, ...
                                    x_cent, k_deriv, eps_pxx, eps_px, ...
                                    eps_pxx2, eps_pxxx, k_pow, ramp, eps)
%INPUTS:
         
    %bdry - function handle for left an right bdry Dirichlet function 
            
    %ghostNodes - left and right ghost node for wider higher order stencil
    
    %nt - Final time
    
    %p - Initial Condition
    
    %dt - Timestep
    
    %dx - Spatial Step
    
    %x_coord - x position to store pressure at every time for pressure
                %versus time plot to check for oscillations
    
    %x_cent - position vector
    
    %k_deriv - discretization of k_x
    
    %eps_pxx - parameter for turning on or off anti-diffusion
    
    %eps_px - parameter for turning on or off advective error term
    
    %eps_pxx2 - parameter for turning on or off second order anti-diffusion
                %for first order schemes
              
%OUTPUTS:
    %p - pressure for plotting
    
    %p_t - pressure at specific x-coordinate
    
    %t - time vector
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %Plot pressure at this x-coordinate for each time to get pressure vs
    %time plot
    ind = find(x_cent < x_coord+1e-6 & x_cent > x_coord - 1e-6);
    if (length(ind) > 1)
        ind = ind(1);
    elseif (isempty(ind)) %make sure xcoord exists
        ind = 1;
        fprintf('Ind warning\n')
    end
    
    %compute time vector
    t = 0: dt: nt;
    p_t = zeros(1,length(t)); %time vector where initial condition is given at t=0
    %total number of gridpoints including the boundaries
    N = length(p(2:end-1)) + 2;
    %vector of interior gridpoints for narrow stencil
    i = 2:N-1;
    %power of k
    %k_pow = 3; %built in diffusion for htis power
    %Time loop: already have it at first time so less iteration
    m = 0;
    A = 0;
    B = 0;
    C = 0;
    if (k_pow == -1)
        p_star = 0.5;
        upper = p_star + eps; %band
        lower = p_star - eps;
        k_upper = 1.0;
        k_lower = 0.01;
        %linear case compute slope of line
        if (eps ~= 0)
            if (strcmp(ramp, 'linear'))
                m = (k_upper - k_lower) / (2 * eps); %upper - lower (p_star cancels)
                b = k_upper - m * upper;
            elseif(strcmp(ramp, 'arctan'))
                C = (k_upper + k_lower) / 2.0;
                A = (k_upper - C) * (4 / pi);
                B = 1 / eps;
            elseif(strcmp(ramp, 'sine'))
                C = (k_upper + k_lower) / 2.0;
                A = k_upper - C;
                B = pi / (2*eps);
            end
        end
    end
    k_p = zeros(length(i), 1);
    k_pp = zeros(length(i), 1);
    k_ppp = zeros(length(i), 1);
    for n = 1:length(t)-1
       %store at time
       p_t(n) = p(ind);
       if (k_pow ~= -1)
            k = (p.^k_pow)';
        else %Jakolein discontinous nonlinear k case
            if (eps ~= 0)
                if (strcmp(ramp, 'linear'))
                    k = m * p' + b;
                elseif(strcmp(ramp, 'arctan'))
                    k = A * atan(B * (p - 0.5))' + C;
                elseif(strcmp(ramp, 'sine'))
                    k = A * sin(B * (p - 0.5) )' + C;
                end
            end
            %Only define in ramp outside cap to values
            k(p <= lower) = k_lower;
            k(p > upper) = k_upper;
        end
       
       %Compute boundary terms
       p(1) = bdry{1}(t(n));
       p(end) = bdry{2}(t(n));
       
       ghostNodes_t = [ghostNodes{1}(t(n)) ghostNodes{2}(t(n)) ...
                       ghostNodes{1}(t(n)) ghostNodes{2}(t(n))];
       %Fourth Order First Derivative
       dp_dx = central_4thOrder(p, ghostNodes_t, dx);
       dp_dx = (p(i) - p(i-1)) / dx;
       
       %Fourth Order Second Derivative
       d2p_dx2 = central_2ndDeriv_4thOrder(p, ghostNodes_t, dx);
       %d2p_dx2 = (p(i+1) - 2* p(i) + p(i-1)) / dx^2;
       
       %Fourth Order Third Derivative
       d3p_dx3 = central_3rdDeriv_4thOrder(p, ghostNodes_t, dx);
       
       %Fourht Order Fourth Derivative
       d4p_dx4 = central_4thDeriv_4thOrder(p, ghostNodes_t, dx);
       
        if (strcmp(ramp, 'linear'))
                k_p = m * ones(length(i),1);
               elseif(strcmp(ramp, 'sine'))
                    k_p = A * B * cos(B * (p(i) - 0.5));
                    k_pp = -A * B^2 * sin(B * (p(i) - 0.5));
                    k_ppp = -A * B^3 * cos(B * (p(i) - 0.5));
                elseif (strcmp(ramp, 'arctan')) %cancel to avoid large inv eps
                    inv = (eps^2 + (p(i) - 0.5).^2).^(-1);
                    k_p = A * eps * inv;
                    k_pp = -2 * A * eps * (p(i) - 0.5) .* inv.^2;
                    k_ppp = 8 * A * eps * (p(i) - 0.5).^2 .* inv.^3 ...
                            - 2 * A * eps *inv.^2;
        end
         const_ind = p(i) < lower | p(i) > upper; %constant outside range
          k_p(const_ind) = 0.0; %not working to for p_xx
            k_pp(const_ind) = 0.0; %not working to for p_xx
            k_ppp(const_ind) = 0.0;%otherwise travel too far only apply locally
       
       %Assumes k_pow = 3 for tau and truncation arguments
       if (strcmp(k_deriv, 'Central_2nd')) %No oscillations:
           %O(dx^2) in space with positive diffusion, O(dx^2) in time with
           %some antidiffusion
           dk_dx = (k(i+1) - k(i-1))' / (2 * dx);
           %cancel positive diffusive term in mod eqtn to create
           %oscillations using tau
           if (k_pow ~= -1)
               tau = -eps_pxx * dx^2 * 3 * p(i) .* dp_dx.^2 .* d2p_dx2;
           else
               if (strcmp(ramp, 'linear'))
                k_p = m * ones(length(i), 1);
                k_pp = 0.0;
                k_ppp = 0.0;
               elseif(strcmp(ramp, 'sine'))
                    k_p = A * B * cos(B * (p(i) - 0.5));
                    k_pp = -A * B^2 * sin(B * (p(i) - 0.5));
                    k_ppp = -A * B^3 * cos(B * (p(i) - 0.5));
                elseif (strcmp(ramp, 'arctan')) %cancel to avoid large inv eps
                    inv = (eps^2 + (p(i) - 0.5).^2).^(-1);
                    k_p = A * eps * inv;
                    k_pp = -2 * A * eps * (p(i) - 0.5) .* inv.^2;
                    k_ppp = 8 * A * eps * (p(i) - 0.5).^2 .* inv.^3 ...
                            - 2 * A * eps *inv.^2;
               end
               %semidiscrete
               tau = -dx^2 /6 * (k_ppp .* dp_dx.^4 + ...
                       3 * k_pp .* dp_dx .* d2p_dx2 + ...
                       k_p .* dp_dx .* d3p_dx3);
%                  tau = tau + dp_dx.^2 .* d2p_dx2 .* ...
%                               (7 * dt / 2 * (k_p.^2 + k(i)' .* k_pp));
%                  tau = tau - d2p_dx2.^2 .* k_p .* (-2 * k(i)' * dt);
%                  tau = tau - d3p_dx3 .* dp_dx .* k_p .* (-3 * k(i)' * dt);
                  const_ind = p(i) < lower | p(i) > upper; %constant outside range
                  tau(const_ind) = 0.0; %not working to for p_xx
%                  tau = tau - d4p_dx4 .* k(i)' .* (-k(i)' * dt/ 2);
                    tau = 0.0;
           end
           %exactly cancels positive second order diffusion in modified
           %eqtn resulting in negative antidiffusion from p_tt: -105/2dtp^4p_x^2p_xx
           
       elseif (strcmp(k_deriv, 'Central_4th')) %No oscillations
           %O(dx^4) in space, O(dx^2) in time some some antidiff
           dk_dx = central_4thOrder(k, [ghostNodes{1}(t(n))^k_pow ...
                                    ghostNodes{2}(t(n))^k_pow], dx);
           %Add second order antidiffusion to create oscillations on
           %leading order in Mod Eqn truncation analysis
           tau = -eps_pxx * dx^2 * 3 * p(i) .* dp_dx.^2 .* d2p_dx2;
           tau = 0;
       elseif (strcmp(k_deriv,'D_plus')) %Oscillations
           %O(dx) in space with antidiffusion, O(dx^2) in time with some
           %antidiff
           dk_dx = (k(i+1) - k(i))' / dx;
           %Subtract term to cancel to remove oscillations:cancel first order terms
           if (k_pow ~= -1)
               tau = -eps_pxx * dx * (3/2) * p(i).^2 .* dp_dx .* d2p_dx2...
                     -eps_px * dx * 3 * p(i) .* dp_dx.^3;% ...
                 %+ dx^2 *3 * p(i) .* dp_dx.^2 .* d2p_dx2;
           else
               if (strcmp(ramp, 'linear'))
                k_p = m;
                k_pp = 0.0;
                k_ppp = 0.0;
               elseif(strcmp(ramp, 'sine'))
                    k_p = A * B * cos(B * (p(i) - 0.5));
                    k_pp = -A * B^2 * sin(B * (p(i) - 0.5));
                    k_ppp = -A * B^3 * cos(B * (p(i) - 0.5));
                elseif (strcmp(ramp, 'arctan')) %cancel to avoid large inv eps
                    inv = (eps^2 + (p(i) - 0.5).^2).^(-1);
                    k_p = A * eps * inv;
                    k_pp = -2 * A * eps * (p(i) - 0.5) .* inv.^2;
                    k_ppp = 8 * A * eps * (p(i) - 0.5).^2 .* inv.^3 ...
                            - 2 * A * eps *inv.^2;
               end
               %semidiscrete
               tau = -dx / 2 * (k_pp .* dp_dx.^3 + k_p .* d2p_dx2 .* dp_dx);
                tau = tau - dx^2 /6 * (eps_px * k_ppp .* dp_dx.^4 + ...
                       eps_pxx * 3 * k_pp .* dp_dx .* d2p_dx2 + ...
                       1 * k_p .* dp_dx .* d3p_dx3);
                 % tau = tau + dp_dx.^2 .* d2p_dx2 .* ...
                %               (7 * dt / 2 * (k_p.^2 + k(i)' .* k_pp));
                % tau = tau - d2p_dx2.^2 .* k_p .* (-2 * k(i)' * dt);
                 % tau = tau - d3p_dx3 .* dp_dx .* k_p .* (-3 * k(i)' * dt);
                  const_ind = p(i) < lower | p(i) > upper; %constant outside range
                  tau(const_ind) = 0.0; %not working to for p_xx
                % tau = tau - d4p_dx4 .* k(i)' .* (-k(i)' * dt/ 2);
           end
       elseif (strcmp(k_deriv, 'D_minus')) %No oscillations
           %O(dx) in space with diffusion, O(dx^2) in time with some
           %antidiff
           dk_dx = (k(i) - k(i-1))' / dx;
           %To make oscillations add antidiffusion
           if (k_pow ~= -1)
               tau =  eps_pxx * (dx * (3/2) * p(i).^2 .* dp_dx) .* d2p_dx2 ...
                    + eps_px * dx * 3*p(i).*dp_dx.^3  ...
                    - eps_pxx2 * dx^2 *3 * p(i) .* dp_dx.^2 .* d2p_dx2;
           else
               if (strcmp(ramp, 'linear'))
                k_p = m;
                k_pp = 0.0;
                k_ppp = 0.0;
               elseif(strcmp(ramp, 'sine'))
                    k_p = A * B * cos(B * (p(i) - 0.5));
                    k_pp = -A * B^2 * sin(B * (p(i) - 0.5));
                    k_ppp = -A * B^3 * cos(B * (p(i) - 0.5));
                elseif (strcmp(ramp, 'arctan')) %cancel to avoid large inv eps
                    inv = (eps^2 + (p(i) - 0.5).^2).^(-1);
                    k_p = A * eps * inv;
                    k_pp = -2 * A * eps * (p(i) - 0.5) .* inv.^2;
                    k_ppp = 8 * A * eps * (p(i) - 0.5).^2 .* inv.^3 ...
                            - 2 * A * eps *inv.^2;
               end
               %semidiscrete
               tau = dx / 2 * (k_pp .* dp_dx.^3 + k_p .* d2p_dx2 .* dp_dx);
               tau =  tau - dx^2 /6 * (k_ppp .* dp_dx.^4 + ...
                       3 * k_pp .* dp_dx .* d2p_dx2 + ...
                       1 * k_p .* dp_dx .* d3p_dx3);
                  %tau = tau + dp_dx.^2 .* d2p_dx2 .* ...
                  %             (7 * dt / 2 * (k_p.^2 + k(i)' .* k_pp));
                  %tau =  - d2p_dx2.^2 .* k_p .* (-2 * k(i)' * dt);
                 % tau = tau - d3p_dx3 .* dp_dx .* k_p .* (-3 * k(i)' * dt);
                  const_ind = p(i) < lower | p(i) > upper; %constant outside range
                  tau(const_ind) = 0.0; %not working to for p_xx
                  %tau = tau - d4p_dx4 .* k(i)' .* (-k(i)' * dt/ 2);
                  tau = 0;
           end
       elseif (strcmp(k_deriv, 'explicit')) %Oscillations
           %O(dx^4) in space, O(dx^2) in time
           dk_dx = k_pow*p(i).^(k_pow-1).*dp_dx; %explicitly compute
           dk_dx = k_p .* (p(i+1) - p(i)) / (dx);
           %Add leading order in truncation term dx^2 positive diffusion to
           %remove oscillations
           tau = eps_pxx * 3* dx^2  * p(i) .* dp_dx.^2 .* d2p_dx2; 
           tau = 0;
       else
           fprintf('Derivative Type Not Defined!\n')
           return;
       end
       %Time step update
       p(i) = p(i) + dt * (dk_dx .* dp_dx + k(i)'.*d2p_dx2 + tau);
       %try rk3 tvd
%        p1 = p;
%        p1(i) = p(i) + dt * (dk_dx .*dp_dx + k(i)'.*d2p_dx2);
%        dp_dx = central_4thOrder(p1, ghostNodes_t, dx);
%        
%        %Fourth Order Second Derivative
%        d2p_dx2 = central_2ndDeriv_4thOrder(p1, ghostNodes_t, dx);
%        p(i) = 0.5 * (p(i) + p1(i)) + 0.5*dt * (dk_dx .*dp_dx + k(i)'.*d2p_dx2) + dt*tau;
        %third order rk
%         p1 = p;
%         p1(i) = p(i) + dt * (dk_dx .*dp_dx + k(i)'.*d2p_dx2);
%         p2 = p;
%         dp_dx = central_4thOrder(p1, ghostNodes_t, dx);
%         d2p_dx2 = central_2ndDeriv_4thOrder(p1, ghostNodes_t, dx);
%         p2(i) = 3 / 4 * p(i) + 1/4 * p1(i) + dt /4 * (dk_dx .*dp_dx + k(i)'.*d2p_dx2);
%         dp_dx = central_4thOrder(p2, ghostNodes_t, dx);
%         d2p_dx2 = central_2ndDeriv_4thOrder(p2, ghostNodes_t, dx);
%         p(i) = 1/3*p(i) + 2/3*p2(i) + ...
%                2/3*dt *(dk_dx .*dp_dx + k(i)'.*d2p_dx2) + dt*tau;
%     end
    %compute boundaries at final time for Forward Euler
    p(1) = bdry{1}(t(end));
    p(end) = bdry{2}(t(end));
    p_t(end) = p(ind);
    end
end

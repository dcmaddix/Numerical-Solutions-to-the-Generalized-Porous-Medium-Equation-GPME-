%plot barenblatt solution
function [U, r] = barenblatt_pme(m,t,d, x,y, radial) %d is dimensions
   %paper assumes u_t = (u^m)_xx need to add +1 to m for our version
   m = m + 1;
   k = (m+1)^(-1);
   %zhang09 specific for 1D
   U = t^(-k)*(pos(1-k*(m-1)/(2*m)*abs(x).^2/(t^(2*k)))).^(1/(m-1));
   alpha = d / (d*(m-1)+2);
   C = 1; %arbitrty constant > 0
   beta = alpha / d;
   k = alpha * (m-1) / (2*m*d);
   if (radial)
       if (~radial)
           U = zeros(length(x),length(y));
       else
           U = zeros(length(x)*length(y),1);
       end
       count = 1;
       r = zeros(size(U));
       for i = 1:length(x)
           for j = 1:length(y)
               r(count) = norm([x(i,j) y(i,j)]);
               if (radial)
                U(count) = t^(-alpha)*(pos(C-k*r(count)^2*t^(-2*beta))).^(1/(m-1)); 
                count = count + 1;
               else
                  U(i,j) = t^(-alpha)*(pos(C-k*r(count)^2*t^(-2*beta))).^(1/(m-1));
               end
           end
       end
   end
   %U = t^(-alpha)*(pos(C-k*abs(x).^2*t^(-2*beta))).^(1/(m-1)); %vazquez dim independent form
   %interface at position where 0 between free boundary
   x_star = sqrt((2*m) / (k*(m-1)))*t^k;
   %in general solve for x in t = c|x|^(d(m-1)+2)
end
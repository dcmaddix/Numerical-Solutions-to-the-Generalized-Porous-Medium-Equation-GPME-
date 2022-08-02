function u = computeBC(u,IC) 
%this function computes the boudnary conditions
%second order stencil imples 4 bc u(1), u(2), u(end), u(end-1)
%statuonary shock case covered with Dirichlet BC fixed at initial
%position-unchanged
if (IC == 3) %fix leftside to be 1 to upper bound high of otherwise unbounded shock
   u(2) = 1; %u(1) = 1 remains unchanged same 
elseif (IC == 2) %periodic case for sin(pi*x) initil data
   u(1) = u(end - 3);%copy value from second rightmost unknown to ghost cell on left
   u(2) = u(end-2); %copy value from rightmost unknown to second ghost cell on left
   u(end) = u(4); %copy value from second leftmost unknown to last ghost cell on right
   u(end - 1) = u(3); %copy value from leftmost unknown to second to last ghost cell on right
end
end
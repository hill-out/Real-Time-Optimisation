function [finalC, dt, newC, newCost, newCons] = CSTRdyn(u, c0, stru, t)
% -------------------------------------------------------------------------
% runs the dynamic CSTR
%
% u     - [1x2]     - flowrates of A and B
% c0    - [1x4]     - initial reactor concentrations
% stru  - struct    - structure of parameters
% t     - double    - time step
%
% finalC    - [4x1]     - final concentration
% dt        - [nx1]     - time steps used by ode45
% newC      - [4xn]     - instant concentrations
% newCost   - [1xn]     - instant cost
% newCons   - [2xn]     - instant constraints
%
% -------------------------------------------------------------------------
[dt, newC] = ode45(@(t,c)(CSTRode(u, c, stru)),[0 t], c0);
newC = newC';
finalC = newC(:,end);

if nargout > 3
    newCost = stru.cost(u, newC);
    newCons = stru.cons([], newC);
end

end
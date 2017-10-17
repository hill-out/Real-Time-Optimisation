function [phi] = CSTR_cost(u, w, cIn, reactOrder, kVal, V)
% CSTR_cost takes the input for the CSTR and calc's the cost

uSize = size(u);
vol = u(1,:);

if uSize(1) == 1
    cSol = CSTR(reactOrder, kVal, cIn, u, V);
else
    cSol = u(2,:);
end

J = (cSol(3)^2*(vol(1)+vol(2))^2/(vol(1)*cIn(1))-w*(vol(1)^2+vol(2)^2));
phi = -J;

end
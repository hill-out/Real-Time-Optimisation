function [G] = CSTR_Qcon(u, dH, Qmax, cIn, reactOrder, kVal, V)
% CSTR_cost takes the input for the CSTR and calc's the cost

uSize = size(u);
vol = u(1,:);

if uSize(1) == 1
    cSol = CSTR(reactOrder, kVal, cIn, u, V);
else
    cSol = u(2,:);
end

Q = V*(kVal(1)*cSol(1)*cSol(2)*dH(1)+kVal(2)*cSol(2)^2*dH(2));
G = Q/Qmax-1;

end
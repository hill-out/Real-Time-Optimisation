function [G] = CSTR_Dcon(u, Dmax, cIn, reactOrder, kVal, V)
% CSTR_cost takes the input for the CSTR and calc's the cost

uSize = size(u);
vol = u(1,:);

if uSize(1) == 1
    if uSize(2) == 2
        u = [u(1),u(2),0,0];
    end
    
    cSol = CSTR(reactOrder, kVal, cIn, u, V);
else
    cSol = u(2,:);
end

D = cSol(4)/sum(cSol);
G = D/Dmax-1;

end
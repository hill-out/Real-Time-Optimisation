function phi = phi_from_Xu(X,u)
% Calculates the phi
% -----------------------------
% X         Mass fraction
% u         Inputs
% 
% phi       Cost
% -----------------------------

phi = -(1143.38*X(4)*(u(1)+u(2)) + 25.92*X(5)*(u(1)+u(2)) ...
    - 76.23*u(1) - 114.34*u(2));

end
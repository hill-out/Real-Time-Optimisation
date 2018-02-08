function g = g_from_Xu(X,u)
% Calculates the g
% -----------------------------
% X         Mass fraction
% u         Inputs
% 
% g         Constraints
% -----------------------------

g(1) = X(1) - 0.09;
g(2) = X(6) - 0.6;

end
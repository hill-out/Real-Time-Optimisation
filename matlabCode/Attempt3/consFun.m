function [cons] = consFun(~, c, stru)
% -------------------------------------------------------------------------
% calculates the cost
% 
% u     - [1x2]     - inlet flowrates
% c     - [1x4]     - concentrations
% stru  - struct    - structure of parameters
%
% cons  - double    - constraints at the current conditions
%
% -------------------------------------------------------------------------

% calculates the constraints (G) [2x1]
cons = zeros(2,size(c,2));

% constraint 1: Q/Qmax - 1 < 0
% Q = sum(r*dH*V)

cn = repmat(c,1,1,size(stru.k,2));
o = repmat(permute(stru.o,[1,3,2]),1,size(c,2),1);
baseM = permute(cn.^(o),[2,3,1]);
r = (repmat(stru.k,size(cn,2),1,1).*prod(baseM,3));

Q = sum(r.*repmat(stru.dH,size(r,1),1).*stru.V,2);
cons(1,:) = Q./stru.Qmax-1;

% constraint 1: D/Dmax - 1 < 0
D = c(4,:)./sum(c);
cons(2,:) = D./stru.Dmax - 1;

end
function r = reactionRate(c, stru)
% -------------------------------------------------------------------------
% calculates the rate of reaction for each reactant for a irrevesable CSTR
% 
% c         - [1x4]     - current concentrations
% stru      - struct    - structure of parameters
%
% r         - [1x4]     - rate of change of all substances
%
% -------------------------------------------------------------------------

% calculate the rate
c = repmat(c,1,size(stru.k,2));
baseM = c.^(stru.o);
r = (stru.k.*prod(baseM,1));

end
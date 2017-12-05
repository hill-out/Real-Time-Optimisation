function r = reactantRate(c, stru)
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
r = stru.s*reactionRate(c,stru)';

end
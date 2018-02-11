function [ J, G, xp ] = plant_test( y )
% y = [Xa (scaled), Fb]
% J is the cost function used by Marchetti (specified in DYCOPS paper)
% G is ??

global COST_PLANT

%plant function (simulated reality),
%CLOSED LOOP example: takes in the controlled variable!

solver_options = optimset('Display', 'off', 'Diagnostics', 'off','TolFun',1e-8,'MaxIter',1e3);

%controller equations:
Fb = y(2) + 2; %offset between setpoint value of Fb and value implemented
Fa = Fb/2.4; % Fa is proportional to Fb 
Xa = 0.75*y(1); %1.5*y(1);  % control error for Xa*100

% y = set points for [xa*100, Fb] xa must be between 0 and 0.27
% x = [Tr, xb, xc, xp, xe, xg]


x0 = [99.944 0.3813    0.0179    0.0977    0.2705    0.1126]; %optimal states
[xp,fval,exitflag] = fsolve(@(x)plantbalances(x, Xa, Fb, Fa),x0,solver_options);
if(exitflag ~= 1), display('fsolve terminated for wrong reason in plant'), J = Inf; G = [0,0]; return; end
if(any(xp(2:end)<0)), display('negative concentrations computed for the plant!'), keyboard, end

J = -( 1143.38*xp(4)*(Fb+Fa) + 25.92*xp(5)*(Fb+Fa) - 76.23*Fa - 114.34*Fb );
G = [-.6 + xp(6),... % xg <= .6
     Xa - 9]; % Xa*100 <= 9, or Xa <= .09 



COST_PLANT(end+1) = J; %save all plant costs for a realistic performance check

end

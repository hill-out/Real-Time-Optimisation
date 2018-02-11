function [ J, G,Gm_only, y ] = model(u)

global E_1 E_2;

solver_options = optimset('Display', 'off', 'Diagnostics', 'off','TolFun',1e-8,'MaxIter',1e3);

%constant
Fa = u(1); %kg/s
Fb = u(2);
% u = [Fa, Fb, Tr]
% y = [xa*100, Fb] xa must be between 0 and 0.27
% xm = [xa, xb, xp, xe, xg]


x0 = [0.08746 , 0.38962, 0.29061, 0.10945, 0.10754]; %optimal states, Marchetti
[xm,fval,exitflag] = fsolve(@(x)modelbalances_aprox(x,u,E_1, E_2),x0,solver_options);
if(exitflag ~= 1), display('fsolve terminated for wrong reason in model'), J = Inf; G = [0,0]; Gm_only = []; y = [0,0]; return; end

y = [xm(1)*100, u(2)]; %[xa*100, Fb]
J = -( 1143.38*xm(3)*(Fb+Fa) + 25.92*xm(4)*(Fb+Fa) - 76.23*Fa - 114.34*Fb );
G = [-.6 + xm(5),... % xg <= ..
    100*xm(1) - 9]; % Xa*100 <= 16
Gm_only = [];% J-1000; %model profit must be > ? ! Constraints that apply only to the model-based optimization problem
            %should not be active upon convergence!
end


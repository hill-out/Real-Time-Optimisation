clear all
close all
clc

global COST_PLANT E_1 E_2;

%% global uncertain model params, TO DELETE
E_1 =8100; %8100;  %8077.6;
E_2 = 12300;%12300; % %12438.5;

optprob.uL = [2, 8.5, 78];      %lower limit for u
optprob.uU = [15, 22, 92];    %upper limit for u
n_p = 10;
comp = 1;
for i=1:n_p+1,
    u1(i) = optprob.uL(1) + (i-1)* (optprob.uU(1)-optprob.uL(1))/n_p;
    for j = 1:n_p+1,
        u2(j) = optprob.uL(2) + (j-1)* (optprob.uU(2)-optprob.uL(2))/n_p;
        for k=1:n_p+1,
           u3(k) = optprob.uL(3) + (k-1)* (optprob.uU(3)-optprob.uL(3))/n_p;
           u = [u1(i), u2(j), u3(k)];
           [ J, G, Gm_only, y ] = model(u);
           Jc(i,j,k) = J;
           G1(i,j,k) = G(1);
           G2(i,j,k) = G(2);
           y1(i,j,k) = y(1);
           y2(i,j,k) = y(2);
           ureg1(comp) = y(1);
           ureg2(comp) = y(2);
           Jcv(comp) = J;
           G1cv(comp) = G(1);
           G2cv(comp) = G(2);
           comp = comp + 1;
        end
    end
end


%% optimal values
% model x = [xa, xb, xp, xe, xg]
% plant x = [Tr, xb, xc, xp, xe, xg]
% u = [Fa,Fb,Tr]
% y = [xa*100, Fb]
xstar_m = [0.1359    0.4121    0.0995    0.2604    0.6891]; %optimal states for model
ustar_m = [1.4260    3.3781   70.6995]; 
xstar_p = [100.5365    0.3897    0.0168    0.0980    0.2705    0.1118]; %optimal states for plant (unconstrained)
ystar_p = [11.3200   10.7456];
% Marchetti's optimal values for the ? [Fb, Tr] = [4.787, 89.70];

%% for MA
%optprob: optimization problem structure
n_u = 3;  %input dimension
optprob.u0 = [5, 10, 89.70]; %inital point for numerical optimization
optprob.cost = @model_cost; % cost function handle
optprob.constraint = @model_constraint; 
% optprob.constraint = @model_constraint;   % constraint function handle, NO
% CONSTRAINTS in this problem
optprob.uL = [2, 4, 60];      %lower limit for u
optprob.uU = [15, 30, 120];    %upper limit for u

%fmincon options
optprob.fmincon_options = optimset('Algorithm','interior-point','Display', 'on', 'Diagnostics', 'off', 'GradObj', 'off', ...
                    'GradConstr', 'off', 'DerivativeCheck', 'off', ... 
                    'LargeScale', 'off', 'Hessian', 'bfgs', ... 
                    'Diagnostics', 'on', 'TolX', 1e-8, ... 
                    'TolFun', 1e-8, 'TolCon', 1e-5, ... 
                    'MaxFunEval', 700, 'MaxIter', 500,'MaxSQPIter',80);
                
%MAoptions  
MAoptions.MAtype = 'Method B';
MAoptions.maxIT = 20;
MAoptions.uTol = 1e-3;
MAoptions.b = 0.3;
MAoptions.step_decrease = 0;
MAoptions.BFGS = 0;  %set this to 0 for now (Hessian modifiers under development..)


%plant function
MAplant = @plant;



%% optimize model
%myG = [];
myG = optprob.constraint;
[ustar, J_star_model,iout,fmincon_output,lagrange_mult,grad]=...
        fmincon( optprob.cost, ustar_m, [], [], [], [], optprob.uL, optprob.uU, myG, optprob.fmincon_options);
J_star_model
    
           [ J_starmodel, G_starmodel, Gm_only_model, y_star_model ] = model(ustar);
           J_starmodel
           y_star_model
           
refU = y_star_model;
param.refU = refU;
param.ureg1 = ureg1;
param.ureg2 = ureg2;
param.cvJ = Jcv;
param.cvJinit = J_starmodel;
param.eigmin1 = 0.03;
param.eigmin2 = 0.03;
param.testG = 0;

if param.testG == 0
    uxinit = [0 0 0 0 0];
    [uregJ] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
else
    uxinit = [0 0];
    [uregJ] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
end

[uregJ] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));


param.cvJ = G1cv;
param.cvJinit = G_starmodel(1);

param.testG = 0;
if param.testG == 0
    uxinit = [0 0 0 0 0];
    [uregG1] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
else
    uxinit = [0 0];
    [uregG1] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
    uregG1(3) = 0;
    uregG1(4) = 0;
    uregG1(5) = 0;
end
% 

param.cvJ = G2cv;
param.cvJinit = G_starmodel(2);
param.testG = 0;
if param.testG == 0
    uxinit = [0 0 0 0 0];
    [uregG2] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
else
    uxinit = [0 0];
    [uregG2] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
    uregG2(3) = 0;
    uregG2(4) = 0;
    uregG2(5) = 0;
end

uregJ
uregG1
uregG2
 

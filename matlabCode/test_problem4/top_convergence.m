clear all;
close all;
clc;
%convergence forcing plots for IFAC 2014

%% Williams Otto reactor from Marchetti2009Thesis (also in Marchetti2010 JPC paper), 
% It is assumed that Tr is adjusted by a controller to ensure a certain
% amount of x_a is present in the reactor. The 2 controlled variables (y) are ? 
% (y in the MAGeneric file refers to the SET POINT for the controlled variables)
% , the 3 manipulated variables (u) are ?. The mismatch is both
% parametric mismatch and control error in following the set points
% (implemented in the plant function). The closed loop MA converges to the
% optimal solution no problem. :)
% Notation should be standardised to be the same as in the DYCOPS paper.
% This version has constraints, it is an extension of the DYCOPS paper. Now
% there is a constraint on the concentration of xg. This is much more
% restrective for the model than the plant, with the result that the plant
% has no active constraints at the optimum, so it's not such a tricky
% constrained problem.. still the closed-loop adaptation with constraints
% still seems to work fine.

global COST_PLANT E_1 E_2;

%% global uncertain model params, TO DELETE
E_1 = 7900; %8100; %8077.6;
E_2 = 12500; %12300; %12438.5;

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
optprob.uL = [.5, 2.5, 60];      %lower limit for u
optprob.uU = [15, 30, 120];    %upper limit for u

%fmincon options
optprob.fmincon_options = optimset('Algorithm','interior-point','Display', 'on', 'Diagnostics', 'off', 'GradObj', 'off', ...
                    'GradConstr', 'off', 'DerivativeCheck', 'off', ... 
                    'LargeScale', 'off', 'Hessian', 'bfgs', ... 
                    'Diagnostics', 'on', 'TolX', 1e-8, ... 
                    'TolFun', 1e-8, 'TolCon', 1e-5, ... 
                    'MaxFunEval', 700, 'MaxIter', 500,'MaxSQPIter',80);
                
%MAoptions  
MAoptions.MAtype = 'Method A';
MAoptions.maxIT = 15;
MAoptions.uTol = 1e-3;
MAoptions.b = 1;
MAoptions.step_decrease = 0;
MAoptions.BFGS = 0;  %set this to 0 for now (Hessian modifiers under development..)

%plant function
MAplant = @plant;



%% optimize model
%myG = [];
myG = optprob.constraint;
[ustar, J_star_model,iout,fmincon_output,lagrange_mult,grad]=...
        fmincon( optprob.cost, ustar_m, [], [], [], [], optprob.uL, optprob.uU, myG, optprob.fmincon_options);
    ustar
    J_star_model

    keyboard
%% optimize plant
 
 y0 = [5,10];
 yL = [6,1]; % [xa*100, Fb]
 yU = [16,25];
 
[ystar, J_star_plant ,iout,fmincon_output,lagrange_mult,grad]=...
        fmincon( @plant_cost, y0, [], [], [], [], yL, yU, @plant_constraint, optprob.fmincon_options);
    ystar
    J_star_plant
    
    

%% plot plant cost function and constraint
%generate points
xa_vals = 0.1:2:20;
Fb_vals = 1:2:32;
[xa_grid,Fb_grid] = meshgrid(xa_vals, Fb_vals);
[m,n] = size(xa_grid);
%
for i = 1:m
    for j = 1:n
        y = [xa_grid(i,j), Fb_grid(i,j)];
        J(i,j) = plant_cost( y );
    end
end
% plot the cost
%%
conthandle = figure,
J(J==inf)=NaN;
[C,h] = contour(xa_grid,Fb_grid,-J,10); %plot plant cost function
set(h,'ShowText','on','LabelSpacing',3000) %,'TextListMode','auto','TextStep',get(h,'LevelStep')*4);
dummy = -100:20:240;
set(h,'LevelList',dummy);
set(h,'TextList', dummy(1:3:end))
title('J'), xlabel('xa (sp)'), ylabel('Fb (sp)')
hold on; 
%plot([6,6],[min(Fb_vals),max(Fb_vals)],'k--') %plot active constraint for plant
%
%plot the constraint function
% for i = 1:length(Fb_vals)
%     xa_constraint(i) = fsolve(@(xa)plant_constraint([xa,Fb_vals(i)]), 10);
% end
% plot(xa_constraint, Fb_vals, 'k')
%
% figure,
% [C,h] = contour(X,Y,10000*G,10); %plot plant constraint function
% colorbar;

%% MA single case
% COST_PLANT = []; %reset
% MAoutput{1} = GenericMA( optprob, MAplant, MAoptions, conthandle );
% display(MAoutput{1}.term_reason);
%%
% costhandle = figure
% hold on
% plot(-COST_PLANT,'r')
%%

%% testing GenericMA
COST_PLANT = []; %reset
E_1vals = [7900];
E_2vals = [12300];
mycolor{1} = {'b','r','g'};
mycolor{2} = {'r','r--','g--'};
%%
costhandle = figure;
step_decrease = [0.5, 0];

for j=1:2
    MAoptions.step_decrease = step_decrease(j);
    
    for i=1:length(E_1vals) %3 cases
        E_1 = E_1vals(i);
        E_2 = E_2vals(i);
        COST_PLANT = []; %reset
        MAoutput{j}{i} = GenericMA( optprob, MAplant, MAoptions, conthandle );
        display(MAoutput{j}{i}.term_reason);

        figure(costhandle)
        hold on
        plot(-COST_PLANT,mycolor{j}{i})
    end

end


%% plotting results

cost2handle = figure;
constrainthandle = figure;
% figure
% plot(MAoutput.y_vals)

for j=1:2

    for i=1:length(E_1vals)
        figure(cost2handle)
        hold on
        plot(1:length(MAoutput{j}{i}.cost_vals),-MAoutput{j}{i}.cost_vals,mycolor{j}{i})
        xlabel('k'), ylabel('\Phi_p')
        %
        figure(conthandle)
        hold on
        plot(MAoutput{j}{i}.y_vals(:,1),MAoutput{j}{i}.y_vals(:,2),mycolor{j}{i},'LineWidth',2)


        figure(constrainthandle)
        hold on
        plot(1:length(MAoutput{j}{i}.cost_vals),MAoutput{j}{i}.constraint_vals, mycolor{j}{i})
        xlabel('k'), ylabel('constraints')
        %legend('1','2')
    end

end
%     for i = 1:(length(MAoutput.cost_vals)-1)
%         y = MAoutput.y_vals(i,:);
%         y_next = MAoutput.y_vals(i+1,:);
%         arrow(y,y_next,10,'BaseAngle',60)
%     end





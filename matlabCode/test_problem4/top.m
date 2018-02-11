clear all;
close all;
clc;

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
E_1 = 8100; %8100; %8077.6;
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
optprob.uL = [2, 4, 60];      %lower limit for u
optprob.uU = [15, 30, 120];    %upper limit for u

%fmincon options
optprob.fmincon_options = optimset('Algorithm','interior-point','Display', 'off', 'Diagnostics', 'off', 'GradObj', 'off', ...
                    'GradConstr', 'off', 'DerivativeCheck', 'off', ... 
                    'LargeScale', 'off', 'Hessian', 'bfgs', ... 
                    'Diagnostics', 'on', 'TolX', 1e-8, ... 
                    'TolFun', 1e-8, 'TolCon', 1e-5, ... 
                    'MaxFunEval', 700, 'MaxIter', 500,'MaxSQPIter',80);
                
%MAoptions  
MAoptions.MAtype = 'Method B';
MAoptions.maxIT = 30;
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
    ustar
    J_star_model

    %keyboard
%% optimize plant
 
 y0 = [5,10];
 yL = [6,1]; % [xa*100, Fb]
 yU = [20,25];
 
[ystar, J_star_plant ,iout,fmincon_output,lagrange_mult,grad]=...
        fmincon( @plant_cost, y0, [], [], [], [], yL, yU, @plant_constraint, optprob.fmincon_options);
    ystar
    J_star_plant
    
    

%% plot plant cost function and constraint
%generate points
xa_vals = 6:1:18;
Fb_vals = 1:1:25;
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
conthandle = figure('name','profit contour');
hold on
J(J==inf)=NaN;
myp = fill([12,18,18,12],[0,0,25,25],.8*[1 1 1]); %,'FaceAlpha',.5);
grid off
[C,h] = contour(xa_grid,Fb_grid,-J,10); %plot plant cost function
%colormap(gray);
set(h,'ShowText','on','LabelSpacing',3000) %,'TextListMode','auto','TextStep',get(h,'LevelStep')*4);
dummy = -350:40:240;
set(h,'LevelList',dummy);
set(h,'TextList', dummy(1:3:end))
xlabel('X_{A,\rm{s}} x 100 (-)'), ylabel('F_{B,\rm{s}} (kg.s^{-1})')
hold on; 
%plot([12,12],[min(Fb_vals),max(Fb_vals)],'k--') %plot active constraint for plant

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
% %%
% costhandle = figure
% hold on
% plot(-COST_PLANT,'r')
%%

%% testing GenericMA
MAoptions.epsilon0 = [0, 0]; % initial backoff
COST_PLANT = []; %reset
E_1vals = [8050 8100 8100]; % 8050
E_2vals = [12500 12500 12300]; % 12500
mycolor{1} = {'b','r','g'};
mycolor{2} = {'b--','r--','g--'};

costhandle = figure;

method = {'Method A', 'Method B'};
for j=1:2
    MAoptions.MAtype = method{j};
    
    for i=1:3 %3 cases
        E_1 = E_1vals(i);
        E_2 = E_2vals(i);
        COST_PLANT = []; %reset
        MAoutput{j}{i} = GenericMA( optprob, MAplant, MAoptions, conthandle );
        display(MAoutput{j}{i}.term_reason);

        figure(costhandle)
        hold on
        plot(-COST_PLANT,mycolor{j}{i})
        %keyboard
    end

end
%% convex MA
for compt = 1:3,
E_1 = E_1vals(compt); %8100; %8077.6;
E_2 = E_2vals(compt); %12300; %12438.5;

% optprob.uL = [2, 8.5, 80];      %lower limit for u
% optprob.uU = [15, 20, 95];    %upper limit for u

optprob.uL = [2, 8.5, 78];      %lower limit for u
optprob.uU = [10, 22, 92];    %upper limit for u

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

xstar_m = [0.1359    0.4121    0.0995    0.2604    0.6891]; %optimal states for model
ustar_m = [1.4260    3.3781   70.6995]; 
xstar_p = [100.5365    0.3897    0.0168    0.0980    0.2705    0.1118]; %optimal states for plant (unconstrained)
ystar_p = [11.3200   10.7456];
% Marchetti's optimal values for the ? [Fb, Tr] = [4.787, 89.70];

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
%MAoptions.maxIT = 30;
MAoptions.uTol = 1e-3;
MAoptions.b = 0.3;
MAoptions.step_decrease = 0;
MAoptions.BFGS = 0;  %set this to 0 for now (Hessian modifiers under development..)

%plant function
MAplant = @plant;

myG = optprob.constraint;
[ustar, J_star_model,iout,fmincon_output,lagrange_mult,grad]=...
        fmincon( optprob.cost, ustar_m, [], [], [], [], optprob.uL, optprob.uU, myG, optprob.fmincon_options);

    
           [ J_smodel, G_smodel, Gm_o_model, y_s_model] = model(ustar);
            J_starmodel{compt}{1} = J_smodel;
            G_starmodel{compt}{1} = G_smodel(1);
            G_starmodel{compt}{2} = G_smodel(2);
            for k=1:2,
                y_starmodel{compt}{k} = y_s_model(k);
            end
                
           
refU = [y_starmodel{compt}{1};y_starmodel{compt}{2}];
param.refU = refU;
param.ureg1 = ureg1;
param.ureg2 = ureg2;
param.cvJ = Jcv;
param.cvJinit = J_starmodel{compt}{1};
param.eigmin1 = 0.05;
param.eigmin2 = 0.05;
param.testG = 0;

if param.testG == 0
    uxinit = [0 0 0 0 0];
    [buff] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
    for k = 1:5
        uregJ{compt}{k} = buff(k);
    end
else
    uxinit = [0 0];
    [buff] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
    for k = 1:5
        if k <=2
            uregJ{compt}{k} = buff(k);
        else
            uregJ{compt}{k} = 0;
        end
    end
    
end

eig([uregJ{compt}{3} uregJ{compt}{4}; uregJ{compt}{4} uregJ{compt}{5}])

param.cvJ = G1cv;
param.cvJinit = G_starmodel{compt}{1};

param.testG = 1;

if param.testG == 0
    uxinit = [0 0 0 0 0];
    [buff] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
   for k = 1:5
       uregG1{compt}{k} = buff(k);
   end

else
    uxinit = [0 0];
    [buff] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
    for k = 1:5
        if k <=2
            uregG1{compt}{k} = buff(k);
        else
            uregG1{compt}{k} = 0;
        end
    end
end

%eig([uregG1{compt}{3} uregG1{compt}{4}; uregG1{compt}{4} uregG1{compt}{5}])

param.cvJ = G2cv;
param.cvJinit = G_starmodel{compt}{2};
param.testG = 1;
if param.testG == 0
    uxinit = [0 0 0 0 0];
    [buff] = fmincon(@(u) costregconv(u, param),uxinit,[],[],[],[],[],[],@(u) constrregconv(u, param),optimset('Display','on','Algorithm','interior-point'));
        for k = 1:5
            uregG2{compt}{k} = buff(k);
        end
else
    uxinit = [0 0];
    [buff] = fminunc(@(u) costregconv(u, param),uxinit,optimset('Display','on','Algorithm','interior-point'));
    for k = 1:5
        if k <=2
            uregG2{compt}{k} = buff(k);
        else
            uregG2{compt}{k} = 0;
        end
    end

end
%eig([uregG2{compt}{3} uregG2{compt}{4}; uregG2{compt}{4} uregG2{compt}{5}])


end

K_filter = 0.85;
K_filter2 = 0.5;
optprob.cL = [1, 3];      %lower limit for u
optprob.cU = [20, 25];    %upper limit for u
for j = 1:3,
    for k = 1:2,
        % model inverse y_1 et y_2 par rapport a G1 et G2  
%         if k+1 == 2
%             bidon = k+1;
%         else
%             bidon = k-1;
%         end
        cs{j}{1}(k) = y_starmodel{j}{k};
    end
end

    
for j = 1:3,
    E_1 = E_1vals(j);
    E_2 = E_2vals(j);
%     for k = 1:2
%         ymo(k) = cs{j}{k};
%     end
%     [Jp_MAc{j}{1}, Gp_MA] = plant(ymo);
%     
%     for k =1:2
%         Gp_MAc{j}{k} = Gp_MA(k);
%     end
    
    % yref = [cs{j}{1}(1)/0.75;cs{j}{1}(2)-2];

    param.alp = [uregJ{j}{1} uregJ{j}{2} uregJ{j}{3} uregJ{j}{4} uregJ{j}{5}];
    param.alp1 = [uregG1{j}{1} uregG1{j}{2} uregG1{j}{3} uregG1{j}{4} uregG1{j}{5}];
    param.alp2 = [uregG2{j}{1} uregG2{j}{2} uregG2{j}{3} uregG2{j}{4} uregG2{j}{5}];
    yref = [cs{j}{1}(1);cs{j}{1}(2)];
    
    param.cvJinit = J_starmodel{j}{1};
    param.G1init = G_starmodel{j}{1};
    param.G2init = G_starmodel{j}{2};   
    refUinit = yref;
    param.refUinit = refUinit;
    param.refU = yref;
    RTOiter = MAoptions.maxIT;
for i = 1:RTOiter,
   
%     if i==1,
%         Jmodcv{j}{i} = J_starmodel{j}{1};
%         G1modcv{j}{i} = G_starmodel{j}{1};
%         G2modcv{j}{i} = G_starmodel{j}{2};
%         gradJmodcv{j}{i} = [uregJ{j}{1}; uregJ{j}{2}];
%         gradG1modcv{j}{i} = [uregG1{j}{1}; uregG1{j}{2}];
%         gradG2modcv{j}{i} = [uregG2{j}{1}; uregG2{j}{2}];
%     else
        Jmodcv{j}{i} = J_starmodel{j}{1}  + [uregJ{j}{1}; uregJ{j}{2}]'*[cs{j}{i}'-refUinit]+0.5*[cs{j}{i}'-refUinit]'*[uregJ{j}{3} uregJ{j}{4};uregJ{j}{4} uregJ{j}{5}]*[cs{j}{i}'-refUinit];
        G1modcv{j}{i} = G_starmodel{j}{1}  + [uregG1{j}{1}; uregG1{j}{2}]'*[cs{j}{i}'-refUinit]+0.5*[cs{j}{i}'-refUinit]'*[uregG1{j}{3} uregG1{j}{4};uregG1{j}{4} uregG1{j}{5}]*[cs{j}{i}'-refUinit];
        G2modcv{j}{i} = G_starmodel{j}{2}  + [uregG2{j}{1}; uregG2{j}{2}]'*[cs{j}{i}'-refUinit]+0.5*[cs{j}{i}'-refUinit]'*[uregG2{j}{3} uregG2{j}{4};uregG2{j}{4} uregG2{j}{5}]*[cs{j}{i}'-refUinit];
        gradJmodcv{j}{i} = [uregJ{j}{1}; uregJ{j}{2}] + [uregJ{j}{3} uregJ{j}{4};uregJ{j}{4} uregJ{j}{5}]*[cs{j}{i}'-refUinit];
        gradG1modcv{j}{i} = [uregG1{j}{1}; uregG1{j}{2}] + [uregG1{j}{3} uregG1{j}{4};uregG1{j}{4} uregG1{j}{5}]*[cs{j}{i}'-refUinit];
        gradG2modcv{j}{i} = [uregG2{j}{1}; uregG2{j}{2}] + [uregG2{j}{3} uregG2{j}{4};uregG2{j}{4} uregG2{j}{5}]*[cs{j}{i}'-refUinit];
%     end
    
    [Jp, Gp] = plant(cs{j}{i}');
    JplantMAcv{j}{i} = Jp;
    mJplantMAcv{j}{i} = -Jp;
    G1plantMAcv{j}{i} = Gp(1);
    G2plantMAcv{j}{i} = Gp(2);
    delt = 0.0000001;
    cs1plus =  [(1+delt)*cs{j}{i}(1); cs{j}{i}(2)];
    cs2plus =  [cs{j}{i}(1); (1+delt)*cs{j}{i}(2)];
    [Jpplus1, Gpplus1] = plant(cs1plus');
    [Jpplus2, Gpplus2] = plant(cs2plus');
%    [grad1] = finite_diff(cs{j}{i}', @plant_cost,del);
    gradJp{j}{i} = [(Jpplus1-Jp)/(delt*cs{j}{i}(1));(Jpplus2-Jp)/(delt*cs{j}{i}(2))];
    col1gradGp = [(Gpplus1(1)-Gp(1))/(delt*cs{j}{i}(1));(Gpplus2(1)-Gp(1))/(delt*cs{j}{i}(2))];
    col2gradGp = [(Gpplus1(2)-Gp(2))/(delt*cs{j}{i}(1));(Gpplus2(2)-Gp(2))/(delt*cs{j}{i}(2))];
    %    [grad2] = finite_diff(cs{j}{i}', @plant_constraint,del);
    %gradGp{j}{i} = grad2;
    gradGp{j}{i} = [col1gradGp col2gradGp];
%    epsJ{j}{i} = JplantMAcv{j}{i} - Jmodcv{j}{i}; 
    if i==1
        epsJ{j}{i} = K_filter2*(JplantMAcv{j}{i} - Jmodcv{j}{i});
        epsG{j}{i} = K_filter2*([G1plantMAcv{j}{i} - G1modcv{j}{i}; G2plantMAcv{j}{i} - G2modcv{j}{i}]);
        lambdaJ{j}{i} = K_filter2*(gradJp{j}{i} - gradJmodcv{j}{i});
        lambdaG{j}{i} = K_filter2*(gradGp{j}{i} - [gradG1modcv{j}{i} gradG2modcv{j}{i}]);
    else
        epsJ{j}{i} = K_filter*(JplantMAcv{j}{i} - Jmodcv{j}{i})+(1-K_filter)*epsJ{j}{i-1};
        epsG{j}{i} = [K_filter 0;0 K_filter2]*[G1plantMAcv{j}{i} - G1modcv{j}{i}; G2plantMAcv{j}{i} - G2modcv{j}{i}]+([1-K_filter 0;0 1-K_filter2])*epsG{j}{i-1};
        lambdaJ{j}{i} = K_filter*(gradJp{j}{i} - gradJmodcv{j}{i})+(1-K_filter)*lambdaJ{j}{i-1};
        lambdaG{j}{i} = K_filter*(gradGp{j}{i} - [gradG1modcv{j}{i} gradG2modcv{j}{i}])+(1-K_filter)*lambdaG{j}{i-1};
    end
    myG = @constraintscv;
%     if i>1,
%     param.epscost = K_filter*epsJ{j}{i}+(1-K_filter)*epsJ{j}{i-1};
%     param.epsconst = K_filter*epsG{j}{i}+(1-K_filter)*epsG{j}{i-1};
%     param.lambdacost = K_filter*lambdaJ{j}{i}+(1-K_filter)*lambdaJ{j}{i-1};
%     param.lambdaconst =K_filter*lambdaG{j}{i}+(1-K_filter)*lambdaG{j}{i-1};
%     param.cvJinit = J_starmodel{j}{1};
%     param.G1init = G_starmodel{j}{1};
%     param.G2init = G_starmodel{j}{2};
%     else
     param.epscost = epsJ{j}{i};
     param.epsconst1 = epsG{j}{i}(1);
     param.epsconst2 = epsG{j}{i}(2);
     param.lambdacost = lambdaJ{j}{i};
     if i==1,
        param.lambdaconst1 = col1gradGp - gradG1modcv{j}{i};
        param.lambdaconst2 = col2gradGp - gradG2modcv{j}{i};    
     else
        param.lambdaconst1 = K_filter*(gradGp{j}{i}(:,1) - gradG1modcv{j}{i}) + (1-K_filter)*(gradGp{j}{i-1}(:,1) - gradG1modcv{j}{i-1});
        param.lambdaconst2 = K_filter*(gradGp{j}{i}(:,2) - gradG2modcv{j}{i}) + (1-K_filter)*(gradGp{j}{i-1}(:,2) - gradG2modcv{j}{i-1});
     end

%     end
    
    ustar_m = cs{j}{i};
    [ustar]= fmincon(@modelcostcv, ustar_m, [], [], [], [], optprob.cL, optprob.cU, myG, optprob.fmincon_options, param);
    ustar
    cs{j}{i+1} = ustar;
    % K_filter*ustar+(1-K_filter)*cs{j}{i};
    yref = [cs{j}{i}(1);cs{j}{i}(2)];
    param.refU = yref;

end
end


%% plotting results

cost2handle = figure;
constrainthandle = figure;
constrainthandle2 = figure;

% figure
% plot(MAoutput.y_vals)
for k=1:RTOiter,
    comptRTO(k) = k;
    RTOcostpt1(k) = mJplantMAcv{1}{k};
    RTOcostpt2(k) = mJplantMAcv{2}{k};
    RTOcostpt3(k) = mJplantMAcv{3}{k};
    RTOXGpt1(k) = G1plantMAcv{1}{k};
    RTOXGpt2(k) = G1plantMAcv{2}{k};
    RTOXGpt3(k) = G1plantMAcv{3}{k};
    RTOXApt1(k) = G2plantMAcv{1}{k};
    RTOXApt2(k) = G2plantMAcv{2}{k};
    RTOXApt3(k) = G2plantMAcv{3}{k};
    u1pt1(k) = cs{1}{k}(1);
    u2pt1(k) = cs{1}{k}(2);
    u1pt2(k) = cs{2}{k}(1);
    u2pt2(k) = cs{2}{k}(2);
    u1pt3(k) = cs{3}{k}(1);
    u2pt3(k) = cs{3}{k}(2);
end
for j=1:2

    for i=1:3
        figure(cost2handle)
        hold on
        grid off
        plot(1:length(MAoutput{j}{i}.cost_vals),-MAoutput{j}{i}.cost_vals,mycolor{j}{i},'LineWidth',2)
        plot([1,RTOiter], 210*[1,1],'k:','LineWidth',2)
%         if j == 2,
%             if i==1,
%                 figure(cost2handle)
%                 hold on
%                 grid off
%                 plot(comptRTO,RTOcostpt1,'b-.','LineWidth',2)
%             elseif i==2,
%                 figure(cost2handle)
%                 hold on
%                 grid off
%                 plot(comptRTO,RTOcostpt2,'r-.','LineWidth',2)
%             elseif i==3,
%                 figure(cost2handle)
%                 hold on
%                 grid off
%                 plot(comptRTO,RTOcostpt3,'g-.','LineWidth',2)
%             end    
%         end
        xlabel('RTO iteration, k (-)'), ylabel('Profit, -\Phi_{\rm{p}} (-)')
        if i == 3
        end
        figure(conthandle)
        grid off
        hold on
        plot(MAoutput{j}{i}.y_vals(:,1),MAoutput{j}{i}.y_vals(:,2),mycolor{j}{i},'LineWidth',2)
        set(gca,'fontsize',20)
        
        x_vec = 1:length(MAoutput{j}{i}.cost_vals);
        figure(constrainthandle)
        grid off
        hold on
        plot(x_vec,MAoutput{j}{i}.constraint_vals(:,1), mycolor{j}{i},'LineWidth',2)
        set(gca,'fontsize',20)
        plot(x_vec, zeros(size(x_vec)),'k--')
        set(gca,'fontsize',20)
        xlabel('RTO iteration, k (-)'), ylabel('X_G^{\rm{max}}-X_G (-)')
        
        
        figure(constrainthandle2)
        hold on
        grid off
        plot(1:length(MAoutput{j}{i}.cost_vals),MAoutput{j}{i}.constraint_vals(:,2), mycolor{j}{i},'LineWidth',2)
        set(gca,'fontsize',20)
        plot(x_vec, zeros(size(x_vec)),'k--')
        set(gca,'fontsize',20)
        xlabel('RTO iteration, k'), ylabel('X_A^{\rm{max}}-X_A (-)')
        %legend('1','2')
        
         if j == 2,
            if i==1,
                figure(cost2handle)
                hold on
                grid off
                plot(comptRTO,RTOcostpt1,'b-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(conthandle)
                grid off
                hold on
                plot(u1pt1,u2pt1,'b-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle)
                grid off
                hold on
                plot(comptRTO,RTOXGpt1,'b-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle2)
                hold on
                grid off
                plot(comptRTO,RTOXApt1,'b-.','LineWidth',2)
            elseif i==2,
                figure(cost2handle)
                hold on
                grid off
                plot(comptRTO,RTOcostpt2,'r-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(conthandle)
                grid off
                hold on
                plot(u1pt2,u2pt2,'r-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle)
                grid off
                hold on
                plot(comptRTO,RTOXGpt2,'r-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle2)
                hold on
                grid off
                plot(comptRTO,RTOXApt2,'r-.','LineWidth',2)
                set(gca,'fontsize',20)
            elseif i==3,
                figure(cost2handle)
                hold on
                grid off
                plot(comptRTO,RTOcostpt3,'g-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(conthandle)
                grid off
                hold on
                plot(u1pt3,u2pt3,'g-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle)
                grid off
                hold on
                plot(comptRTO,RTOXGpt3,'g-.','LineWidth',2)
                set(gca,'fontsize',20)
                figure(constrainthandle2)
                hold on
                grid off
                plot(comptRTO,RTOXApt3,'g-.','LineWidth',2)
                set(gca,'fontsize',20)
            end    
        end
    end
end

% plot letters A,B,C for each case on contour plot
mystr = {'I','II','II'};
for i = 1:3
    figure(conthandle)
    hold on
    plot(ystar(1),ystar(2),'ko','MarkerSize',8,'MarkerFaceColor','k') %plot a circle at the plant optimum%
    myt = text(MAoutput{1}{i}.y_vals(1,1)-.3,MAoutput{1}{i}.y_vals(1,2)...
        ,mystr{i}, 'Color',mycolor{1}{i});
end

%%
% break
set(gcf, 'Position', [1200 500 450 300])
axis tight
%%
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',25)
%     for i = 1:(length(MAoutput.cost_vals)-1)
%         y = MAoutput.y_vals(i,:);
%         y_next = MAoutput.y_vals(i+1,:);
%         arrow(y,y_next,10,'BaseAngle',60)
%     end





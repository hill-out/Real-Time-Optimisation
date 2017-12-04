function [] = dynamicRTO(modelFun,plantFun,gain,u0,c0,tau, method)
% Takes the modelFun (either convex or other) and the plantFun and uses a
% standard multiple unit MA to work out the optimum position of the plant

if nargin < 7 ||isempty(method)
    method = 'MU';
end

if strcmp(method,'MU') % multiple units
    meth = 1;
elseif strcmp(method,'NE') %neighbouring extremals
    meth = 2;
else % incorrect input
    method = 'MU';
    meth = 1;
end

nTau = 0.01; %number of minutes for the time step in the dynamics CSTR
n = tau/nTau; %number of time steps

% Base Case
delta = 50e-3; %change in the u value

if meth == 1
    [base, other, gradJ, gradG] = gradMU(plantFun, c0, repmat(c0,2,1), u0, n, tau, delta);
elseif meth == 2
    [base, other, gradJ, gradG] = gradNE(plantFun, c0, repmat(c0,2,1), u0, 10, 0.1, delta);
else
    error('no method')
end

% Model Steady State
m = cellFuncHandle2vec(modelFun,{[u0(1),u0(2)]},ones(size(modelFun))')';
mGrad = zeros(1,2,3);

e(1) = base.J(end) - m(1);
e(2:3) = base.G(end,:) - m(2:3);
e = e;
l = [permute(gradJ,[3,2,1]),permute(gradG,[3,2,1])]-permute(mGrad(1,:,:),[2,3,1]);

modifiedFun{1} = @(v)(modelFun{1}(v)+e(1,1)+(v-u0)*(l(:,1)));
modifiedFun{2} = @(v)(modelFun{2}(v)+e(1,2)+(v-u0)*(l(:,2)));
modifiedFun{3} = @(v)(modelFun{3}(v)+e(1,3)+(v-u0)*(l(:,3)));

first = 1;
solved = 0;
u = u0;
j = 2;

while solved == 0 || first == 1
    
    u(j,:) = CSTR_Opt(u(j-1,:), modifiedFun{1}, {modifiedFun{2:end}});
    
    % calc the current plant gradient

    if meth == 1
        [base, other, gradJ, gradG] = gradMU(plantFun, [], [], u(j,:), n, tau, delta, base, other, gradJ, gradG);
    elseif meth == 2
        [base, other, gradJ, gradG] = gradNE(plantFun, [], [], u(j,:), n, tau, delta, base, other, gradJ, gradG);
    else
        error('no method')
    end
    
    % Model Steady State
    m(j,:) = cellFuncHandle2vec(modelFun,{u(j,:)},ones(size(modelFun))')';
    mGrad(j,:,:) = dfunc(modelFun,u(j,:));
    
    newe(1) = base.J(end) - m(j,1);
    newe(2:3) = base.G(end,:) - m(j,2:3);
    newl = [permute(gradJ(end,:,:),[3,2,1]),permute(gradG(end,:,:),[3,2,1])]-permute(mGrad(j,:,:),[2,3,1]);
    
    e(j,1) = newe(1).*gain*0.1 + (1-gain*0.1).*e(j-1,1);
    e(j,2:3) = newe(2:3).*gain + (1-gain).*e(j-1,2:3);
    l(:,1,j) = newl(:,1).*gain*0.1 + (1-gain*0.1).*l(:,1,j-1);
    l(:,2:3,j) = newl(:,2:3).*gain + (1-gain).*l(:,2:3,j-1);
    
    modifiedFun{1} = @(v)(modelFun{1}(v)+e(j,1)+(v-u(j,:))*(l(:,1,j)));
    modifiedFun{2} = @(v)(modelFun{2}(v)+e(j,2)+(v-u(j,:))*(l(:,2,j)));
    modifiedFun{3} = @(v)(modelFun{3}(v)+e(j,3)+(v-u(j,:))*(l(:,3,j)));
    
    first = 0;
    
    converged = all(all([(abs((e(j,:) - e(j-1,:))./e(j-1,:)) < 1e-3);(abs((l(:,:,j,:) - l(:,:,j-1))./l(:,:,j-1)) < 1e-3)]));
    failed = any(u(j,1:2)<1);
    
    if j*tau > 500 || converged || failed
        solved = 1;
    else
       j = j + 1;
    end
end

t = linspace(0,(j)*tau,length(u));
tJG = linspace(0,(j)*tau,length(base.J));
figure
plot(t,u(1:end,:))
figure
plot(tJG,base.J)
hold on
plot(tJG,other.J(:,1,1))
figure
plot(tJG,base.G)

end
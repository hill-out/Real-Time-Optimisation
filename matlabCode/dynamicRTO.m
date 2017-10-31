function [] = dynamicRTO(modelFun,plantFun,gain,u0,c0,tau)
% Takes the modelFun (either convex or other) and the plantFun and uses a
% standard multiple unit MA to work out the optimum position of the plant

% Base Case
base = struct('C', zeros(1001, 4),...
    'J', zeros(1001, 1),...
    'G', zeros(1001, 2));

[bC, bJ, bG] = plantFun(c0,u0,1000,100);

base.C = bC;
base.J = bJ;
base.G = bG;

% Other Units
delta = 50e-2;
other = struct('C', zeros(1001, 4, numel(u0)),...
    'J', zeros(1001, 1, numel(u0)),...
    'G', zeros(1001, 2, numel(u0)));

for i = 1:numel(u0)
    uNew = u0;
    uNew(i) = uNew(i) + delta;
    
    [oC, oJ, oG] = plantFun(c0,uNew,1000,100);
    
    other.C(:,:,i) = oC;
    other.J(:,:,i) = oJ;
    other.G(:,:,i) = oG;
end

% Grad
gradJ = (other.J(end,:,:) - repmat(base.J(end,:),1,1,numel(u0)))/(2*delta);
gradG = (other.G(end,:,:) - repmat(base.G(end,:),1,1,numel(u0)))/(2*delta);

% Model Steady State
m = cellFuncHandle2vec(modelFun,{[u0(1),u0(2)]},ones(size(modelFun))');

e(1) = base.J(end) - m(1);
e(2:3) = base.G(end,:)' - m(2:3);
e = e;
l = [permute(gradJ,[3,2,1]),permute(gradG,[3,2,1])];

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

    [bC, bJ, bG] = plantFun(base.C(end,:),u(j,:),1000,tau);
    
    base.C((j-1)*1000+1:j*1000+1,:) = bC;
    base.J((j-1)*1000+1:j*1000+1,:) = bJ;
    base.G((j-1)*1000+1:j*1000+1,:) = bG;
    
    % Other Units
    newC0 = other.C(end,:,:);
    clear other oC oJ oG
    
    other = struct('C', zeros(1001, 4, numel(u0)),...
        'J', zeros(1001, 1, numel(u0)),...
        'G', zeros(1001, 2, numel(u0)));
    
    for i = 1:numel(u0)
        uNew = u(j,:);
        uNew(i) = uNew(i) + delta;

        [oC, oJ, oG] = plantFun(newC0(1,:,i),uNew,1000,tau);

        other.C(:,:,i) = oC;
        other.J(:,:,i) = oJ;
        other.G(:,:,i) = oG;
    end
    
    % Grad
    gradJ = (other.J(end,:,:) - repmat(base.J(end,:),1,1,numel(u0)))/(2*delta);
    gradG = (other.G(end,:,:) - repmat(base.G(end,:),1,1,numel(u0)))/(2*delta);
    
    % Model Steady State
    m = cellFuncHandle2vec(modelFun,{u(j,:)},ones(size(modelFun))');
    mGrad = dfunc(modelFun,u(j,:));
    
    newe(1) = base.J(end) - m(1);
    newe(2:3) = base.G(end,:)' - m(2:3);
    newl = [permute(gradJ,[3,2,1]),permute(gradG,[3,2,1])]-mGrad;
    
    e(j,:) = newe.*gain + (1-gain).*e(j-1,:);
    l(:,:,j) = newl.*gain + (1-gain).*l(:,:,j-1);
    
    modifiedFun{1} = @(v)(modelFun{1}(v)+e(j,1)+(v-u(j,:))*(l(:,1,j)));
    modifiedFun{2} = @(v)(modelFun{2}(v)+e(j,2)+(v-u(j,:))*(l(:,2,j)));
    modifiedFun{3} = @(v)(modelFun{3}(v)+e(j,3)+(v-u(j,:))*(l(:,3,j)));
    
    first = 0;
    
    converged = all(all([((e(j,:) - e(j-1,:))./e(j-1,:) < 1e-3);((l(:,:,j,:) - l(:,:,j-1))./l(:,:,j-1) < 1e-3)]));
    
    if j*tau > 2000 || converged
        solved = 1;
    else
       j = j + 1;
    end
end

plot(linspace(0,(j)*tau,length(u)),u)

end
function [] = dynamicRTO(modelFun,modelGrad,plantFun,plantGrad,gain,u0)
% Takes the modelFun (either convex or other) and the plantFun and uses a
% standard MA to work out the optimum position of the plant

% Initial run
if isempty(plantGrad)
    l = dfunc(plantFun,u0,1e-2);
else
    l = cellFuncHandle2vec(plantGrad,{[u0(1),u0(2),0,0]},ones(size(plantGrad))');
end

p = cellFuncHandle2vec(plantFun,{[u0(1),u0(2),0,0]},ones(size(plantFun))');
m = cellFuncHandle2vec(modelFun,{[u0(1),u0(2)]},ones(size(modelFun))');

e = (p-m)';

modifiedFun{1} = @(v)(modelFun{1}(v)+e(1,1)+(v-u0)*(l(:,1)));
modifiedFun{2} = @(v)(modelFun{2}(v)+e(1,2)+(v-u0)*(l(:,2)));
modifiedFun{3} = @(v)(modelFun{3}(v)+e(1,3)+(v-u0)*(l(:,3)));
modifiedGrad = [];

first = 1;
solved = 0;
u = u0;
i = 2;
while solved == 0 || first == 1
    
    u(i,:) = CSTR_Opt(u0, modifiedFun{1}, {modifiedFun{2:end}});
    % calc the current plant gradient
    if isempty(plantGrad)
        pGrad = dfunc(plantFun,u(i,:));
    else
        pGrad = cellFuncHandle2vec(plantGrad,{u(i,:)},ones(size(plantGrad))');
    end
    
    % calc the current model gradient
    if isempty(modifiedGrad)
        mGrad = dfunc(modelFun,u(i,:));
    else
        mGrad = cellFuncHandle2vec(modelGrad,{u(i,:)},ones(size(modelGrad))',2);
    end
    
    % calc the current value
    p(:,i) = cellFuncHandle2vec(plantFun,{u(i,:)},ones(size(plantFun))');
    m(:,i) = cellFuncHandle2vec(modelFun,{u(i,:)},ones(size(modelFun))');
    
    e(i,:) = (p(:,i) - m(:,i))'.*gain + (1-gain).*e(i-1,:);
    l(:,:,i) = (pGrad - mGrad).*gain + (1-gain).*l(:,:,i-1);
    
    modifiedFun{1} = @(v)(modelFun{1}(v)+e(i,1)+(v-u(i,:))*(l(:,1,i)));
    modifiedFun{2} = @(v)(modelFun{2}(v)+e(i,2)+(v-u(i,:))*(l(:,2,i)));
    modifiedFun{3} = @(v)(modelFun{3}(v)+e(i,3)+(v-u(i,:))*(l(:,3,i)));
    
    i = i + 1;
    first = 0;
    
    if i == 10
        solved = 1;
    end
end


end
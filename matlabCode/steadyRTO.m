function [] = steadyRTO(modelFun,modelGrad,plantFun,plantGrad,modelPara,plantPara,gain,u0)
% Takes the modelFun (either convex or other) and the plantFun and uses a
% standard MA to work out the optimum position of the plant

% Initial run
if isempty(plantGrad)
    l = dfunc(plantFun,u0,0);
else
    l = cellFuncHandle2vec(plantGrad,ones(size(plantGrad))*u(i-1),ones(size(plantGrad)));
end

p = cellFuncHandle2vec(plantFun,ones(size(plantFun))*u(i-1),ones(size(plantFun)));
m = cellFuncHandle2vec(plantFun,ones(size(modelFun))*u(i-1),ones(size(modelFun)));


first = 1;
u = u0;
i = 2;
while solved == 0 || first == 1
    
    
    % calc the current plant gradient
    if isempty(plantGrad)
        pGrad = dfunc(plantFun,u0,0);
    else
        pGrad = cellFuncHandle2vec(plantGrad,ones(size(plantGrad))*u(i-1),ones(size(plantGrad)));
    end
    
    % calc the current model gradient
    if isempty(modelGrad)
        mGrad = dfunc(modelFun,u0,0);
    else
        mGrad = cellFuncHandle2vec(modelGrad,ones(size(modelGrad))*u(i-1),ones(size(modelGrad)));
    end
    
    % calc the current value
    p = cellFuncHandle2vec(plantFun,ones(size(plantFun))*u(i-1),ones(size(plantFun)));
    m = cellFuncHandle2vec(plantFun,ones(size(modelFun))*u(i-1),ones(size(modelFun)));
    
    
    
    i = i + 1;
    first = 0;
end


end
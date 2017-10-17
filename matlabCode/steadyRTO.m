function [] = steadyRTO(modelFun,modelGrad,plantFun,plantGrad,modelPara,plantPara,gain,u0)
% Takes the modelFun (either convex or other) and the plantFun and uses a
% standard MA to work out the optimum position of the plant

first = 1;
u = u0;
i = 2;
while solved == 0 || first == 1
    
    
    %calc the current plant gradient
    if isempty(plantGrad)
        pGrad = dfunc(plantFun,u0,0);
    else
        pGrad = plantGrad(plantPara);
    end
    
    %calc the current model gradient
    if isempty(modelGrad)
        mGrad = dfunc(modelFun,u0,0);
    else
        mGrad = modelGrad(modelPara);
    end
    
    
    i = i + 1;
    first = 0;
end


end
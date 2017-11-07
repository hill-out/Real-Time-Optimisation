%dynamictest

load('matlab.mat')

dynamicRTO(conv_funArray,plant_fun,0.5,model_opt,plant_c0,60,'MU');
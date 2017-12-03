%dynamictest

load('matlab.mat')

dynamicRTO(conv_funArray,plant_fun,0.2,model_opt,plant_c0,30,'MU');
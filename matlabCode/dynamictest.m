%dynamictest

load('dyn.mat')

dynamicRTO(conv_funArray,plant_fun,0.8,model_opt,model_cSolOpt,30,'MU');
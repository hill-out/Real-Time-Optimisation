load('matlab.mat')

steadyRTO(conv_funArray,conv_dfunArray,{plant_phi,plant_G{:}},[],0.5,model_opt);


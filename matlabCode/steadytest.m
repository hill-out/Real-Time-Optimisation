load('matlab.mat')

steadyRTO(conv_funArray,conv_dfunArray,{plant_phi,plant_G{:}},[],0.9,model_opt);


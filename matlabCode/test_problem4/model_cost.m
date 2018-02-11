function [ J, y] = model_cost( u )
%cost (model, takes in an n_u row vector and returns a scalar cost and an
%n_y row vector of model outputs
    [ J, G, Gm_only, y ] = model(u);
end

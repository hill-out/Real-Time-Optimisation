function [ G, Geq, Gm_only, y ] = model_constraint( u )
    
    [ J, G,Gm_only, y ] = model(u);
    Geq = [];
end

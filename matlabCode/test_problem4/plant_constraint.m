function [G, Geq] = plant_constraint( y )
    [ J, G ] = plant( y );
    Geq = [];
    


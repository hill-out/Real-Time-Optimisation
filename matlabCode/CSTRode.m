function dc = CSTRode(c, c_in, q_in, k, V, rO, H)
% the ode of a CSTR
%
% c - Current concentrations [nx1] (mol/L)
% cIn - Inlet concentrations [nx1] (mol/L)
% k - Reaction constants [1xm] L/(mol min)
% V - Volume (L)
% rO - Reaction order [nxm] 
%
% dc - gradient of c

%% Gradient of c
%%
% 
% $$\frac{dc}{dt} = \frac{q_{in}}{V}\times(c-c_{in}) + r$$
% 


    function r = rate(c,k,rO)
        
    end

end


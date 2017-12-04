function [u1, u2] = consArea(cons,a)
% -------------------------------------------------------------------------
% plots uA vs uB for which cons = 0
%
% cons      - struct        - structure of parameters
%
% -------------------------------------------------------------------------

u1 = 10:0.5:20;
x0 = 30;

for i = 1:numel(u1)
    
    u2(i) = fzero(@(x)(consX(cons,x)), 0);
    x0 = u2(i);
end

plot(u1,u2)

    function [c] = consX(f,x)
        c = f([u1(i),x]);
        c = c(a);
    end

end
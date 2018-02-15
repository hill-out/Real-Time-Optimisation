function [convVal] = convCalc(x,y)
% Finds the convex function values that fits this equation:
% 
% min((y(x)-convFun(convVal,x)).^2)
% where: convFun(convVal,x) = m0 + m1*x + 0.5*m2*x^2 (quadratic)
% -------------------------------------------------------------------------
% x         Function inputs [nx x m]
% y         Function outputs [ny x m]
% 
% convVal   Convex parameters [1 x 7]
% -------------------------------------------------------------------------

nx = size(x,1);
ny = size(y,1);
m = size(y,2);

convVal = fminsearch(@lsr, [1,1,1,1,1,1]);

    function [square] = lsr(v)
        f = convFun(v, x);
        square = (y - f).^2;
        while numel(square)>1
            square = sum(square);
        end
    end




end
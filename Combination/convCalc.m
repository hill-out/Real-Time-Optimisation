function [convVal] = convCalc(x,y,x_opt,y_opt,a)
% Finds the convex function values that fits this equation:
% 
% min((y(x)-convFun(convVal,x)).^2)
% where: convFun(convVal,x) = m0 + m1*x + 0.5*m2*x^2 (quadratic)
% -------------------------------------------------------------------------
% x         Function inputs [nx x m]
% y         Function outputs [ny x m]
% x_opt     Inputs at optimum [nx x 1]
% y_opt     Outputs at optimum [ny x 1]
% a         Initial guess [1 x 7] (optional)
% 
% convVal   Convex parameters [1 x 7]
% -------------------------------------------------------------------------

if nargin < 5
    a = [1, 1, 1, 0, 1];
end

options = optimset('Display','final','MaxFunEvals',5000,'MaxIter',5000);
convVal = fminsearch(@lsr, a, options);

    function [square] = lsr(v)
        f = convFun([y_opt, v], bsxfun(@minus,x,x_opt));
        square = (y - f).^2;
        while numel(square)>1
            square = sum(square);
        end
    end
end
function [dfun] = finDiff(fun, val, d)
% Find the gradient using a finite differences method
% ---------------------------------------------------
% fun       The function to differentiate
% val       The value to find the derivative at
% pm        the relative change
% 
% dfun      The derivative
% ---------------------------------------------------

f0 = fun(val);
dfun = zeros(numel(val),numel(f0));
dval = diag(d.*val);

for i = 1:numel(val)
    vali = val + dval(i,:);
    dfun(i,:) = (fun(vali)-f0)/dval(i,i);
end

end


function [df] = dfunc(f,u,e)
% funds the gradient of a function (f) at point (u) using finite the
% differences method

if nargin < 3
    e = 1e-2;
end

df = zeros(numel(u),numel(f));
for i = 1:numel(u)
    
    uUp = u;
    uUp(i) = uUp(i) + e;

    uDown = u;
    uDown(i) = uDown(i) - e;
    
    
    for j = 1:numel(f)
        df(i,j) = (f{j}(uUp)-f{j}(uDown))/(2*e);
    end
end

end
function [g] = gradFD(fun, nF, x0, dx)
% -------------------------------------------------------------------------
% calculates the gradient of f at x0 using finite differences
% 
% fun       - fun_hand      - function of x
% nF        - double        - number of outputs of fun
% x0        - double        - input to x (can be vector)
% dx        - double        - size of cahnge of x
% 
% g         - double        - gradient of f at x      
%
% -------------------------------------------------------------------------

% find f0
f0 = fun(x0);

% run for each change
f = zeros(nF, numel(x0));

for i = 1:numel(x0)
    % change each value of x
    x = x0;
    x(i) = x(i) + dx;
    
    % run f at changed x
    f(:, i) = fun(x);
end

g = (f - repmat(reshape(f0,[],1),1,nF))./dx;

end
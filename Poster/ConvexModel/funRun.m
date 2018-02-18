function [f, X] = funRun(fun, u)
% finds the values of a function
% --------------------------------------------------------
% fun       function to approximate
% u	        Range of u values to approximate at
% 
% f         Values of fun at x
% X         Mass fractions at x
% --------------------------------------------------------

% set up u

xGuess = [0.1, 0.5, 0.1, 0.1, 0.1, 0.1];
n = size(u,2);

X = zeros(6,n);
f = zeros(1,n);

% run CSTR and calc fun for all u
for i = 1:n
    X(:,i) = CSTRmodel(u(:,i),xGuess);
    
    f(1,i) = fun(u(:,i)',X(:,i)');
    xGuess = X(:,i);
end

end
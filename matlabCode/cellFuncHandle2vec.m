function [funcOut] = cellFuncHandle2vec(cellFunc,para,paraNum,o)
% takes a cell of function handles that take para{paraNum} and outputs
% their outputs as a vector (used in fmincon)

if nargin < 4
    o = 1;
end

funcOut = zeros(length(cellFunc),o);

for i = 1:numel(cellFunc)
    funcOut(i,:) = cellFunc{i}(para{paraNum(i,:)});
end

end
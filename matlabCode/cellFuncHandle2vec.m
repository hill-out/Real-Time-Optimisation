function [funcOut] = cellFuncHandle2vec(cellFunc,para,paraNum)
% takes a cell of function handles that take para{paraNum} and outputs
% their outputs as a vector (used in fmincon)

funcOut = zeros(size(cellFunc));

for i = 1:numel(cellFunc)
    funcOut(i) = cellFunc{i}(para{paraNum(i,:)});
end

end
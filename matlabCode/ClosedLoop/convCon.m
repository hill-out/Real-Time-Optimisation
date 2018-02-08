function [clt,ceq] = convCon(a,b,c,in,inOpt)
% calculates the values of the constratints
% --------------------------------------------------------
% a         Zeroth order term (1xng)
% b         First order term (nxng)
% c         Second order term (n^2xng)
% in        Inputs to the convex approximation (nx1)
% inOpt     Model optimum inputs (nx1)
% 
% clt       Less than constraint (ngx1)
% ceq       Equal to constraint (ngx1)
% --------------------------------------------------------

n = numel(in);
ng = numel(a);

a = reshape(a,1,1,ng);
b = permute(reshape(b,1,n,ng),[3,2,1]);                  % make a row
c = permute(reshape(c,n,n,ng),[3,2,1]); % make a square
in = reshape(in,n,1);                   % make a column
inOpt = reshape(inOpt,n,1);             % make a colummn

for i = 1:ng
    clt(i) = a(:,:,i) + b(:,i)'*(in-inOpt) + 0.5*(in-inOpt)'*c(:,:,i)*(in-inOpt);
end
clt = reshape(clt,ng,1);
ceq = [];

end
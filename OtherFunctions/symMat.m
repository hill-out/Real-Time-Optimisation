function mat = symMat(n)
% makes a nxn symmetric matrix
% some how works?

m = sum(1:n);

mat = zeros(n);

for i = 1:m
    j = i;
    l = 0;
    while j>0
        k = j;
        j = j - n + l;
        l = l + 1;
    end
    j = j + n;
    mat(j,l) = i;
    mat(l,j) = i;
end

end
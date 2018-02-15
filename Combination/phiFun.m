function phi = phiFun(u,X)
% calculates phi from u and x
F = u(:,1)+u(:,2);
phi = -(1143.38.*X(:,4).*F + 25.92.*X(:,5).*F - 76.23.*u(:,1) - 114.34.*u(:,2));
end
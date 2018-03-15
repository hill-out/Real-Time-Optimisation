function [y] = trendNoise(n,pow)
% creates trended noise
% ---------------------
% n         number of points
% pow       power of noise
%
% y         noise
% ---------------------

y = randn/pow;

for i = 2:n
    y(i) = y(i-1)*5/6 + randn/pow;
end

end
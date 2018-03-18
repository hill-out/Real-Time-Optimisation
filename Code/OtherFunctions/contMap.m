function [fig] = allPlots(X_As, F_Bs, fig)
% run the plant of choice until steady for a map of points to get contours

if nargin < 2
    F_Bs = linspace(3,12,41);
    if nargin < 1
        X_As = linspace(0.08,0.12,21);
    end
end

[a,b] = meshgrid(X_As, F_Bs);

rRange = [a(:),b(:)];

Kp = -1000;
T0 = 120;
xGuess = [0.08, 0.37, 0.1, 0.25, 0.1, 0.1];

for i = 1:size(a,1)
    for j = 1:size(a,2)
        r = [a(i,j),b(i,j)];
        up = plantController2(r,xGuess,Kp,T0)';
        [~,c] = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 10000],[up, xGuess]);
        up0 = c(end,1:3);
        Xp0 = c(end,4:end);
        
        phip(i,j) = phiFun(up0,Xp0);
        g1p(i,j) = g1Fun(up0,Xp0);
        g2p(i,j) = g2Fun(up0,Xp0);
        
        xGuess = Xp0;
    end
end

try figure(fig.rr) % r1 vs. r2
    hold on
catch
    fig.rr = figure('name','r1 vs. r2');
end
contour(X_As,F_Bs,phip,20,'ShowText','on')
contour(X_As,F_Bs,g1p,[0,0])

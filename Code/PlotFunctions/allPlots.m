function [fig] = allPlots(rpi, rp0, base, algo,tau, K, fig)
% Plots rr, phit, ut and gt
% -----------------------------------
% rpi       Iterative plant setpoints
% rp0       First setpoint from model
% base      Struct of plant time data
% algo      Algorithm for closed loop
% tau, K    Freq and gain of solution
% fig       Struct of figures to plot
% -----------------------------------
set(0,'DefaultFigureWindowStyle','docked')
dataName = sprintf('%s / \x03c4 = %4.0f / K = %5.2f',algo,tau,K);
try figure(fig.rr) % r1 vs. r2
    hold on
catch
    fig.rr = figure('name','r1 vs. r2');
    
end
plot([rp0(1), rpi(:,1)'],[rp0(2), rpi(:,2)'],'DisplayName',dataName)
legend('-DynamicLegend');
legend show

try figure(fig.phit) % phi vs t
    hold on
catch
    fig.phit = figure('name','phi vs t');
end
plot(base.t,-base.phip,'DisplayName',dataName)
legend('-DynamicLegend');
legend show

try figure(fig.ut) % u vs t
catch
    fig.ut = figure('name','u vs t');
end
subplot(1,3,1)
hold on
plot(base.t,base.u(:,1),'DisplayName',dataName)
legend('-DynamicLegend');
legend show
subplot(1,3,2)
hold on
plot(base.t,base.u(:,2),'DisplayName',dataName)
legend('-DynamicLegend');
legend show
subplot(1,3,3)
hold on
plot(base.t,base.u(:,3),'DisplayName',dataName)
legend('-DynamicLegend');
legend show

try figure(fig.gt) % g vs t
    hold on
catch
    fig.gt = figure('name','g vs t');
end
plot(base.t,base.g1p,'DisplayName',dataName)
legend('-DynamicLegend');
legend show

end
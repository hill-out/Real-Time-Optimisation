clear
close all
allTau = ones(10,1)*20;
allK = 0.45:0.05:0.55;

for i = 1:3
    tau = allTau(i);
    K = allK(i);
    T_closedRR
    clearvars -except fig allTau allK tFinal
end

optionu = optimoptions('fmincon','Display','off');
xGuess = [0.09, 0.36, 0.1, 0.25, 0.1, 0.1];

T0 = 120;
Kp = -1000;
rp_opt   = [0.10605,6.05325];
up_opt   = plantController2(rp_opt,xGuess,Kp,T0)';
[~,a]    = ode15s(@(t,y)closedPlantODE(t,y,Kp), [0 100000],[up_opt, xGuess]);
up_opt   = a(end,1:3);
Xp_opt   = a(end,4:end);
phip_opt = phiFun(up_opt,Xp_opt);
g1p_opt  = g1Fun(up_opt,Xp_opt);
g2p_opt  = g2Fun(up_opt,Xp_opt);

figure(fig.rr)
plot(rp_opt(1),rp_opt(2),'rx','LineWidth',4,'MarkerSize',15)

figure(fig.phit)
plot([0 tFinal],-[phip_opt,phip_opt],'r--')

figure(fig.ut)
subplot(1,3,1)
plot([0 tFinal],[up_opt(1),up_opt(1)],'r--')
subplot(1,3,2)
plot([0 tFinal],[up_opt(2),up_opt(2)],'r--')
subplot(1,3,3)
plot([0 tFinal],[up_opt(3),up_opt(3)],'r--')

figure(fig.gt)
plot([0 tFinal],[g1p_opt,g1p_opt],'r--')
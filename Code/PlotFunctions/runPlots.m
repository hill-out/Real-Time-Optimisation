T_closedUU
T_closedUR
T_closedRR

figure(fig.rr)
plot(rp_opt(1),rp_opt(2),'rx')

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
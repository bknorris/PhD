%calculate bed shear stress from near-bed Reynolds stress. HTA1 and VTA.

clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
rho = 1029; %kg/m^3
time = RS.day1.vpro2.time;
tau1 = abs((RS.day1.vpro1.uw(:,5)./100)*rho);
tau2 = abs((RS.day1.vpro2.uw(:,5)./100)*rho);
tau3 = abs((RS.day1.vpro3.uw(:,5)./100)*rho);
tau4 = abs((RS.day4.vpro1.uw(:,5)./100)*rho);

plot(time,tau1),hold on
plot(time,tau2,'r')
plot(time,tau3,'g')
plot(RS.day4.vpro1.time,tau4,'c')

disp(['VP1 bed shear stress (average): ' num2str(nanmean(tau1)) ' N/m^2'])
disp(['VP2 bed shear stress (average): ' num2str(nanmean(tau2)) ' N/m^2'])
disp(['VP3 bed shear stress (average): ' num2str(nanmean(tau3)) ' N/m^2'])
disp(['VTA VP1 bed shear stress (average): ' num2str(nanmean(tau4)) ' N/m^2'])
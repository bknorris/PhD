%Calculate wave boundary layer thickness

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\AD5117_12March2015wvs')
time = wave.time2;
Aw = wave.Uwrms*sqrt(2)./wave.omegar; %wave orbital excursion
REw = wave.Uwrms*sqrt(2).*Aw./1E-6; %wave reynolds number
%Fredsoe and Deigaard say that boundary layer flow is laminar if REw <
%5x10^4
id = find(REw < 5E4);
Re = NaN(length(Aw),1);
Re(id) = REw(id);
Odel = (1E-6./wave.omegar).^(1/2); %order of magnitude laminar wave boundary layer thickness (Luhar et al 2010).

% figure
% plot(time,REw),hold on
% plot(time,Re,'r')
% datetick('x','dd HH:MM:SS')
% ylabel('Wave Reynolds #')

%basically, no times during the experiment were below the threshold of 5E4
%for wave RE numbers.
% Kn = 0.1; %order of magnitude estimate for Nikuradse's roughless length (see: Augustijn et al. 2008)
%Julia thinks Kn should be the d50, i.e. those other studies are looking at
%more regional lengthscales (corresponding to the vegetation, 0.1m is quite
%large)
Kn = 100E-6; %from Fricke et al. 2017 (d50 of sediment from SW during March)
Wdel = 0.09*Kn*(Aw./Kn).^0.82; %eqn. from Fredsoe & Deigaard, 1992
C4 = [0.072 0.896 0.111];
C5 = [-0.25 -0.469 -0.246];

Wdel2 = zeros(3,length(Aw));
for i = 1:3
    Wdel2(i,:) = Aw.*C4(i).*(Aw./Kn).^(C5(i));
end

figure
c = lines(4);
p = zeros(4,1);
p(1) = plot(time,Wdel.*1000,'color',c(1,:));hold on
for i = 1:3
    p(i+1) = plot(time,Wdel2(i,:).*1000,'color',c(i+1,:));
end
text = {'F & D';'Jonsson';'Sleath';'JSF'};
legend(p,text)
datetickzoom('x','dd HH:MM:SS')
ylabel('\delta_w (mm)')
title('Wave Boundary Layer for F2F3')

%compare velocity profile to 3 VPs
% load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\F2F3_2Vels.mat')
% start = datenum(2015,03,11,15,10,00);
% stop = datenum(2015,03,11,15,15,00);
% fn = fieldnames(dat);
% vph = [61,61,61];
% figure
% pp = zeros(3,1);
% c = lines(3);
% for i = 1:3
%     time2 = dat.(fn{i}).time;
%     id = find(time2 >= start & time2 <= stop);
%     z = dat.(fn{i}).z1(id,:);
%     rb = vph(i)-40-linspace(0,30,35);
%     xprof = nanmean(z);
%     
%     %%%
%     pp(i) = plot(xprof,rb,'color',c(i,:));
%     hold on
% end
    
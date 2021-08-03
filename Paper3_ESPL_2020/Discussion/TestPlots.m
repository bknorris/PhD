clear

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat')
area = {'mud';'fringe';'forest'};
figure(1), hold on
for k = 1:length(area)
s(k) = subplot(1,3,k);
%load variables
taub = [dat.(area{k}).wave.taub; dat.(area{k}).ig.taub];
Hs = [dat.(area{k}).wave.sigh; dat.(area{k}).ig.sigh];
h = [dat.(area{k}).wave.depth; dat.(area{k}).ig.depth];
eps = [dat.(area{k}).wave.eps; dat.(area{k}).ig.eps];
netbd = abs(detrend([dat.(area{k}).wave.deltbd; dat.(area{k}).ig.deltbd]));
%filter netbd around 0
zid = find(netbd<0.002 | netbd > 0.08); %find changes less than 2 mm or greater than 80 mm
taub(zid) = []; %remove pts
Hs(zid) = [];
h(zid) = [];
eps(zid) = [];
netbd(zid) = [];

%calculate quartiles of net BLE change
q1 = median(netbd(find(netbd<median(netbd))));
q3 = median(netbd(find(netbd>median(netbd))));
IQR = q3-q1;

% determine extreme Q3 outliers (e.g., x > Q1 + 3*IQR)
if 0.005 > q3+3*IQR %greater than 5 mm change
    iy = find(netbd>0.005);
else
    iy = find(netbd>q3+3*IQR);
end
netbed = netbd(iy);
bedshear = taub(iy);
gamma = Hs(iy)./h(iy);
plot(netbed,gamma,'o')
end
% plot3(x(id),y(id),z(id),'o')
% grid on
set(s,'ylim',[0 0.4],'xlim',[0 0.1])
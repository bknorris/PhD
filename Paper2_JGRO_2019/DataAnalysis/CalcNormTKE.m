%%Load the Averaged TKE files. Find max of the absolute value of TKE and 
%extract the height of that measurement in the profile. 
%Save these data.
clear
close all
figdir = 'c:\Users\Bnorr\Documents\GradSchool\Papers\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';
ddir = 'd:\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\';
wdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\WPF\';
dn = {'day1';'day2';'day3';'day4'};
files = {'HTA_1TKE_bdadj';'HTA_2TKE';'HTA_4TKE';'VTA_2TKE_bdadj'};
wfiles = {'WpfV5108_HTA_1';'Wpf5116_HTA_1';...
    'WpfHR3_HTA_2';'Wpf5116_HTA_2';...
    'WpfHR3_HTA_4';'Wpf5116_HTA_4';...
    'WpfV5109_VTA';'Wpf5116_VTA'};
wf = [1 2;3 4;5 6;7 8];
delx = [65 56 56 31]; %m
theta = [-54 0 0 0]; %deg
rho = 1011.4; %kg/m^3
for i = 1:4
    load([ddir files{i} '.mat'])
    w1 = load([wdir wfiles{wf(i,1)}]);
    w2 = load([wdir wfiles{wf(i,2)}]);
    fn = fieldnames(Stat);
    for ii = 1:3
        %Average data into 10-minute chunks
        n = length(Stat.(fn{ii}).time);
        fs = 0.1; %Hz, epsilon calculated every 10 sec
        avt = 10*60*fs;
        idx = [1 avt:avt:n];
        emax = zeros(length(idx)-1,1);
        emean = zeros(length(idx)-1,1);
        eid = zeros(length(idx)-1,1);
        time = zeros(length(idx)-1,1);
        h = w2.WPF.h;
        wpf = abs((w1.WPF.F-w2.WPF.F)./(delx(i).*cosd(theta(i))));
        if i == 1 || (i == 4 && ii == 1)
            bins = 1:5;
        else
            bins = 13:17;
        end
        for j = 1:length(idx)-1
            ids = idx(j):idx(j+1);
            E = (Stat.(fn{ii}).z1.E(bins,ids)+...
                    Stat.(fn{ii}).z2.E(bins,ids))./2;
            enorm = (nanmean(E,2).*mean(h(ids))*rho)./mean(wpf(ids));
            emean(j) = mean(nanmean(E,2));
            [~,ex] = max(enorm);
            emax(j) = enorm(ex);eid(j) = bins(ex);
            time(j) = Stat.(fn{ii}).time(idx(j));
        end
        mydata.(dn{i}).(fn{ii}).time = time;
        mydata.(dn{i}).(fn{ii}).emean = emean;
        mydata.(dn{i}).(fn{ii}).emax = emax;
        mydata.(dn{i}).(fn{ii}).eid = eid;
    end
end
savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\TKE\';
save([savedatdir 'NormalizedTKE'],'mydata','-v7.3')
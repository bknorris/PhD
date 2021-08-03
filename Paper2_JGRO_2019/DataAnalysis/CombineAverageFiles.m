%Load the Averaged Vectrino files, the Averaged ADV files, and the Averaged
%Duet. Compute variables of interest
clear
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\';
files = {'HTA_1Vave.mat';'HTA_2Vave.mat';'HTA_4Vave.mat';...
    'VTA_2vp1Vave.mat';'VTA_2vp2Vave.mat';'VTA_2vp3Vave.mat'};
dates = {'day1';'day2';'day3';'day4';'day4';'day4'};
VecAve = struct();
for i = 1:length(files)
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    if length(fn) > 1
        g = 1:3;
    else
        g = 1;
    end
    for ii = g
        VecAve.(dates{i}).(fn{ii}).time = Avgs.(fn{ii}).time;
        VecAve.(dates{i}).(fn{ii}).Uc = Avgs.(fn{ii}).Uc;
        VecAve.(dates{i}).(fn{ii}).Uw = Avgs.(fn{ii}).Uw;
    end
end
load('d:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\HTA_ADVave.mat')
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\Duet_140315VTA2wvs.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 07:01:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 13:19:00'};
mydata = struct();
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
phi = [0.07 0.03 1.6E-4 0.03]; %phi, from table 2 (z/hc)
 for i = 1:4
    hc_i = hc(i);
    if i < 4
        fn = fieldnames(VecAve);
        t1 = Vave.(fn{i}).time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        Uwinf = Vave.(fn{i}).Uw(id(1):id(2));
        Ucinf = Vave.(fn{i}).Uc(id(1):id(2));
        Hratio = Vave.(fn{i}).h(id(1):id(2))./hc_i;
        dfn = fieldnames(VecAve.(fn{i}));
        for ii = 1:3
            t2 = VecAve.(fn{i}).(dfn{ii}).time;
            Uwc = nanmean(VecAve.(fn{i}).(dfn{ii}).Uw);
            Ucc = nanmean(VecAve.(fn{i}).(dfn{ii}).Uc);
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            Uwc = Uwc(id(1):id(2));
            Ucc = Ucc(id(1):id(2));
            t2 = t2(id(1):id(2));
            Uwratio = Uwc./Uwinf;
            Ucratio = Ucc./Ucinf;
            %quality control, remove very large values
            mu = mean(Uwratio);
            stdv = std(Uwratio);
            nid = find(Uwratio < mu+3*stdv & Uwratio > mu-3*stdv);
            Uwratio = Uwratio(nid);
            mu = mean(Ucratio);
            stdv = std(Ucratio);
            nid = find(Ucratio < mu+3*stdv & Ucratio > mu-3*stdv);
            Ucratio = Ucratio(nid);
            Hr = Hratio(nid);
            time = t2(nid);
            %save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = time;
            mydata.(fn{i}).(dfn{ii}).Uwratio = Uwratio;
            mydata.(fn{i}).(dfn{ii}).Ucratio = Ucratio;
            mydata.(fn{i}).(dfn{ii}).Hratio = Hr;
            mydata.(fn{i}).(dfn{ii}).phi = phi(i);
        end
    else
        t1 = wave.time2;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        Hratio = wave.h(id(1):id(2))./hc_i;
        %%%
        t2 = VecAve.day4.vpro3.time;
        [~,id(1)] = min(abs(t2-start));
        [~,id(2)] = min(abs(t2-stop));
        Uwinf = nanmean(VecAve.day4.vpro3.Uw(:,id(1):id(2)));
        Ucinf = nanmean(VecAve.day4.vpro3.Uc(:,id(1):id(2)));
        t2 = VecAve.day4.vpro1.time;
        [~,id(1)] = min(abs(t2-start));
        [~,id(2)] = min(abs(t2-stop));
        Uwc = nanmean(VecAve.day4.vpro1.Uw(:,id(1):id(2)));
        Ucc = nanmean(VecAve.day4.vpro1.Uc(:,id(1):id(2)));
        Uwratio = Uwc(1:length(Uwinf))./Uwinf;
        Ucratio = Ucc(1:length(Ucinf))./Ucinf;
        Hratio = Hratio(1:length(Ucinf));
        %quality control, remove very large values
        mu = mean(Uwratio);
        stdv = std(Uwratio);
        nid = find(Uwratio < mu+2*stdv & Uwratio > mu-2*stdv);
        Uwratio = Uwratio(nid);
        mu = mean(Ucratio);
        stdv = std(Ucratio);
        nid = find(Ucratio < mu+2*stdv & Ucratio > mu-2*stdv);
        Ucratio = Ucratio(nid);
        Hr = Hratio(nid);
        time = t2(nid);
        %save variables to structure
        mydata.(fn{i}).time = time;
        mydata.(fn{i}).Uwratio = Uwratio;
        mydata.(fn{i}).Ucratio = Ucratio;
        mydata.(fn{i}).Hratio = Hr;
        mydata.(fn{i}).phi = phi(i);
    end
end
savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\';
save([savedatdir 'CurrentHPhi_data'],'mydata','-v7.3')
%Calculate sediment parameters: tau_cr from Wu et al. 2017, and U_cr from Van Rijn (1993)
clear
rdir = 'g:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'SedimentData.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%n','delimiter',',');
rexp = rfile{1}; %experiment name
rdate = rfile{2}; %date of experiment
rfld = rfile{3}; %area of forest
rd50 = rfile{4}; %sediment d50 (m)
rd90 = rfile{5}; %sediment d90 (m)
rsand = rfile{6}; %percent sand in sample (%)
rmud = rfile{7}; %percent mud in sample (%)
tcr = zeros(length(rexp),1);
for i = 1:length(rexp)
    if rmud(i)/100 < 0.15
        type  = 'muddy';
    else
        type = 'sandy';
    end
    tcr(i) = critstresswu(rsand(i)/100,rmud(i)/100,rd50(i),type);
end
tcr(tcr>2.1) = NaN; %maximum on graph (Smith et al. 2015)
id = strcmp(rfld,'mud');
fprintf('Average mudflat critical shear stress: %0.2f Pa\n',nanmean(tcr(id)))
id = strcmp(rfld,'fringe');
fprintf('Average fringe critical shear stress: %0.2f Pa\n',nanmean(tcr(id)))
id = strcmp(rfld,'forest');
fprintf('Average forest critical shear stress: %0.2f Pa\n',nanmean(tcr(id)))
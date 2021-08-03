%Create dat files to be read into the other VP processing scripts. Unlike
%AutomatedQC_Vpros.m, this script does not quality control, it only
%combined files from the three VPs into two files: dat (for data) and sen
%(for ancillary sensor data- SNR, Amp, Cor). Note: QC is now run at the
%beginning of data processing ("runVecPro_v2.m")
%
% This script was written by Benjamin K Norris, 2017
% University of Waikato, New Zealand
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
year = '2015';
var = {'vpro1';'vpro2';'vpro3'}; %substructure variable names

%find files
fdir = 'd:\Projects\Mekong_W2015\Data\Vectrino\';
tier1 = dir(fdir);
isub = [tier1(:).isdir];
tier1 = {tier1(isub).name}';
tier1 = tier1(3:end);
toprocess = [1];
for i = toprocess
    subfdir = [fdir tier1{i} '\'];
    if strcmp(year,'2014')
        files = regexprep(tier1{i},'[^a-zA-Z\d\s:]','_');
    else
        files = tier1{i};
    end
    tier2 = dir([subfdir files '*']);tier2 = {tier2.name};
    for ii = 1:3
        disp(['Loading ' subfdir tier2{ii}])
        load([subfdir tier2{ii}])
        
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        dat.(var{ii}).time = VPRO.Data.Time;
        dat.(var{ii}).x = VPRO.Data.Profiles_VelX;
        dat.(var{ii}).y = VPRO.Data.Profiles_VelY;
        dat.(var{ii}).z1 = VPRO.Data.Profiles_VelZ1;
        dat.(var{ii}).z2 = VPRO.Data.Profiles_VelZ2;
        dat.(var{ii}).rb = VPRO.Data.Profiles_Range;
        dat.(var{ii}).sr = VPRO.Config.sampleRate;
        dat.(var{ii}).bdtime = VPRO.Data.BottomCheck_HostTimeMatlab;
        dat.(var{ii}).bdist = VPRO.Data.BottomCheck_BottomDistance;
        %%%
        sen.(var{ii}).time = VPRO.Data.Time;
        sen.(var{ii}).Amp1 = VPRO.Data.Profiles_AmpBeam1;
        sen.(var{ii}).Amp2 = VPRO.Data.Profiles_AmpBeam2;
        sen.(var{ii}).Amp3 = VPRO.Data.Profiles_AmpBeam3;
        sen.(var{ii}).Amp4 = VPRO.Data.Profiles_AmpBeam4;
        sen.(var{ii}).Cor1 = VPRO.Data.Profiles_CorBeam1;
        sen.(var{ii}).Cor2 = VPRO.Data.Profiles_CorBeam2;
        sen.(var{ii}).Cor3 = VPRO.Data.Profiles_CorBeam3;
        sen.(var{ii}).Cor4 = VPRO.Data.Profiles_CorBeam4;
        clear VPRO
        
    end
    %save out to data file
    if strcmp(year,'2014')
        spath = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper3\VPs\';
        dname = [files '_Vels'];sname = [files '_Sen'];
    elseif strcmp(year,'2015')
        spath = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
        filen = regexp(tier2{ii},'(.*?)(?=\_)','tokens');
        dname = [char(filen{1}) '_Vels'];sname = [char(filen{1}) '_Sen'];
    end
    disp(['Saving ' dname ', ' sname  ', velocity & sensor files'])
    save([spath dname],'-struct','dat','-v7.3')
    save([spath sname],'-struct','sen','-v7.3')
    
    clear dat sen
end
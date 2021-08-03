%this script loads and QCs all Vectrino data currently un-QC'd
clear
maindir = 'D:\Projects\Mekong_F2014\Data\Vectrino\';
savedatdir = 'D:\Projects\Mekong_F2014\DataAnalysis\Paper2\';
dateofexperiment = '2014'; %'2015'
folders = dir(maindir);folders = {folders.name};
toprocess = 9; %folders to process
var = {'vpro1';'vpro2';'vpro3'};                                            %substructure variable names
tread = fopen([savedatdir 'cropbottom.txt']);                               %Load VPRO processing commands
crop = textscan(tread,'%n','delimiter',',');
crop = crop{1};
if strcmp(dateofexperiment,'2014')
    cr = reshape(crop,[3,8]);
elseif strcmp(dateofexperiment,'2015')
    cr = reshape(crop,[3,10]);
end

for i = 1:length(toprocess)                                                 %loop through directories
    fdir = [maindir folders{toprocess(i)} '\'];
    fname = dir([fdir 'VP*_*.mat']);
    fname = {fname.name};
    
    for ii = 1:length(fname)                                                %loop through instrument files
        disp(['Loading ' fdir fname{ii}])
        load([fdir fname{ii}])
        disp('Running QC processing')
        tic
        
        %VecPro basic QC
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,1,1,0);        %run QC on data
        downwardfacing = cr(ii,i);                                          %optionally run bottom tracking for near-bed instruments
        if downwardfacing == 1
            VPRO.Data = fixbbvpvels(VPRO.Data,2);
        end
        disp(['Basic Quality Controls completed in: ' num2str(toc/60) ' minutes'])
        %%%% extract time of average
        if strcmp(dateofexperiment,'2014')
            tc = datenum(0,0,0,6,0,0);                                      %2014 used NZT as standard time
            dat.(var{ii}).time = VPRO.Data.Profiles_HostTimeMatlab-tc;
        elseif strcmp(dateofexperiment,'2015')
            tc = datenum(0,0,0,7,0,0);                                      %2015 used ICT as standard time
            dat.(var{ii}).time = VPRO.Data.Profiles_HostTimeMatlab+tc;

        end
        
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        dat.(var{ii}).time = VPRO.Data.Profiles_HostTimeMatlab;
        dat.(var{ii}).x = VPRO.Data.Profiles_VelX;
        dat.(var{ii}).y = VPRO.Data.Profiles_VelY;
        dat.(var{ii}).z1 = VPRO.Data.Profiles_VelZ1;
        dat.(var{ii}).z2 = VPRO.Data.Profiles_VelZ2;
        dat.(var{ii}).rb = VPRO.Data.Profiles_Range;
        dat.(var{ii}).sr = VPRO.Config.sampleRate;
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'xb');
        dat.(var{ii}).beam1 = VPRO.Data.Profiles_VelBeam1;
        dat.(var{ii}).beam2 = VPRO.Data.Profiles_VelBeam2;
        dat.(var{ii}).beam3 = VPRO.Data.Profiles_VelBeam3;
        dat.(var{ii}).beam4 = VPRO.Data.Profiles_VelBeam4;
        dat.(var{ii}).SNR1 = VPRO.Data.Profiles_SNRBeam1;
        dat.(var{ii}).SNR2 = VPRO.Data.Profiles_SNRBeam2;
        dat.(var{ii}).SNR3 = VPRO.Data.Profiles_SNRBeam3;
        dat.(var{ii}).SNR4 = VPRO.Data.Profiles_SNRBeam4;
        clear VPRO
    end
    %save out to data file
    disp(['Saving ' folders{toprocess(i)} ' velocity file'])
    filen = [folders{toprocess(i)} 'Vels'];
    save([savedatdir filen],'dat','-v7.3')
    
    clear dat
end

%just for now, delete this section once the data is all processed!

maindir = 'D:\Projects\Mekong_W2015\Data\Vectrino\';
savedatdir = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper2\';
dateofexperiment = '2015'; %'2015'
folders = dir(maindir);folders = {folders.name};
toprocess = 5; %folders to process
var = {'vpro1';'vpro2';'vpro3'};                                            %substructure variable names
tread = fopen([savedatdir 'cropbottom.txt']);                               %Load VPRO processing commands
crop = textscan(tread,'%n','delimiter',',');
crop = crop{1};
if strcmp(dateofexperiment,'2014')
    cr = reshape(crop,[3,8]);
elseif strcmp(dateofexperiment,'2015')
    cr = reshape(crop,[3,10]);
end

for i = 1:length(toprocess)                                                 %loop through directories
    fdir = [maindir folders{toprocess(i)} '\'];
    fname = dir([fdir 'VP*_*.mat']);
    fname = {fname.name};
    
    for ii = 1:length(fname)                                                %loop through instrument files
        disp(['Loading ' fdir fname{ii}])
        load([fdir fname{ii}])
        disp('Running QC processing')
        tic        
        
        %some VP files are too large. Crop
        if toprocess(i) == 6
            start = datenum(2015,03,13,18,50,00);
            stop = datenum(2015,03,13,23,55,00);
            VPRO.Data = cropvp(VPRO.Data,start,stop);
        elseif toprocess(i) == 7
            start = datenum(2015,03,14,07,00,00);
            stop = datenum(2015,03,14,15,00,00);
            VPRO.Data = cropvp(VPRO.Data,start,stop);
        elseif toprocess(i) == 12
            start = datenum(2015,03,09,04,30,00);
            stop = datenum(2015,03,09,07,30,00);
            VPRO.Data = cropvp(VPRO.Data,start,stop);
        else
            %%%% extract time of average
            if strcmp(dateofexperiment,'2014')
                tc = datenum(0,0,0,6,0,0);                                      %2014 used NZT as standard time
                VPRO.Data.Time = VPRO.Data.Profiles_HostTimeMatlab-tc;
            elseif strcmp(dateofexperiment,'2015')
                tc = datenum(0,0,0,7,0,0);                                      %2015 used ICT as standard time
                VPRO.Data.Time = VPRO.Data.Profiles_HostTimeMatlab+tc;
            end
        end
        %VecPro basic QC
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,1,1,0);        %run QC on data
        downwardfacing = cr(ii,i);                                          %optionally run bottom tracking for near-bed instruments
%         if downwardfacing == 1
%             VPRO.Data = fixbbvpvels(VPRO.Data,2);
%         end
        disp(['Basic Quality Controls completed in: ' num2str(toc/60) ' minutes'])

        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        dat.(var{ii}).time = VPRO.Data.Profiles_HostTimeMatlab;
        dat.(var{ii}).x = VPRO.Data.Profiles_VelX;
        dat.(var{ii}).y = VPRO.Data.Profiles_VelY;
        dat.(var{ii}).z1 = VPRO.Data.Profiles_VelZ1;
        dat.(var{ii}).z2 = VPRO.Data.Profiles_VelZ2;
        dat.(var{ii}).rb = VPRO.Data.Profiles_Range;
        dat.(var{ii}).sr = VPRO.Config.sampleRate;
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'xb');
        dat.(var{ii}).beam1 = VPRO.Data.Profiles_VelBeam1;
        dat.(var{ii}).beam2 = VPRO.Data.Profiles_VelBeam2;
        dat.(var{ii}).beam3 = VPRO.Data.Profiles_VelBeam3;
        dat.(var{ii}).beam4 = VPRO.Data.Profiles_VelBeam4;
        dat.(var{ii}).SNR1 = VPRO.Data.Profiles_SNRBeam1;
        dat.(var{ii}).SNR2 = VPRO.Data.Profiles_SNRBeam2;
        dat.(var{ii}).SNR3 = VPRO.Data.Profiles_SNRBeam3;
        dat.(var{ii}).SNR4 = VPRO.Data.Profiles_SNRBeam4;
        clear VPRO
        if toprocess(i) == 6 || toprocess(i) == 7 || toprocess(i) == 12
                %The VTA experiments are too large to save as one big file.
                %Save each vectrino individually. 
                disp(['Saving ' folders{toprocess(i)} ' velocity file'])
                filen = [var{ii} folders{toprocess(i)} 'Vels'];
                save([savedatdir filen],'dat','-v7.3')
                clear dat
        end
    end
    if any(toprocess(i) == setxor(toprocess,[6 7 12]))
        %save out to data file
        disp(['Saving ' folders{toprocess(i)} ' velocity file'])
        filen = [folders{toprocess(i)} 'Vels'];
        save([savedatdir filen],'dat','-v7.3')
        clear dat
    end
end
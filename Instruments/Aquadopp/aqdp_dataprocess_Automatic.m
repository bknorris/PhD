% <<< Change these values >>>
metdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Weather\';
dirname = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\Raw\'; %file local directory
fname = 'HR3_15March2015'; %file name of the aquadopp files
settings.heading = 18; %compass heading of instrument measured in the field
settings.ornt = 'DOWN'; %instrument orientation
settings.HAB = 700; %instrument Height Above Bottom (in mm)
settings.lat  = 9.565444411; %instrument latitude in dms (dd,mm,ss)
settings.lon = 106.2925574; %instrument longitude in dms (dd,mm,ss)
settings.magdec = 0.27; %magnetic declination to be applied
settings.depdate = datenum(2015,3,13,15,30,00); %date deployed (incl. time)
settings.recdate = datenum(2015,3,14,12,00,00); %date recovered (incl. time)

%% Load the data
cd(dirname)
files = dir([fname '*']);
filelist = {files.name};varname = filelist{1};
newvarname = regexprep(varname,'.a1','');

DATA = aqdp_read(newvarname); %read the data
fn = fieldnames(settings);
for i = 1:length(fn)
    DATA.metadata.(fn{i}) = settings.(fn{i});
end
%these are used in the code:
ornt = settings.ornt;heading = settings.heading;magdec = settings.magdec;
start = settings.depdate;stop = settings.recdate;

%crop fields to start and stop times
disp('Clipping datafile to the length of deployment')
int = find(DATA.datenum >= start & DATA.datenum <= stop);
fn = fieldnames(DATA);
if isfield(DATA,'cor1') %indices of fieldnames differ between LR and HR
    for i = [3:22 25:27]; %clip all fields to start/stop times
        DATA.(fn{i}) = DATA.(fn{i})(int,:);
    end
else
    for i = [3:17 20:22]; %clip all fields to start/stop times
        DATA.(fn{i}) = DATA.(fn{i})(int,:);
    end
end

%run the pressure adjustment
if 1 %turn off (set to 0) if you don't want to run this step
    %met is a structure with fields:
    %met.datetime
    %met.pressure
load([metdir 'Mekong2015_metstation.mat'])
t = DATA.temperature;
s = str2double(DATA.metadata.salinity);
ss = DATA.sound_spd;
p = DATA.pressure;

met.pressure = met.pressure.*0.0001; %convert Pa to dBar
ind = find(met.datetime >= start & met.datetime <= stop);
mettime = met.datetime(ind);
metpres = met.pressure(ind);
DATA.pressure = correctatmp(mettime,metpres,DATA.datenum,p,t);
end

%speed of sound adjustment for velocimeters
[DATA.vel_b1,DATA.vel_b2,DATA.vel_b3] = adjsound(t,s,ss,DATA.vel_b1,DATA.vel_b2,DATA.vel_b3);
clear t s ss p

%go to the unwrapping step. Unwrap velocities burst-by-burst,
%bin-by-bin. First check for bad correlations:
nlin = 1E2;maxbr = 1E3;maxg = 1E4; %cmgbridge settings
spb = DATA.metadata.spb;
int = spb:spb:length(DATA.burst);
remd = length(DATA.burst);
intv = [1 int remd];
for j = 1:length(intv)-1 %burst counter
    for k = 1:DATA.nCells; %bin counter
        v1 = DATA.vel_b1(intv(j):intv(j+1),k);
        v2 = DATA.vel_b2(intv(j):intv(j+1),k);
        v3 = DATA.vel_b3(intv(j):intv(j+1),k);
        
        if strcmp(ornt,'UP') == 1
            v1=fliplr(v1);
            v2=fliplr(v2);
            v3=fliplr(v3);
        end
        
        %Unwrap using Julia's code
        vwrap=(max(v1(:))-min(v1(:)))*0.5;
        w1=unwrap_w_prof1(v1,vwrap);
        w1=fliplr(w1);
        vwrap=(max(v1(:))-min(v2(:)))*0.5;
        w2=unwrap_w_prof1(v2,vwrap);
        w1=fliplr(w1);
        vwrap=(max(v3(:))-min(v3(:)))*0.5;
        w3=unwrap_w_prof1(v3,vwrap);
        w1=fliplr(w1);
        
        ccrit = 70;
        bc1 = find(DATA.cor1(intv(j):intv(j+1),k)<=ccrit);
        bc2 = find(DATA.cor2(intv(j):intv(j+1),k)<=ccrit);
        bc3 = find(DATA.cor3(intv(j):intv(j+1),k)<=ccrit);
        w1(bc1) = NaN;w2(bc2) = NaN;w3(bc3) = NaN;
        w1 = cmgbridge(w1,nlin,maxbr,maxg);
        w2 = cmgbridge(w2,nlin,maxbr,maxg);
        w3 = cmgbridge(w3,nlin,maxbr,maxg);
        DATA.wvel_b1(intv(j):intv(j+1),k) = w1;
        DATA.wvel_b2(intv(j):intv(j+1),k) = w2;
        DATA.wvel_b3(intv(j):intv(j+1),k) = w3;
        clear w1 w2 w3 v1 v2 v3 bc1 bc2 bc3
    end
end

disp('Unwrapping complete')

%Rotate velocities to ENU, then to true N
[u,v,w] = aqdpbeam2enu(DATA.transfm,DATA.nCells,DATA.wvel_b1,DATA.wvel_b2,DATA.wvel_b3,heading,DATA.pitch,DATA.roll,ornt);
[u,v] = mag2truenorth(u,v,magdec);
DATA.u = u;DATA.v = v;DATA.w = w;

%run surface tracking if the instrument is pointing up
if strcmp(ornt,'UP')
    DATA = aqdpstrack(DATA,DATA.metadata.lat);
end

% fields = {'vel_b1','vel_b2','vel_b3','wvel_b1','wvel_b2','wvel_b3'};
% DATA = rmfield(DATA,fields);

%Pad bursts
DATA = padbursts(DATA);

aqdp = DATA;
clear DATA
save([dirname fname],'aqdp','-v7.3')
disp(['file ' name '.mat saved'])
function [burstData, puvData, aqdpMeta] = wad2puv(wadFile, whdFile, aqdpMeta, userMeta);
% wad2puv.m  A function to run PUV analysis on wave data from a Nortek
%            Aquadopp.
%
%   usage:  [burstData, puvData] = wad2puv(wadFile, aqdpMeta, userMeta);
%
%       where:  burstData - a structure with the following fields
%                       burst - a vector of burst numbers
%                          dn - time, ML datenum
%                           p - pressure, m
%                           u - east velocity, m/s
%                           v - north velocity, m/s
%                           w - vertical velocity, m/s
%                           a - amplitude, counts
%               puvData - a structure with the following fields
%                            burst - a vector of burst numbers 
%                               dn - time, ML datenum
%                                p - burst-mean pressure, m
%                             hsig - significant wave height, m
%                            peakF - wave frequency at the peak, Hz
%                          peakDir - wave direction at the peak, deg
%                       peakSpread - wave spreading at the peak, deg
%                               Su - surface elevation spectra based on
%                                    velocity, m^2/Hz
%                               Sp - surface elevation spectra based on
%                                    pressure, m^2/Hz
%                              Dir - wave direction, deg
%                           Spread - wave spreading
%                                F - center frequency of each frequency
%                                    band
%                               dF - bandwidth of each frequency band
%                wadFile - the name of the .wad file output by AquaPro
%                whdFile - the name of the .whd file output by AquaPro
%               aqdpMeta - a structure containing Aquadopp setup
%                          parameters
%               userMeta - a structure containing user-defined metadata
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Dependencies:
%   beam2enu.m
%   gmean.m
%   wds.m
%   hs.m

% C. Sullivan   10/27/05,   version 1.2
% Provide user additional feedback regarding code execution.  Add version
% number.
% C. Sullivan   06/06/05,   version 1.1
% Velocities collected in beam coordinates must be converted to geographic
% coordinates prior to PUV analysis. Alert user that PUV parameters in
% user-defined metadata file are being applied. If user enters an erroneous
% input (anything non-numeric including carriage return) for either
% 'first_good' or 'last_good', they are now asked to re-enter their answers.
% C. Sullivan   06/02/05,   version 1.0
% This function combines elements of George Voulgaris' function split_wad.m
% and puv_nortek.m. The idea is to read the Aquadopp's .wad file, split the
% pressure and velocity data into individual bursts, and perform PUV
% analysis burst by burst for only good bursts. It is assumed colums 8, 9,
% and 17 in the .wad file are empty. The pressure and velocity data is
% stored in arrays that are dimensioned [nSamples x nBursts]. The user is
% asked which bursts to perform the PUV analysis on. The mfiles for the PUV
% analysis are wds.m and hs.m, which were provided by Nortek. The data from
% PUV analysis is stored in arrays that are dimensioned [1 x nGoodBursts]
% for wave parameters and [nFreq x nGoodBursts] for wave spectra.


version = '1.2';

disp(' ')

% ---- User and istrument metadata -------------------------------------- %
nSamples = aqdpMeta.Wave__Number_of_samples;                 %number of samples
nums = ~isletter(aqdpMeta.Wave__Sampling_rate); 
sampleDt = 1 / str2num(aqdpMeta.Wave__Sampling_rate(nums)); %sampling rate, sec
nums = ~isletter(aqdpMeta.Wave__Interval); 
burstDt = str2num(aqdpMeta.Wave__Interval(nums));           %wave interval, sec
nums = ~isletter(aqdpMeta.Wave__Cell_size);
cell_ht = str2num(aqdpMeta.Wave__Cell_size(nums));          %height of velocity
                                                            % cell above bed, m
xducer_ht = userMeta.sensor_height;                         %height of sensor
                                                            % head above
                                                            % bed, m
coord_sys = aqdpMeta.Coordinate_system;                     %coordinate system
                                                            % of data
MagDev = userMeta.magnetic_variation*pi/180;                %magnetic variation,
                                                            % radians
                                                            
% ---- PUV parameters --------------------------------------------------- %
lf = userMeta.lf; %low frequency cutoff, Hz  
maxfac = userMeta.maxfac; %maximum value of factor scaling pressure to waves
minspec = userMeta.minspec; %minimum spectral level for computing direction and spreading, mm^2/Hz
Ndir = userMeta.Ndir; %direction offset, deg (includes compass error and misalignment of cable probe relative to case. (the offset for the Aquadopp Profiler is 0)
nf = userMeta.nf; %nominal number of output frequencies

parms=[lf maxfac minspec Ndir];

% ---- Load data -------------------------------------------------------- %
disp(['Loading Aquadopp pressure and velocity data from ',wadFile,...
      '. This may take a few minutes.'])
data = load(wadFile);
mo = data(:,1);
day = data(:,2);
yy = data(:,3);
hh = data(:,4);
mi = data(:,5);
ss = data(:,6);
p = data(:,7);   % meters
v1 = data(:,10);  % beam1|X|East velocity (m/s)
v2 = data(:,11);  % beam2|Y|North velocity (m/s)
v3 = data(:,12);  % beam3|Z|UP velocity (m/s)
a1 = data(:,14); %Amplitude Beam 1  (counts)
a2 = data(:,15); %Amplitude Beam 2  (counts)
a3 = data(:,16); %Amplitude Beam 3  (counts)

data = load(whdFile);
hdg = data(:,12);
ptch = data(:,13);
roll = data(:,14);

% ---------------- Define time ------------------------------------------ %
%we use a time-based method for defining bursts in case the # of bursts
%recorded to the .wad file does not equal the length(.wad file)/nSamples
dnO = datenum(yy,mo,day,hh,mi,ss);  %time stamp on samples, ML datenum
dnOff = datenum(yy(1),1,1,0,0,0);   
dn = dnO - dnOff;                   %<---  Why is this done?

clear nums data mo day yy hh mi ss 

% ---------------- Split data ------------------------------------------- %
%with a time-based method for defining bursts, the time between
%bursts should be > (burstDt - (nSamples*sampleDt)). this is only valid
%when (nSamples*sampleDt) < burstDt
time_btw_bursts = (burstDt - (nSamples * sampleDt)) / (3600*24); %days 

%Istart (a vector) indicates the indices in the .wad file
%at which bursts begin. however be sure to include 0 as the
%first value in Istart or else you will miss the first burst!
burst_start_ind = find(diff(dn) > time_btw_bursts); %<--- why not just use dnO?
burst_start_ind = [0; burst_start_ind];

%number of bursts is equal to the length of burst_start_ind
nBursts = length(burst_start_ind);

%sometimes instrument is stopped during a burst
if length(dn) < (burst_start_ind(end)+nSamples)
    nBursts = nBursts-1;
end

%pre-allocate the structure burstData. pressure, velocities, amplitudes,
%and dates will be written to this structure as fields that are
%dimensioned [nSamples x nBursts]
burstData.burst = nan*ones(1, nBursts);
burstData.dn = nan*ones(nSamples, nBursts);
burstData.p = nan*ones(nSamples, nBursts);
burstData.u = nan*ones(nSamples, nBursts);
burstData.v = nan*ones(nSamples, nBursts);
burstData.w = nan*ones(nSamples, nBursts);
burstData.a = nan*ones(nSamples, nBursts);

%loop through the bursts and assign all samples in each burst to
%burstData
for b=1:nBursts;
  range = [ burst_start_ind(b)+1 : burst_start_ind(b)+nSamples ];
  burstData.burst(1,b) = b;
  burstData.dn(:,b) = dnO(range);
  burstData.p(:,b) = p(range);
  
  vv1 = v1(range); %all velocity samples in burst b
  vv2 = v2(range); 
  vv3 = v3(range); 
    
  switch coord_sys
      case 'BEAM'
          %if data was collected in BEAM coordinates, convert data to
          %geographic coordinates (and correct for magnetic variation).
          %do the conversion for each burst, sample by sample
          if mod(b,10)==0
              disp(['Converting velocities in burst ',num2str(b),' from BEAM to ',...
                    'ENU coordinates, and correcting for declination. ',...
                    num2str(toc/60),' minutes elapsed'])
          end
          for ii = 1:nSamples
              vbeam = [vv1(ii) vv2(ii) vv3(ii)];
              venu=beam2enu(vbeam,hdg(b),ptch(b),roll(b),0);
              East(ii) = venu(1); %East
              North(ii) = venu(2); %North
              Up(ii) = venu(3); %Up
          end
              
              %correct for magnetic variation
              East_mag =  East*cos(MagDev) + North*sin(MagDev);
              North_mag = -East*sin(MagDev) + North*cos(MagDev);
              Up_mag = Up;
              
              burstData.u(:,b) = East_mag';
              burstData.v(:,b) = North_mag';
              burstData.w(:,b) = Up_mag';
      otherwise
          burstData.u(:,b) = vv1; %East
          burstData.v(:,b) = vv2; %North
          burstData.w(:,b) = vv3; %Up
          
  end
              
    burstData.a(:,b) = a1(range)+a2(range)+a3(range);
  
end

%clear variables you no longer need
clear p v1 v2 v3 a1 a2 a3 

% ---- Define good bursts to include in PUV analysis -------------------- %
%display a figure with burst-averaged pressure and velocities.  ask the
%user if, based on the figure, they would like to include all bursts in the
%PUV analysis.  if the answer is no, ask the user to enter the the first
%and last good bursts. all bursts between the first and last good bursts
%will be included in PUV analysis.
figure
subplot(211)
plot(burstData.burst, gmean(burstData.p))
ylabel('Burst-mean pressure (m)')
xlabel('Burst number')
subplot(212)
plot(burstData.burst, gmean(burstData.u))
hold on
plot(burstData.burst, gmean(burstData.v),'r')
plot(burstData.burst, gmean(burstData.w),'g')
ylabel('Burst-mean Velocity (m/s)')
xlabel('Burst number')
legend('Eastward Velocity','Northward Velocity','Vertical Velocity')
disp(' ')
disp('Please view the figure (feel free to zoom!) and look for bursts which might be bad')
disp('Bad bursts might include out-of-water bursts, for example.')
answer = input(['Would you like to include all bursts (bursts 1:',...
              num2str(nBursts),') in PUV analysis? y/n :  '],'s');
if ~strcmp(answer,'Y') & ~strcmp(answer,'y') & ...
   ~strcmp(answer,'N') & ~strcmp(answer,'n')
        disp(['Valid answers are either yes ("y") or no ("n")'])
        answer = input(['Would you like to include all bursts (bursts 1:',...
              num2str(nBursts),') in PUV analysis? y/n :  '],'s');
   
end
if strcmp(answer,'Y') || strcmp(answer,'y')
    first_good = 1;
    last_good = nBursts;
    disp(['Including bursts 1 through ',num2str(nBursts),' in PUV analysis'])
elseif strcmp(answer,'N') || strcmp(answer,'n')
    first_good = input('Please view the figure and enter the first good burst:  ');
    if isempty(first_good)
        disp('Valid answers are numeric only')
        first_good = input('Please view the figure and enter the first good burst:  ');
    end
    last_good = input('Please view the figure and enter the last good burst:  ');
    if isempty(last_good)
        disp('Valid answers are numeric only')
        last_good = input('Please view the figure and enter the last good burst:  ');
    end
    disp(['Including bursts ',num2str(first_good),' to ',num2str(last_good), ...
          ' in PUV analysis'])
    disp(' ')
end
goodBursts = [first_good:last_good];
nGoodBursts = length(goodBursts);

% ---- Do PUV Analysis -------------------------------------------------- %
disp(' ')
disp(['Performing PUV analysis.  This may take a few minutes'])
disp(' ')
disp('Users should have specified PUV parameters that suit their ')
disp('project needs in their metadata file.')
disp(['PUV parameter values defined the metadata file include:'])
disp(['     lf = ',num2str(lf),' Hz'])
disp(['     maxfac = ',num2str(maxfac)])
disp(['     minspec = ',num2str(minspec),' mm^2/Hz'])
disp(['     Ndir = ',num2str(Ndir),' deg'])
disp(' ')
%pre-allocate the structure puvData. results from the PUV analysis will
%will be written to this structure as fields
puvData.burst = nan*ones(1, nGoodBursts);
puvData.dn = nan*ones(1, nGoodBursts);
puvData.p = nan*ones(1, nGoodBursts);
puvData.Hsig = nan*ones(1, nGoodBursts);
puvData.peakF = nan*ones(1, nGoodBursts);
puvData.peakDir = nan*ones(1, nGoodBursts);
puvData.peakSpread = nan*ones(1, nGoodBursts);
puvData.Su = nan*ones(nf, nGoodBursts); %<--- use nF but the # of output freq.
                                   %     is not necessarily nF.
puvData.Sp = nan*ones(nf, nGoodBursts);
puvData.Dir = nan*ones(nf, nGoodBursts);
puvData.Spread = nan*ones(nf, nGoodBursts);
puvData.F = nan*ones(nf, nGoodBursts);
puvData.dF = nan*ones(nf, nGoodBursts);

for b=1:nGoodBursts;
    
   burst = goodBursts(b);
   if mod(b,10) == 0
       disp(['Finished PUV analysis on burst ',num2str(burst),'. ',...
              num2str(toc/60),' minutes elapsed'])
   end
   dn = gmean(burstData.dn(:,burst)); %set time to time at center of burst
   p = gmean(burstData.p(:,burst));
   
   [Su,Sp,Dir,Spread,F,dF] = WDS(burstData.u(:,burst),burstData.v(:,burst),...
                                 burstData.p(:,burst),sampleDt,nf,xducer_ht,...
                                 cell_ht,parms);
                             
   [Hsig,peakF,peakDir,peakSpread] = HS(Su,Sp,Dir,Spread,F,dF);
   
   theVars = fieldnames(puvData);
          
   for v=1:length(theVars)
       eval(['[nf_out, col] = size(',theVars{v},');'])
       if nf_out == 1
           eval(['puvData.',theVars{v},'(1,b) = ',theVars{v},';'])
       else
           eval(['puvData.',theVars{v},'(1:nf_out,b) = ',theVars{v},';'])
           if b == nGoodBursts & nf_out < nf
                %remove extra rows (frequencies) if the number of
                %output frequencies is less than the nominal number
                %of frequencies (nf)
                eval(['puvData.',theVars{v},'(nf_out+1:end,:) = [];'])
                %add nf_out to instrument metatdata.  this value will be
                %used to dimension NetCDF variables
                aqdpMeta.nf_out = nf_out;
           end 
       end
   end

end

disp(['Finished PUV analysis. ',num2str(toc/60),' minutes elapsed'])

return


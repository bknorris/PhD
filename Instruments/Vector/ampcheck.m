function ADV = ampcheck(ADV,metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For Nortek Vectors: perform a check using the amplitudes to determine
%   points below a specified cutoff amplitude
%
%   ADV - a structure ADV containing velocities and sensor data
%   Metadata - user defined metadata. For a  complete list of required 
%   metadata fields, see adv_dataprocess.m
%
%   Contains a partial adaptation of vec2ncBurstNative.m
%   Written by Kurt J Rosenberger
%   krosenberger@usgs.gov
%   USGS Pacific Coastal Marine Science Center
%
%   This script was compiled by Benjamin K Norris, 2015
%   University of Waikato, New Zealand
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check amplitudes - Nortek lists the noise floor as 70 counts for the
%vector


disp('Screening data by predetermined Amplitude cutoff for this frequency system....')
if isfield(metadata,'amplitude_cutoff')
    cutoff = metadata.amplitude_cutoff;
    disp(['User selected an amplitude cutoff of ',num2str(cutoff),' for screening data'])
else
    cutoff = 70;
end
disp(['Amplitude of ',num2str(cutoff),' counts used for data screening']);

ind = find(ADV.Amp1<cutoff | ADV.Amp2<cutoff | ADV.Amp3<cutoff);
ADV.U(ind) = NaN;
ADV.V(ind) = NaN;
ADV.W(ind) = NaN;

% calculate the statistics; for time, we will use the middle of the burst
flds = fieldnames(ADV);
for i = 1:length(flds)
    fld = flds{i};
    if length(ADV.(fld))>1; %i.e. only work on timeseries
        STATS.(fld) = nanmean(ADV.(fld));
        if any(strcmp(fld,{'U';'V';'W';'pressure'}))
            stdnm = [fld '_STD'];STATS.(stdnm) = nanstd(ADV.(fld))';
            maxnm = [fld '_MAX'];STATS.(maxnm) = nanmax(ADV.(fld))';
            minnm = [fld '_MIN'];STATS.(minnm) = nanmin(ADV.(fld))';
        end
    end
end
ADV.Stats = STATS;
end
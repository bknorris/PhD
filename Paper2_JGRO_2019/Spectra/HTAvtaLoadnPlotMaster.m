%Load and process files from the HTA and VTA. Hoo boy this is gonna be a
%big one.

%you're locked on to BBC Radio One and One Xtra

clear
close all
savedatdir = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';

%first run the overview plots script
if 0
    plotPressureRecord
end

%load the data, save into structures in a file location
if 1
    load([savedatdir 'HTAtimes.mat'])
    load([savedatdir 'VTAtimes.mat'])
    LoadnQCdataHTAvta
end

%Generate slope t-s from the turbulent f bands of spectra
if 1
    HTASlopeFromSpec
    VTASlopeFromSpec
end

%run the spectral analysis
if 1
%     RunSpectra
%     AnimateSlopeFromSpec
    EvolutionarySpectra
end
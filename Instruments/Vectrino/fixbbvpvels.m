function Data = fixbbvpvels(Data,rmv)
%function to remove below-bed vectrino velocities for instruments  that 
%were mounted close to bed level. Removal is completed by applying the 
%bed level bottom check estimate. Any velocity falling below this value 
%for a given time stamp is set to NaN.
%
% Inputs: the Data structure of any VecPro file. rmv is a flag to denote
% the removal type. if rmv = 1, set all below-bed velocities to NaN. If
% rmv = 2, set all below-bed velocities to 0.
% Outputs: the Data structure with velocity fields cropped to bed level.
%
% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand
if rmv == 1
    rmf = NaN;
elseif rmv == 2
    rmf = 0;
end

disp('User requests sub bed-level velocity measurements to be removed')
%first need to refine the bottom dist parameter
bdist = Data.BottomCheck_BottomDistance;
bdmax = runningmax(bdist,100);
bdmax = my_running_median(bdmax,500); %despike
bdmax = smooth(bdmax,100,'sgolay');
%need to make bdmax the same length as the t-s
bd = spline(Data.BottomCheck_HostTimeMatlab,bdmax,Data.Profiles_HostTimeMatlab);
if length(Data.Time) ~= length(Data.Profiles_VelBeam1)
    bd = bd(1:length(Data.Profiles_VelBeam1));
end
dhts = repmat(Data.Profiles_Range,length(bd),1);

for ii = 1:length(bd)
    indx = dhts(ii,:) >= bd(ii);
    Data.Profiles_VelBeam1(ii,indx) = rmf;
    Data.Profiles_VelBeam2(ii,indx) = rmf;
    Data.Profiles_VelBeam3(ii,indx) = rmf;
    Data.Profiles_VelBeam4(ii,indx) = rmf;
end
disp('Sub bed-level velocity measurements have been removed')
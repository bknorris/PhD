%since I can't program in R, I'll need to organize the R output data. I've
%saved the 'mss' variable from the photogrammetric statistical analysis as
%a .mat file called "RENAMEstats.mat". This should be renamed to the
%Quadrat number per deployment. The .mat file contains fields "layer#"
%that correspond to the number of elevation slices per point cloud
%analysis. This code is designed to read these in, concatenate them into a
%structure, and save out a properly formatted file that can be read into
%Matlab.

clear
ddir = 'd:\Projects\Mekong_F2014\Images\Quadrats\Reconstructions\Q5.10\';
name = 'Q5_10';
l = 0.005; %step length between layers used during the R GUI step in m (default is 0.01)
scalar = 'Y'; %[Y/N] set to 'Y' if the quadrat dimensions were not 1x1 (m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([ddir name 'stats.mat'],'file')
    movefile([ddir 'RENAMEstats.mat'],[ddir name 'stats.mat'])
end
load([ddir name 'stats.mat'])

%get the names of the variables for concatenation
data = who('layer*');n = length(data);
DATA = struct();
idx = zeros(n,3);
%layers 0-9 have only six alphanumeric digits
%layers 10-99 have seven alphanumeric digits
%layers 100-999 have eight alphanumeric digits
for i = 1:n
    if length(char(data{i})) == 6
        idx(i,1) = i;
    elseif length(char(data{i})) == 7
        idx(i,2) = i;
    elseif length(char(data{i})) == 8
        idx(i,3) = i;
    end
end
fname = {'Rradii';'ID';'Layer';'Layerminh';'Layermaxh';'Rcenter';...
    'Clusterlabel';'Closestlabel';'Cbandw';'Datapts';'Scaled';'Scaledby'};
DATA.info = {'Rradii is a matrix of size N x 8, containing the 8 circular radii for each of the N roots';...
'ID is a vector of length N containing a unique ID for each root';...
'Layer represents the number of the current layer. It is 1 for layer1, 2 for layer2, etc';...
'Layerminh and Layermaxh are the lower and upper heights of the layer, respectively, in meters';...
'Rcenter is an N x 2 matrix that contains the coordinates of the root centers';...
'Cbandw is the bandwidth used for the mean-shift clustering';...
'Datapts are the data points located within the layer';...
'Scaled and Scaledby are the scaling factors that are used if you specified a quadrat that was not a 1x1 meter'};
l = 0:l:l*n;
%for layers1-9
counter = 1;
for i = 1:3
    for j = 1:n
        if idx(j,i) == 0 %zeros exist where there aren't indexes, skip these and continue to the next iteration
            continue;
        else
            tmp1 = eval(data{j});
            oldfname = fieldnames(tmp1);
            for ii = 1:length(oldfname)
                if isempty(tmp1.(oldfname{ii})) %R isn't always saving some variables out,
                    %so this section replaces this missing info with the correct info
                    if strcmp(oldfname{ii},'cluster.slice')
                        [tmp2.(fname{ii})] = counter;
                    elseif strcmp(oldfname{ii},'cluster.hm')
                        [tmp2.(fname{ii})] = l(counter);
                    elseif strcmp(oldfname{ii},'cluster.hM')
                        [tmp2.(fname{ii})] = l(counter+1);
                    elseif strcmp(oldfname{ii},'scaled')
                        [tmp2.(fname{ii})] = scalar;
                    end
                else
                    [tmp2.(fname{ii})] = tmp1.(oldfname{ii}){:};
                end
            end
            counter = counter+1;
            DATA.(data{j}) = tmp2;
            clearvars tmp1 tmp2
        end
    end
end
clearvars -except DATA ddir nn name
nname = [name 'stats_fixed.mat'];
save([ddir nname],'DATA','-mat','-v7.3')

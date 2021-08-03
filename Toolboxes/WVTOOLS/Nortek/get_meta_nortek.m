function [userMeta, aqdpMeta] = get_meta_nortek(metaFile, hdrFile);
% get_meta_nortek.m  A function to load user-defined metadata and
%                    instrument setup parameters for a Nortek Aquadopp
%                    Profiler (AP).
%
%   usage:  [userMeta, aqdpMeta] = get_meta_nortek(metaFile, hdrFile);
%
%       where:  metaFile - a string specifying the ascii file in which
%                          metadata is defined, in single quotes
%                          excluding the .txt file extension
%               hdrFile - a string specifying the ascii file in which
%                         Aquadopp setup parameters are defined, in single
%                         quotes including the extension .hdr
%               userMeta - a structure with user-defined metadata
%               aqdpMeta - a structure with Aquadopp setup parameters
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   10/27/05,   version 1.1
% Metadata file must be a .txt file.  The .txt extension is not required in
% the input 'metaFile'. Provide the user additional feedback regarding code
% execution.
% C. Sullivan   06/02/05,   version 1.0
% This function is coded to read a user-defined metadata file with specific
% formatting.  See the toolbox documentation for more information.


version = '1.1';

% Gather user-defined metadata
disp(' ')
disp('Reading user metadata')
[atts, defs] = textread([metaFile,'.txt'],'%s %63c','commentstyle','shell');
defs = cellstr(defs);
for i = 1:length(atts)
    theAtt = atts{i}(:)';
    theDef = defs{i}(:)';
    %deblank removes trailing whitespace
    theAtt = deblank(theAtt);
    theDef = deblank(theDef);
    %check for and replace spaces in
    %the attributes with underscores
    f1 = find(isspace(theAtt));
    f2 = strfind(theAtt,'-');
    f = union(f1,f2);
    if ~isempty(f)
        theAtt(f) = '_';
    end
    %attribute definitions read in as characters; convert to
    %numbers where appropriate
    theDefNum = str2num(theDef);
    if ~isempty(theDefNum) 
        theDef = theDefNum;
        eval(['userMeta.',theAtt,' = theDef;'])
    else
        eval(['userMeta.',theAtt,'=''',theDef,''';'])
    end
end


% Gather instrument setup parameters
disp('Reading Aquadopp setup parameters')
fid = fopen(hdrFile,'r');
while 1
    tline = fgetl(fid);
    if ~isempty(tline) & ~strcmp(tline(1),'[') & ~strcmp(tline(1),'-') & ...
       ~strcmp(tline,'User setup') & ~strcmp(tline,'Hardware configuration') & ...
       ~strcmp(tline,'Head configuration') & ~strcmp(tline,'Data file format')
        %deblank removes trailing whitespace
        theAtt = deblank(tline(1:38));
        theDef = deblank(tline(39:end));
        %check for and replace spaces in
        %the parameters with underscores
        f = find(isspace(theAtt));
        if ~isempty(f)
            theAtt(f) = '_';
        end
        %check for and replace dashes in
        %the parameters with underscores
        f = strfind(theAtt,'-');
        if ~isempty(f)
            theAtt(f) = [];
        end
        %parameter values read in as characters; convert to
        %numbers where appropriate
        theDefNum = str2num(theDef);
        if ~isempty(theDefNum) & ~strcmp(theAtt(1:4),'Time')
            theDef = theDefNum;
            eval(['aqdpMeta.',theAtt,' = theDef;'])
        else
            eval(['aqdpMeta.',theAtt,'=''',theDef,''';'])
        end
    end
    %only read up to data file format in the .hdr file
    if strcmp(tline,'Head configuration'), break, end
end
fclose(fid);

return
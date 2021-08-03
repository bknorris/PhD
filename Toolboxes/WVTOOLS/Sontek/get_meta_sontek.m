function [userMeta, argnMeta] = get_meta_sontek(metaFile, ctlFile);
% get_meta_sontek.m  A function to load user-defined metadata and 
%                    instrument setup parameters from a Sontek Argonaut.
%
%   usage:  [userMeta, argnMeta] = get_meta_sontek(metaFile, ctlFile);
%
%       where:  metaFile - a string specifying the ascii file in which
%                          metadata is defined, in single quotes
%                          excluding the .txt file extension
%               ctlFile - a string specifying the ascii file in which
%                         Argonaut setup parameters are defined, in single
%                         quotes including the file extension .ctl
%                         which contains Argonaut setup parameters
%               userMeta - a structure with user-defined metadata
%               argnMeta - a structure with Argonaut setup parameters
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   11/02/05,   version 1.1
% Metadata file must be a .txt file.  The .txt extension is not required in
% the input 'metaFile'. Provide the user additional feedback regarding code
% execution.
% C. Sullivan   05/31/05,   version 1.0
% This function is coded to read a user-defined metadata file with specific
% formatting.  See the toolbox documentation for more information. It is
% assumed that the comments in the Argonaut's .hdr file occupy no more than
% 3 lines.


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
disp('Reading Argonaut setup parameters')
fid = fopen(ctlFile,'r');
while 1
    tline = fgetl(fid);
    %comments are tricy to read. assume they
    %use 3 lines for now
    if strcmp(tline,'Comments:')
        theDef=[fgetl(fid),'; ',fgetl(fid),'; ',fgetl(fid)];
        eval(['argnMeta.Comments =''',theDef,''';'])
        break
    end
    if ~isempty(tline) & ~strcmp(tline(1:8),'Argonaut') & ~strcmp(tline(1),'-')
        theAtt = tline(1:28);
        theDef = tline(29:end);
        %if the parameter contains units in parenthesis, the
        %units will be written as theAtt.units = units
        fs = strfind(theAtt,'(');
        fe = strfind(theAtt,')');
        if ~isempty(fs) & ~isempty(fe)
            theUnits = theAtt(fs+1:fe-1);
            theAtt=regexprep(theAtt,'\([^)]*\)','');
        end
        %check for and remove dashes in
        %the parameters
        f = strfind(theAtt,'-');
        if ~isempty(f)
            theAtt(f) = [];
        end
        %check for and remove slashes in
        %the parameters
        f = strfind(theAtt,'/'); 
        if ~isempty(f)
            theAtt(f)= [];
        end
        %deblank removes trailing whitespace
        theAtt = deblank(theAtt);
        theDef = deblank(theDef);
        %replace whitespace with underscores in parameter
        %names
        f = isspace(theAtt);
        theAtt(f) = '_';
        %transformation matrix is tricky to read in
        if strcmp(theAtt,'Transformation_Matrix')
            temp = fscanf(fid,'%6f',[3,2]);
            temp = temp';
            temp = num2str(temp);
            theDef = strvcat(theDef, temp);
        end
        %parameter values read in as characters; convert to
        %numbers where appropriate
        theDefNum = str2num(theDef);
        if ~isempty(theDefNum) & isempty(strfind(theAtt,'Time'))
            %parameter value is numeric
            eval(['argnMeta.',theAtt,' = theDefNum;'])
        else
            %parameter value is ascii
            eval(['argnMeta.',theAtt,'=''',theDef,''';'])
        end
        if exist('theUnits')
            eval(['argnMeta.',theAtt,'Units=''',theUnits,''';'])
            clear theUnits
        end
    end
end
fclose(fid);

return
%Determine the pneumatophore statistics (d,n,a and phi) as cylinders for a
%select number of pneumatophores in one of the Quadrats

%Updates: 08/06/2016 redesigned as a utility to run through a selection of
%vegetation stats files and calculate specific statistics using user input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

disp('Running vegetation statistics calculator')
disp('Select a working directory that contains vegetation statistics')
folder = uigetdir;
fname = dir([folder '/' '*stats_fixed.mat']);fname = {fname.name};
vegdat = struct();

%Load the VP positions data
vpdir = '/Volumes/NORRIS/Data/DataReports/';
load([vpdir 'VP_positions.mat'])
if strfind(folder,'2014') > 1
    VPpos.Qname = VPpos.Qname(1:21);
    VPpos.VP_X = VPpos.VP_X(1:21);
    VPpos.VP_Y = VPpos.VP_Y(1:21);
    year = 1;
elseif strfind(folder,'2015') > 1
    VPpos.Qname = VPpos.Qname(22:end);
    VPpos.VP_X = VPpos.VP_X(22:end);
    VPpos.VP_Y = VPpos.VP_Y(22:end);
    year = 2;
end
%%%Rotate quadrat so the top of fig points N
if year == 1
    rotx = [0 0 180 0 0 0 0 0 0 0 0 0 0 0 180 0 0];
    roty = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 180 0];
    rotz = [90 90 0 -90 0 0 0 0 0 0 0 0 0 0 90 90 0];
elseif year == 2
    rotx = [0 180 180 0 0 180 0 180];
    roty = [0 0 0 0 180 0 0 180];
    rotz = [90 0 0 270 90 0 -90 0];
end
for k = 1:length(fname)
    %%%Load the Data
    disp(['Loading ' folder '/' fname{k}])
    load([folder '/' fname{k}])
    qname = regexprep(fname{k},'stats_fixed.mat','');
    DATA = rmfield(DATA,'info'); %remove info field for processing
    
    %%%Load the correct VP data and scatterplot
    vpid = zeros(length(VPpos.Qname),1);
    for b = 1:length(VPpos.Qname)
        vpid(b) = strcmp(VPpos.Qname{b},qname);
    end
    vpx = (VPpos.VP_X.*vpid)';vpx(vpx == 0) = [];
    vpy = (VPpos.VP_Y.*vpid)';vpy(vpy == 0) = [];
    %%%Check to see if the selection is correct with the 'while' statement
    yes = 0;
    while yes ~= 1
        %%%First, plot the entire quadrat as an overview
        f1 = figure(1);
        set(f1,'PaperOrientation','portrait',...
            'position',[200 300   800   700]);
        
        %Create rotation matrix, rotate about z-axis if specified, flip
        %axes if specified
        thx = rotx(k)*pi/180;thy = roty(k)*pi/180;thz = rotz(k)*pi/180;
        Tx = [cos(2*thx) sin(2*thx); sin(2*thx) -cos(2*thx)];
        Ty = [-cos(2*thy) sin(2*thy); sin(2*thy) cos(2*thy)];
        Tz = [cos(thz) -sin(thz); sin(thz) cos(thz)];
        
        %Adjust VP positions based on rotations and reflections
        o = repmat([0.5;0.5],1,length(vpx)); %origin
        vpxy = [vpx;vpy];        
        if rotz(k) ~= 0
            vpxy = Tz*(vpxy-o)+o;
        end
        if rotx(k) ~= 0
            vpxy = Tx*(vpxy-o)+o;
        end
        if roty(k) ~= 0
            vpxy = Ty*(vpxy-o)+o;
        end

        %Plot VP Positions
        h = plot(vpxy(1,:),vpxy(2,:),'^','MarkerSize',12);
        set(h,'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 0.9 0])
        hold on
        l = length(DATA.layer1.Rcenter);
        x = DATA.layer1.Rcenter(:,1)';
        y = DATA.layer1.Rcenter(:,2)';
        
        %Rotate & Reflect coordinates
        o = repmat([0.5;0.5],1,length(x));
        xy = [x;y];        
        if rotz(k) ~= 0
            xy = Tz*(xy-o)+o;
        end
        if rotx(k) ~= 0
            xy = Tx*(xy-o)+o;
        end
        if roty(k) ~= 0
            xy = Ty*(xy-o)+o;
        end

        for i = 1:l
            xc = xy(1,i);
            yc = xy(2,i);
            r = mean(DATA.layer1.Rradii(i,:));
            ang=0:0.01:2*pi;
            xp=r*cos(ang);
            yp=r*sin(ang);
            p = plot(xc+xp,yc+yp);set(p,'linewidth',1.5)
            hold on
        end
        axis equal
        grid on
        set(gca,'Xlim',[0 1],'Ylim',[0 1],'grid',':')
        xlabel('X'),ylabel('Y')
        title(qname)
        hold off

        %%%Draw a rectangle around the desired points (for help: 'help rbbox')
        jj = waitforbuttonpress;
        point1 = get(gca,'CurrentPoint');
        finalRect = rbbox;
        point2 = get(gca,'CurrentPoint');
        point1 = point1(1,1:2);
        point2 = point2(1,1:2);
        p1 = min(point1,point2);
        offset = abs(point1-point2);
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        hold on
        axis manual
        plot(x,y,'r','linewidth',2)
        x = unique(x);xlo = x(1);xhi = x(2);
        y = unique(y);ylo = y(1);yhi = y(2);
        
        %%%Crop DATA layer1 based on answer inputs, then replot
        xid = find(xy(1,:) >= xlo & xy(1,:) <= xhi);
        yid = find(xy(2,:) >= ylo & xy(2,:) <= yhi);
        id = intersect(xid,yid);
        
        figure(2)
        %%%Plot VP Positions
        h = plot(vpxy(1,:),vpxy(2,:),'^','MarkerSize',12);
        set(h,'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 0.9 0])
        hold on
        l = length(id);
        R = zeros(l,1);
        X = zeros(l,1);
        Y = zeros(l,1);
        for i = 1:l
            xc = xy(1,id(i));
            yc = xy(2,id(i));
            r = mean(DATA.layer1.Rradii(id(i),:));
            ang=0:0.01:2*pi;
            xp=r*cos(ang);
            yp=r*sin(ang);
            p = plot(xc+xp,yc+yp);set(p,'linewidth',1.5)
            hold on
            R(i,:) = r;
            X(i,:) = xc;
            Y(i,:) = yc;
        end
        axis equal
        set(gca,'Xlim',[xlo xhi],'Ylim',[ylo yhi])
        xlabel('X'),ylabel('Y')
        title([qname ' Cropped'])
        hold off
        
        %%%Check with the user to affirm the selection is correct
        prompt = 'Is the selection correct? [y/n] ';
        str = input(prompt,'s');
        if strcmp(str,'y')
            yes = 1;
        elseif strcmp(str,'n')
            yes = 0;
            close(figure(1)),close(figure(2))
        end
    end
    
    %distance formulation: find the smallest distance between all of the
    %points, then take the mean
    dist = zeros(length(X),1);
    if length(X) > 1
        for i = 1:length(X)
            x1 = X(i);y1 = Y(i);
            xx = setxor(X,x1);yy = setxor(Y,y1);
            dd = sqrt((xx-x1).^2+(yy-y1).^2);
            dist(i) = min(dd);
        end
        deltaS = mean(mean(dist));
    else
        deltaS = 0;
    end
    %really small 'pneumatophores' are biasing the Phi estimates. Eliminate
    %radii smaller than 4mm for the averages.
    r = R(R >= 4E-3);
    D = mean(r)*2;
    n = numel(id);
    a = n*D;
    if length(X) > 1
        phi = (pi/4)*(D/deltaS)^2;
    else
        phi = (n*pi*(D^2))/4;
    end
    
    %%%Report Results, save results to a structure
    disp(['Average pneumatophore spacing (Delta S): ' num2str(deltaS)])
    disp(['Number of stems in first slice (n): ' num2str(length(id))])
    disp(['Mean stem diameter of first slice (d): ' num2str(D)])
    disp(['Frontal area density of vegetation (a): ' num2str(a)])
    disp(['Volume fraction occupied by vegetation (phi): ' num2str(phi)])
    
    vegdat.name(k,1:length(qname)) = qname;
    vegdat.n(k) = n;
    vegdat.d(k) = D;
    vegdat.deltaS(k) = deltaS;
    vegdat.phi(k) = phi;
    vegdat.x(k,1:2) = [xlo xhi];
    vegdat.y(k,1:2) = [ylo yhi];
    close(figure(1)),close(figure(2))
    clear R r D n a phi X Y
end

%%%Save Vegdat file to folder
titletxt = 'Save File as:';
prompt = 'Enter file name';
sfile = inputdlg(prompt,titletxt,[1, length(titletxt)+20]);
save([folder '\' sfile{:}],'vegdat','-v7.3')
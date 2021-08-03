%Plot U,W velocity pairs to show vertical movement of water across the
%Vectrino array. Will need to load each Vectrino file, then delete the
%unnecesary information 

clear
filename = 'FSS3_Q030715.gif';
%this will take some time to run initially
tic
load VP1_070315.mat
%convert to XYZ coordinates
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
%save out X,Z pairs
V1.X = VPRO.Data.Profiles_VelX;
V1.Z1 = VPRO.Data.Profiles_VelZ1;
V1.Z2 = VPRO.Data.Profiles_VelZ2;
sr = VPRO.Config.sampleRate;
clearvars VPRO

load VP2_070315.mat
%convert to XYZ coordinates
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
%save out X,Z pairs
V2.X = VPRO.Data.Profiles_VelX;
V2.Z1 = VPRO.Data.Profiles_VelZ1;
V2.Z2 = VPRO.Data.Profiles_VelZ2;
clearvars VPRO

load VP3_070315.mat
%convert to XYZ coordinates
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
%save out X,Z pairs
V3.X = VPRO.Data.Profiles_VelX;
V3.Z1 = VPRO.Data.Profiles_VelZ1;
V3.Z2 = VPRO.Data.Profiles_VelZ2;
clearvars VPRO
disp(['Finished loading data in ' num2str(toc/60) ' minutes'])

%compute depth averages
[m,n] = size(V1.X);  
for i = 1:m
    V1.Xav(i,1) = nanmean(V1.X(i,1:n));
    V1.Z1av(i,1) = nanmean(V1.Z1(i,1:n));
    V1.Z2av(i,1) = nanmean(V1.Z2(i,1:n));

    V2.Xav(i,1) = nanmean(V2.X(i,1:n));
    V2.Z1av(i,1) = nanmean(V2.Z1(i,1:n));
    V2.Z2av(i,1) = nanmean(V2.Z2(i,1:n));
    
    V3.Xav(i,1) = nanmean(V3.X(i,1:n));
    V3.Z1av(i,1) = nanmean(V3.Z1(i,1:n));
    V3.Z2av(i,1) = nanmean(V3.Z2(i,1:n));
end
%average Zs
V1.Zav = (V1.Z1av+V1.Z2av)./2;
V2.Zav = (V2.Z1av+V2.Z2av)./2;
V3.Zav = (V3.Z1av+V3.Z2av)./2;

%window
intv = 5; %averaging window in minutes
avt = 60*intv*sr;
ind = avt:avt:m;    
ind = [1 ind m]; %include remainder
p = length(ind)-1;
for i = 1:p
    idx = ind(i):ind(i+1);
    V1.Xdt(i,1) = nanmean(V1.Xav(idx)); %'dt' for depth-time average
    V1.Zdt(i,1) = nanmean(V1.Zav(idx));

    V2.Xdt(i,1) = nanmean(V2.Xav(idx));
    V2.Zdt(i,1) = nanmean(V2.Zav(idx));
    
    V3.Xdt(i,1) = nanmean(V3.Xav(idx));
    V3.Zdt(i,1) = nanmean(V3.Zav(idx));
end
c = jet(p-1);
xs = [-10 10 20]; %location in cm of vectrinos
ys = 1;
scale = 100; %scaling factor for arrows
%lets try something here
for i = 1:p
    f1 = figure;
    set(f1,'PaperOrientation','portrait',...
            'position',[400 200   700   450]);
    set(gcf,'color','w','PaperPositionMode','auto'),hold on
    q(1) = quiver(xs(1),ys,V1.Xdt(i)*scale,V1.Zdt(i)*scale);
    q(2) = quiver(xs(2),ys,V2.Xdt(i)*scale,V2.Zdt(i)*scale);
    q(3) = quiver(xs(3),ys,V3.Xdt(i)*scale,V3.Zdt(i)*scale);
    %gif commands
    drawnow
    frame = getframe(f1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256,'nodither');
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    close(f1)
end
    
    
    

%Calculate volumes from pneumatophore point cloud reconstructions
%V_cyl = pi*r^2*z, where r is the radius of the individual pneumatophore,
%and z is the height, which increases from 0:0.005:z.
clear, close all
ddir = 'd:\Projects\Mekong_F2014\Images\Quadrats\Reconstructions\Q4\Analysis\NoQuad\';
fdir = 'd:\Projects\Mekong_F2014\Figures\Quadrats\';
load([ddir 'Q4stats_fixed.mat'])
DATA = rmfield(DATA,'info'); %remove info field for processing
h = 0.005; %spacing between layers (m)
fn = fieldnames(DATA);n = length(fn);
H = h:h:h*n; %vector of heights to max H

%(x,y,z) points of VecPros in the Quadrat, measured in CloudCompare
%[VP1 VP2 VP3]
% vpx = [0.8033 0.5997 0.4978;0.8033 0.5997 0.4978;0.8033 0.5997 0.4978];
% vpy = [0.245 0.245 0.245;0.245 0.245 0.245;0.254 0.254 0.254;];
% vpz = [0.062 0.063 0.061;0.25 0.25 0.25;0.46 0.46 0.46]; %relative to the top of the quadrat, where the model bed level is initialized

%Rotate the coordinates
rotx = 0;
roty = 0;
rotz = -90;

%Create rotation matrix, rotate about z-axis if specified, flip
%axes if specified
thx = rotx*pi/180;thy = roty*pi/180;thz = rotz*pi/180;
Tx = [cos(2*thx) sin(2*thx); sin(2*thx) -cos(2*thx)];
Ty = [-cos(2*thy) sin(2*thy); sin(2*thy) cos(2*thy)];
Tz = [cos(thz) -sin(thz); sin(thz) cos(thz)];
        
%plot cylinders as volumes to show the process for estimating veg. volumes
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
%generate 'ground' and borders in figure
patch([0 1.5 1.5 0],[0 0 1.5 1.5],[0.85 0.85 0.85])
% patch([0 0 0 0],[1 0 0 1],[0 0 1 1],[0.85 0.85 0.85])
% patch([1 0 0 1],[1 1 1 1],[0 0 1 1],[0.85 0.85 0.85])
hold on

% f(1) = scatter3(vpx(1,:),vpy(1,:),vpz(1,:),200,'^');
% set(f(1),'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 0.9 0])
% f(2) = scatter3(vpx(2,:),vpy(2,:),vpz(2,:),200,'^');
% set(f(2),'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 0 0.8])
% f(3) = scatter3(vpx(3,:),vpy(3,:),vpz(3,:),200,'^');
% set(f(3),'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 0.8 1])
hrange = 0:0.1:n*h;
c = brewermap(length(hrange),'blues');c = flipud(c);
c(end,:) = [0.9 0.9 1];

for i = 1:n
    [nn,~] = size(DATA.(fn{i}).Rcenter);
    for j = 1:nn
        x = DATA.(fn{i}).Rcenter(j,1);
        y = DATA.(fn{i}).Rcenter(j,2);
 
        %Rotate & Reflect coordinates
        o = repmat([0.5;0.5],1,length(x));
        xy = [x;y];
        if rotz ~= 0
            xy = Tz*(xy-o)+o;
        end
        if rotx ~= 0
            xy = Tx*(xy-o)+o;
        end
        if roty ~= 0
            xy = Ty*(xy-o)+o;
        end
        
        r = mean(DATA.(fn{i}).Rradii(j,:));
        lmin = DATA.(fn{i}).Layerminh;
        lmax = DATA.(fn{i}).Layermaxh;
        
        tmp = abs(hrange-lmax);
        [~,idx] = min(tmp);
        
        %plot filled cylinders
        coords = meshgrid(linspace(0,2*pi),linspace(0,2*pi));
        X = r.*cos(coords);X = X+xy(1); %transforms X coord to rcenter(x,~)
        Y = r.*sin(coords);Y = Y+xy(2); %transforms Y coord to rcenter(~,y)
        Z = meshgrid(linspace(lmin,lmax),linspace(lmin,lmax))';
        cyl = surf(X,Y,Z);hold on
        set(cyl,'FaceColor',c(idx,:),'FaceAlpha',1,'EdgeAlpha',0);
        p(1) = patch(X(1,:),Y(1,:),Z(1,:),Z(1,:)); %top of volume 
        p(2) = patch(X(end,:),Y(end,:),Z(end,:),Z(end,:)); %bottom of volume
        set((p),'FaceColor',c(idx,:),'EdgeColor',[0 0 0],'EdgeAlpha',0.5);
        
    end
end
% cb = colorbar('eastoutside');
caxis([0 n])
% set(cb,'position',[0.825 0.1 0.02 0.8],'LineWidth',1.5)
%ylabel(cb,'Slice Number','FontSize',16)
axis equal

% leg = legend([f(3) f(2) f(1)],'\bf\itDay 3','\bf\itDay 2','\bf\itDay 1');
% set(leg,'box','off','position',[0.68 0.75 0.025 0.025],'FontSize',18)
set(gca,'LineWidth',1.5,...
    'FontName','Arial',...
    'FontSize',18,...
    'XLim',[0 1],'YLim',[0 1],'ZLim',[0 1],...
    'XTick',0:0.2:1,...
    'YTick',0:0.2:1)
% zlabel('Height Above Bed (m)','FontSize',18)
% xlabel('Along-shore Distance (m)','FontSize',18)
% ylabel('Cross-shore Distance (m)','FontSize',18)
grid on
hold off
% axis([0 1 0 1 0 1])
set(gca,'position',[0.1 0.1 0.7 0.8])
view(-140,25)
%  view(0,0)
export_fig([fdir 'Q4_Volume_oblique_m6'],'-pdf','-nocrop','-M6')
% export_fig([fdir 'Q2B_instLocVols_oblique'],'-png','-native')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following code was designed to calculate volumes of individual
%pneumatophores. In many datasets, the ID field does not exist for all
%layers, making this method defunct. I'm keeping this code because it could
%be useful for estimating the volume of particular elements in the future.

% centers = double(DATA.layer1.ID); %centers of the radii of the pneumatophores at slice 1
% %we need to find the 'height' of the pneumatophores. Do this by extracting
% %the index where the intersection of the centers of the first slice no
% %longer intersect with the centers of slice(j).
% m = length(centers)-1;
% idx = zeros(m,n);
% idx(:,1) = ones(1,m);
% Z = zeros(1,m);
% for i = 1:m
%     for j = 2:n
%         idx(i,j) = any(intersect(centers(i),DATA.(fn{j}).ID));
%         z = find(idx(i,:)==1,1,'last');
%         %z is the last point where the intersect of centers and slice(j) is
%         %not an empty matrix (i.e. where the pneumatophore ends).
%         if isempty(z)
%             Z(i) = n;
%         else
%             Z(i) = z;
%         end
%     end
% end
% %plot individual pneumatophores
% for i = 1:m
%     for j = 1:Z(i)
%         id = find(intersect(centers(i),DATA.(fn{j}).ID) == DATA.(fn{j}).ID);
%         x = DATA.(fn{j}).Rcenter(id,1);
%         y = DATA.(fn{j}).Rcenter(id,2);
%         r = mean(DATA.(fn{j}).Rradii(id,:));
%         ang=0:0.01:2*pi;
%         xp=r*cos(ang);
%         yp=r*sin(ang);
%         Hm = repmat(H(j),1,length(xp));
%         p = plot3(x+xp,y+yp,Hm);
%         hold on
%     end
% end
% 

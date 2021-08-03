%Plot averaged eddy viscosity for the HTA and VTA (2x3 panel subplots).
clear
load('d:\projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
load('D:\prokects\Mekong_W2015\DataAnalysis\Paper2\U_z_gradient_v3.mat')
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600  800]);
set(gcf,'color','w','paperpositionmode','auto')
s = [5 3 1]; %HTA
p = [6 4 2]; %VTA
sp1 = zeros(3,1); %HTA
sp2 = zeros(3,1); %VTA
dn = fieldnames(RS);
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07 0.416 0.806];
hc = [0.58 0.58 0.58 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:4
    fn = fieldnames(RS.(dn{i}));
    for j = 1:3
        uw = RS.(dn{i}).(fn{j}).uw./100;
        time1 = RS.(dn{i}).(fn{j}).time;
        time2 = slope.(dn{i}).(fn{j}).time';
        %average dudz to the same size as uw
        
        dudz = zeros(length(time1),35);
        for kk = 1:length(time1)-1
            td = find(time2 >= time1(kk) & time2 <= time1(kk+1));
            dudz(kk,:) = nanmean(slope.(dn{i}).(fn{j}).dudz(:,td),2);
        end
        dudz(end,:) = slope.(dn{i}).(fn{j}).dudz(:,end);
        
        %eddy viscosity is Re stress/velocity gradient
        nue = uw./dudz;
        %         nue = (uw)./(dudz);
%         nue(isinf(nue)) = NaN;
        nue = nanmean(nue);
        %         for k = 1:length(time1)
        %             mad = nanmedian(abs(nue(k,:)-nanmedian(nue(k,:)))); %remove outliers
        %             id = find(nue(k,:) > mad | nue(k,:) < -mad);
        %             if ~isempty(id)
        %                 nue(k,id) = NaN;
        %             end
        %         end
        %          nue = nanmean(nue);
        if i ~= 1 && i ~= 4
            nue = fastsmooth(nue,5,1,1);
        end
        z = hab(i,j)-0.04-linspace(0,0.03,35);
        zhc = z./hc(i);
        markx = nue(1:5:end);
        marky = zhc(1:5:end);
        if i < 4
            sp1(i) = subplot(3,2,s(i));
            plot(nue,zhc,line{j},...
                'color',cl(j,:),...
                'linewidth',1.5), hold on
            plot(markx,marky,symb{j},...
                'linewidth',1.5,...
                'markerfacecolor',cl(j,:),...
                'markeredgecolor','k',...
                'markersize',6), hold on
            grid on
        else
            sp2(j) = subplot(3,2,p(j));
            plot(nue,zhc,'-',...
                'color','k',...
                'linewidth',1.5), hold on
            plot(markx,marky,'^',...
                'linewidth',1.5,...
                'markerfacecolor','w',...
                'markeredgecolor','k',...
                'markersize',6)
            grid on
        end
    end
end
break
set(sp1,'xlim',[-0.2 0.2],...
    'xtick',-0.2:0.1:0.2,'gridlinestyle',':')
set(sp2,'xlim',[-0.1 0.1],...
    'xtick',-0.1:0.05:0.1,'gridlinestyle',':')
% set(sp1(3),'xticklabel',[])
% set([sp2(2) sp2(3)],'xticklabel',[])
%HTA
set(sp1(1),'ylim',[0 0.04],...
    'ytick',0:0.01:0.04,...
    'xlim',[-0.5 0.5],...
    'xtick',-0.5:0.25:0.5,...
    'position',[0.13 0.08 0.35 0.25])
set(sp1(2),'ylim',[0.29 0.35],...
    'ytick',0.29:0.02:0.35,...
    'position',[0.13 0.39 0.35 0.25])
set(sp1(3),'ylim',[0.82 0.88],...
    'ytick',0.82:0.02:0.88,...
    'position',[0.13 0.7 0.35 0.25])
title(sp1(3),'HTA')
ylabel(sp1(2),'z/h_c')
xl = xlabel(sp1(1),'$\overline{\nu_e} \quad (m^2s^{-1})$');
set(xl,'Interpreter','latex','fontname','arial')
%VTA
set(sp2(1),'ylim',[0.01 0.05],...
    'ytick',0.01:0.01:0.05,...
    'xlim',[0 0.04],...
    'xtick',0:0.01:0.04,...
    'position',[0.6 0.08 0.35 0.25])
set(sp2(2),'ylim',[0.58 0.625],...
    'ytick',0.58:0.015:0.63,...
    'position',[0.6 0.39 0.35 0.25])
set(sp2(3),'ylim',[1.22 1.28],...
    'ytick',1.22:0.02:1.28,...
    'xlim',[-0.2 0.2],...
    'xtick',-0.2:0.1:0.2,...
    'position',[0.6 0.7 0.35 0.25])
title(sp2(3),'VTA')
xl = xlabel(sp2(1),'$\overline{\nu_e} \quad (m^2s^{-1})$');
set(xl,'Interpreter','latex','fontname','arial')
prettyfigures('text',13,'labels',14,'box',1,'tlength',[0.025 0.025])
% export_fig([savefigdir 'EddyVisProf_HTA_VTA_v5'],'-pdf')

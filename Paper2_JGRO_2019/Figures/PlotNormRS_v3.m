%Try plotting the u* ratios by canopy height
clear, close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\NormFrictionU_v3.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 600 500]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'o','d','p'};symb = repmat(symb,3,1);
sy = repmat({'^'},1,3);symb = [symb;sy];
fn = fieldnames(mydata);
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0.42 0.81];
zmean = NaN(4,3);
usmean = NaN(4,3);
p = NaN(4,3);
cmap = csvread('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Design\Qualitative8_4.csv');
cmap = flipud(cmap)./255;
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).uid));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).usnorm);
        xstd = nanstd(mydata.(fn{i}).(dfn{ii}).usnorm);
        disp([fn{i} ' ' dfn{ii} ' drag coef: ' num2str(xs)])
        cb(i,ii) = xs;
        if xs-xstd < 0
            xsd = xs-abs(xstd-xs);
        else
            xsd = [xstd xstd];
        end
        zhc = vh(id)/hc(i);
        eb = ploterr(xs,zhc,xsd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        usmean(i,ii) = nanmean(mydata.(fn{i}).(dfn{ii}).usnorm);
    end
end
break
%plot mean values
p2(1) = plot(usmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(usmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(usmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(usmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.');
plot(linspace(1E-3,1,10),ones(1,10)*1,...
    '-k','linewidth',1.5)
plot(linspace(1E-3,1,10),ones(1,10)*0.33,...
    '--k','linewidth',1.5)
plot(linspace(1E-3,1,10),ones(1,10)*0.23,...
    ':k','linewidth',1.5)
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        m = length(mydata.(fn{i}).(dfn{ii}).usnorm);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).uid));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).usnorm);
        xstd = nanstd(mydata.(fn{i}).(dfn{ii}).usnorm);
        zhc = vh(id)/hc(i);
        p(i,ii) = plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',8,...
            'markerfacecolor',cl,...
            'linewidth',1);
    end
end
% uwghis = [0.0004392783398308245, 0.8103130755064427
% 0.0002213930943212733, 2.1362799263351704
% 0.0003925601715656993, 3.977900552486183
% 0.000001195857549450563, 5.598526703499083
% -0.00034647979399360374, 6.850828729281766
% -0.0016034855261374548, 8.324125230202576
% -0.0033363628391252648, 9.502762430939228
% -0.007536772619645544, 10.681399631675873
% -0.012819274034743644, 11.712707182320443
% -0.015505807881498485, 14.06998158379374
% -0.01801543453477155, 12.96500920810314
% -0.01360240128195929, 15.32228360957643
% -0.013950316105012236, 16.79558011049724
% -0.012696180430987056, 17.974217311233886
% -0.011788365103282229, 19.15285451197054
% -0.01027472834102668, 20.552486187845304
% -0.010752114674766604, 21.65745856353592
% -0.010320968166272027, 23.27808471454881
% -0.00746518061435188, 24.5303867403315
% -0.00521537394465571, 25.70902394106815
% -0.004827198584104655, 27.034990791896874
% -0.004438943499716977, 28.287292817679567
% -0.0017996858880836748, 29.6132596685083
% -0.0029265025950108783, 30.791896869244944
% -0.001066465762598355, 32.11786372007367
% 0.00009957507195076931, 34.69613259668509
% 0.0010056364752497396, 37.495395948434634
% 0.0018694442451348617, 39.33701657458565];
% ughis = [0.8441558441558445, 0.5904059040590468
% 0.7359307359307348, 2.2140221402214095
% 0.7792207792207773, 3.9114391143911504
% 0.7792207792207773, 5.6826568265682695
% 0.8008658008657985, 6.863468634686356
% 0.887445887445887, 8.339483394833948
% 0.9956709956709933, 9.594095940959413
% 1.1471861471861455, 10.7749077490775
% 1.5151515151515138, 11.808118081180815
% 1.883116883116882, 12.915129151291524
% 2.251082251082247, 14.022140221402218
% 2.5757575757575726, 15.498154981549824
% 2.9004329004328984, 16.75276752767528
% 3.095238095238093, 18.00738007380074
% 3.2683982683982666, 19.18819188191882
% 3.506493506493504, 20.590405904059047
% 3.658008658008656, 21.771217712177126
% 3.874458874458872, 23.173431734317347
% 4.004329004329, 24.428044280442805
% 4.134199134199134, 25.68265682656827
% 4.307359307359304, 27.23247232472325
% 4.3722943722943715, 28.48708487084871
% 4.610389610389609, 29.66789667896679
% 4.567099567099566, 30.996309963099634
% 4.675324675324672, 32.103321033210335
% 4.696969696969694, 34.90774907749078
% 4.588744588744584, 37.63837638376384
% 4.480519480519474, 39.33579335793358];
% uwghis = uwghis(:,1).*(ughis(:,1).^2);
% uwghis = uwghis*1E-4; %convert to m2/s2
% ughis = ughis.*0.01; %convert to m/s
% unorm = sqrt(uwghis(:,1).^2)./(ughis(:,1).^2);
% zhc = (ughis(:,2))./0.139;
% plot(unorm,zhc,'ok','linewidth',1.5)
% hold off;

set(gca,'ylim',[0 1.4],...
    'ytick',0:0.2:1.4,...
    'xscale','log',...
    'xlim',[3E-2 1],...
    'xtick',[1E-2 1E-1 1])
tl = xlabel('$\widehat{u_*} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
ylabel('z/h_{max} [-]')
% leg = legend([p(1,:) p(4,:) p2],{'x = -10 cm';'x = 10 cm';'x = 20 cm';...
%     'z/h_c = 0.04';'z/h_c = 0.60';'z/h_c = 1.25';...
%     'x = -10 cm';'x = 10 cm';'x = 20 cm';'VTA'});
% set(leg,'position',[0.81 0.36 0.1 0.42])
prettyfigures('text',13,'labels',14,'box',1)
savefigdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
export_fig([savefigdir 'NormRS_HTA_VTA_v8'],'-pdf')

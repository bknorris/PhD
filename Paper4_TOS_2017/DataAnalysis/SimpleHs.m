%Generate Hs plots for Julia, F2F2

clear
data = struct(); %preallocate data structure
start = datenum(2015,03,05,23,50,00);stop = datenum(2015,03,06,05,49,00);

%%%Data Pre-Processing%%%
%process ADCP data first; everything will be rotated to principal wave axis
inst{1} = 'D:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\HR3_7March2015.mat';
inst{2} = 'D:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\AD5116_9March2015.mat';
inst{3} = 'D:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\AD5117_9March2015.mat';

%%%Spectra settings
fn = fieldnames(data);
zp = [0.072 0.075 0.076];
zuv = zp;
var = {'m';'fr';'fo'};
lf = 1.2;
hf = 0.05;
win = 180; %seconds (5 minutes)
step = 10; %seconds
for i = 1:3
    disp(['Loading ' inst{i}])
    load(inst{i})
    time1 = aqdp.datenum;
    vid = find(time1 >= start & time1 <= stop);
    U = nanmean(aqdp.u(vid),2);V = nanmean(aqdp.v(vid),2);P = nanmean(aqdp.pressure(vid),2);
    x = (sin(aqdp.metadata.lat/57.29578))^2;
    fs = 8;
    avt = fs*step;
    nwin = fs*win;
    swin = fs*10; %30 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
    P = cmgbridge(P,100,100,1000);
    vname = var{i};
    
    %save variables to structure
    data.(vname).p = P;
    data.(vname).u = U;
    data.(vname).v = V;
    data.(vname).time1 = time1(vid);
    time = time1(vid);
    
    %%%Data Analysis
    nsamp = length(vid);
    ind = [1 avt:avt:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
        end
        %spectra of ADVs
        u = U(idx);
        v = V(idx);
        p = P(idx)+zp(i);
        time2 = time(ind(ii));
        
        g = zeros(length(p),1);h = zeros(length(p),1);
        for j = 1:length(p)
            g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
            h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
        end
        h = mean(h);
        g = mean(g);
        p = detrend(p);
        u = detrend(u);
        v = detrend(v);
        %sig wave height should be calculated with non-directional
        %spectra, i.e. from the pressure t-s
        [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
        
        %wave parameters
        lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        kh = k*h;
        kz = k*zuv(i);
        attn = cosh(kz)./cosh(kh);
        attn(attn<0.2) = 0.2;
        Spp = Cpp./(attn.^2);           %surface elevation spectrum
        m0 = sum(Spp(lfc:hfc)*df);
        Hs = 4*sqrt(m0);                %sig. wave height
        
        %save variables to structure
        data.(vname).time2(ii) = time2;
        data.(vname).Hs(ii) = Hs;
    end
end
save('d:\Projects\Mekong_W2015\DataAnalysis\ToS\F2F2_Hs','data','-v7.3')

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');
subplot(211)
pp(1) = plot(data.m.time2,data.m.Hs,'-.','linewidth',1.5,...
    'Color',[0.6 0.6 0.6]);hold on
pp(2) = plot(data.fr.time2,data.fr.Hs,'--','linewidth',1.5,...
    'Color',[0.4 0.4 0.4]);
pp(3) = plot(data.fo.time2,data.fo.Hs,'-','linewidth',1.5,...
    'Color',[0 0 0]);
leg = legend(pp,{'Mudflat','Fringe','Forest'});
set(leg,'position',[0.86 0.85 0.05 0.05])
ylabel('Hs (m)')
set(gca,'xticklabel',[])
subplot(212)
pq(1) = plot(data.m.time2,data.fr.Hs./data.m.Hs,'-r','linewidth',1.5);hold on
pq(2) = plot(data.m.time2,data.fo.Hs./data.m.Hs,'-b','linewidth',1.5);
leg = legend(pq,{'Hs_f_r/Hs_m_u_d';'Hs_f_o/Hs_m_u_d'});
set(leg,'position',[0.86 0.38 0.05 0.05])
datetick('x','HH:MM','keepticks','keeplimits')
ylabel('Hs_f_r_,_f_o/Hs_m_u_d')
xlabel('Time on 05/03/15')
prettyfigures('labels',14,'fweight','bold','fangle','italic','box',1)
export_fig('d:\Projects\Documents\Writing\TOSpaper\Figures\Hs\F2F2_Hs_compare','-png')
%IG Band spectra - F2F2
%Load Aquadopps collocated with Vectrinos, calculate spectra and plot
%results
clear, close all
fdir = 'd:\Mekong_W2015\Data\Aquadopp\F2F2\';
fileList = {'HR3_7March2015.mat';'AD5116_9March2015.mat';'AD5117_9March2015.mat'};
start = datenum(2015,03,05,15,32,10);
stop = datenum(2015,03,05,16,52,30);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0.1 0.1 0.1];
line = {'-';'-';'-'};
sp = zeros(1,3);
for i = 1:length(fileList)
    disp(['Loading ' fdir fileList{i}])
    load([fdir fileList{i}])
    time1 = aqdp.datenum;
    vid = find(time1 >= start & time1 <= stop);
    hp = aqdp.metadata.HAB/1000; %height of pressure sensor
    P = aqdp.pressure(vid);
    fs = str2double(regexprep(aqdp.metadata.samprate,' Hz',''));
    p = cmgbridge(P,100,100,1000);p = detrend(p);
    %%%
    if mod(length(p),2)
        p = p(1:end-1);
    end
    dt = 1/fs;[np,nt]=size(p);
    Dt=np*dt;
    f=[1:(np/2)]/Dt;            % frequency array up to Nyquist
    w = 0.54-0.46*cos(2*pi*[0:np-1]'/(np-1));
    pf=(fft(p.*w));
    pp=abs(pf(2:(np/2+1),:).^2);
    pp=pp*2/np^2/f(1);
    [F, Cpp, dF, Ns, Ne] = logavg(f,pp,100);
    %sea surface scaling
    mp=mean(p);        % vertical scaling
    F2=F*ones(1,nt);
    k=wavek(F2,mp+hp);
    coshkh=cosh(k.*(hp+mp));        %cosh scaling for water depth
    coshkhz=cosh(k.*hp);            %scaling for pressure sensor elevation
    Sp=Cpp.*(coshkh./coshkhz).^2;
    %02/03/2020 compute wave energy flux for Steve
    kh = k*(hp+mp);
    c = sqrt(9.81.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c; Sp = Sp.*Cg*(9.81)*1025;
    
    sp(i) = plot(F(5:end),Sp(5:end),...
        line{i},...
        'color',cl(i,:),...
        'linewidth',1.5);
    if i == 1;
        %CI's
        DOF=2*(Ne-Ns+1);
        v = DOF;
        alfa = 1-0.95;
        ci = zeros(length(DOF),2);
        for j = 1:length(DOF)
            c = chi2inv([1-alfa/2 alfa/2],v(j));
            ci(j,:) = v(j)./c;
        end
        s2est = std(p);
        ff = sum(w.*w)/np;
        cs = (ff*s2est)/(2*(np/2+1)*(1/fs));
        ci = ci.*cs;
        cf = [Sp-ci(:,1) Sp+ci(:,2)];cf(cf<0) = NaN;
        nid = ~isnan(cf(:,1));
        L = mean(10*log10(Sp(nid))-10*log10(cf(nid,1)));
        U = mean(10*log10(cf(nid,2))-10*log10(Sp(nid)));
        er = ploterr(2,10^-4,10^-4*log10(L)*10,10^-4*log10(U)*10,'.','logy','hhxy',0.1);
        set(er,'linewidth',1.5,'color','k')
    end
    clear aqdp
end
set(gca,...
    'xscale','log',...
    'yscale','log')
text(2.5,10^-4,'\bfx 10')
leg = legend(sp,{'Mudflat';'Fringe';'Forest'});
set(leg,'position',[0.8 0.8 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
xlabel('f [Hz]')
ylabel('Energy Flux [J/m^-^2]')
sfdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper3\Figures\';
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[sfdir 'WaveSpectraExmp_v3'],'-dpdf','-r0')
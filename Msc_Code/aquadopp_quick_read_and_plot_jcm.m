function  aquadopp=aquadopp_quick_read_and_plot_jcm(filename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in currents, backscatter and sen file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[aquadopp.vel_b1]= textread([filename '.v1'],'','delimiter',' /:','headerlines',0);
[aquadopp.vel_b2]= textread([filename '.v2'],'','delimiter',' /:','headerlines',0);
[aquadopp.vel_b3]= textread([filename '.v3'],'','delimiter',' /:','headerlines',0);


[aquadopp.beam1]= textread([filename '.a1'],'','delimiter',' /:','headerlines',0);
[aquadopp.beam2]= textread([filename '.a2'],'','delimiter',' /:','headerlines',0);
[aquadopp.beam3]= textread([filename '.a3'],'','delimiter',' /:','headerlines',0);


[month, day, year, hour, minute, second, aquadopp.error_code, aquadopp.status_code, ...
    aquadopp.batt_voltage, aquadopp.soundspeed, aquadopp.heading, aquadopp.pitch, aquadopp.roll,...
        aquadopp.pressure, aquadopp.temperature, junk1, junk2]= textread([filename '.sen'],'','delimiter',' /:','headerlines',0);

aquadopp.yearday = yearday(year,month,day+hour/24+minute/(24*60)+second/(24*60*60));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in hdr file and check settings/errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hdrfile]= textread([filename '.hdr'],'%s ','delimiter','\t','headerlines',0);

checksumerrors= str2num(hdrfile{4}(39:end));
starttime=hdrfile{5};
endtime=hdrfile{6};
dt=str2num(hdrfile{9}(39:40));
ncells=str2num(hdrfile{10}(39:end));
cellsize=str2num(hdrfile{11}(39:41))/100;
measload=str2num(hdrfile{13}(39:41));
blankdist=str2num(hdrfile{15}(39:42));
wavemeas=hdrfile{18}(39:end);
coordsys=hdrfile{27}(39:end);
serialno=str2num(hdrfile{87}(43:end));


if checksumerrors~=0
    disp('There are checksum errors')
end

if wavemeas~='DISABLED'
    disp('Wave measurements are ON')
end

if coordsys~='BEAM'    
    disp('Coordinate system is NOT beam')
end

aquadopp.rangebins= blankdist+cellsize:cellsize:blankdist+ncells*cellsize; %This should now be correct 2014.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Basic plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(211)
    plot(aquadopp.yearday, aquadopp.pressure,'b')
    title(['Aquadopp ' num2str(serialno) ': pressure'])
    ylabel('P (dbar)')
subplot(212)
    plot(aquadopp.yearday, aquadopp.temperature,'b')
    title('temperature')
    ylabel('T (deg C)')
    xlabel('yearday')

figure(2)
m1= subplot(311);
    plot(aquadopp.yearday,aquadopp.heading,'r')
    set(m1,'Ylim', [0 360],'YTick',[0:90:360])
    title(['Aquadopp ' num2str(serialno) ':heading'])
m2= subplot(312);
    plot(aquadopp.yearday,aquadopp.pitch,'r')
    set(m2,'Ylim', [0 360],'YTick',[0:90:360])
    title('pitch')
m3 = subplot(313);
    plot(aquadopp.yearday,aquadopp.roll,'r')
    set(m3,'Ylim', [0 360],'YTick',[0:90:360])
    title('roll')
    xlabel('yearday')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot backscatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
ax(1) = subplot(311);
    uimagesc(aquadopp.yearday,aquadopp.rangebins, aquadopp.beam1')
    hold on
    plot(aquadopp.yearday,aquadopp.pressure,'k')
    title(['Aquadopp ' num2str(serialno) ': backscatter'])
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
ax(2) = subplot(312);
    uimagesc(aquadopp.yearday,aquadopp.rangebins,aquadopp.beam2')
    hold on
    plot(aquadopp.yearday,aquadopp.pressure,'k')
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
ax(3) = subplot(313);
    uimagesc(aquadopp.yearday,aquadopp.rangebins, aquadopp.beam3')
    hold on
    plot(aquadopp.yearday,aquadopp.pressure,'k')
    xlabel('yearday')
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
    set(ax,'YDir','Normal','YLim',[0 3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(coordsys,'BEAM') 
    figure(4)
    ax(1) = subplot(311);
        uimagesc(aquadopp.yearday,aquadopp.rangebins, aquadopp.vel_b1')
        hold on
        plot(aquadopp.yearday,aquadopp.pressure,'k')
        title(['Aquadopp ' num2str(serialno) ': beam velocities'])
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
    ax(2) = subplot(312);
        uimagesc(aquadopp.yearday,aquadopp.rangebins,aquadopp.vel_b2')
        hold on
        plot(aquadopp.yearday,aquadopp.pressure,'k')
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
    ax(3) = subplot(313);
        uimagesc(aquadopp.yearday,aquadopp.rangebins, aquadopp.vel_b3')
        hold on
        plot(aquadopp.yearday,pressure,'k')
        xlabel('yearday')
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
        set(ax,'YDir','Normal','YLim',[0 3])
end
   
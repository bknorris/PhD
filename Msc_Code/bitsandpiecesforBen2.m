clear
close all


adnumber=1%input('Enter Aquadopp number: ');
fnmonth='May'%input('Enter month: ');
dlnumber=3%input('Enter download number: ');

dirname = ['../Raw/HRaquadopps/' fnmonth '/'];  

    
filename=['HR' num2str(adnumber) 'download' num2str(dlnumber)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in currents and sen file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp('read one file')
[month, day, year, hour, minute, second, aquadopp.burstcounter, aquadopp.ensemblecounter,aquadopp.error_code, aquadopp.status_code, ...
    aquadopp.batt_voltage, aquadopp.soundspeed, aquadopp.heading, aquadopp.pitch, aquadopp.roll,...
        aquadopp.pressure, aquadopp.temperature, junk1, junk2]= textread([dirname filename '.sen'],'','delimiter',' /:','headerlines',0);

aquadopp.yearday = yearday(year,month,day+hour/24+minute/(24*60)+second/(24*60*60));

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in hdr file and check settings/errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Interpolate (with spline) met station data to Aquadopp times and subtract air pressure
%  variations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lose junk when instrument was out of the water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %ind 3 and ind4 mark start/end indices
 


aquadopp.pressure=aquadopp.pressure(ind3:ind4);
aquadopp.temperature=aquadopp.temperature(ind3:ind4);
aquadopp.yearday=aquadopp.yearday(ind3:ind4);
aquadopp.heading=aquadopp.heading(ind3:ind4);
aquadopp.pitch=aquadopp.pitch(ind3:ind4);
aquadopp.roll=aquadopp.roll(ind3:ind4);
aquadopp.burstcounter=aquadopp.burstcounter(ind3:ind4);
aquadopp.ensemblecounter=aquadopp.ensemblecounter(ind3:ind4);
vel_b1=vel_b1(ind3:ind4,:);
vel_b2=vel_b2(ind3:ind4,:);
vel_b3=vel_b3(ind3:ind4,:);
pressure_corrected=pressure_corrected(ind3:ind4);

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Split up record into bursts 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %[vel_b1]= textread([dirname filename '.v1'],'','delimiter',' /:','headerlines',0);
%[vel_b2]= textread([dirname filename '.v2'],'','delimiter',' /:','headerlines',0);
%[vel_b3]= textread([dirname filename '.v3'],'','delimiter',' /:','headerlines',0);


    
    numbersections=aquadopp.burstcounter(end);

    for ii=aquadopp.burstcounter(1):numbersections
        clear ind5 aquadopp_burst burst_filename
        ind5=find(aquadopp.burstcounter==ii);
       

        aquadopp_burst.rangebins=aquadopp.rangebins;
        aquadopp_burst.pressure=aquadopp.pressure(ind5(1):ind5(end));
        aquadopp_burst.pressure_corrected=pressure_corrected(ind5(1):ind5(end));
        aquadopp_burst.temperature=aquadopp.temperature(ind5(1):ind5(end));
        aquadopp_burst.yearday=aquadopp.yearday(ind5(1):ind5(end));
        aquadopp_burst.heading=aquadopp.heading(ind5(1):ind5(end));
        aquadopp_burst.pitch=aquadopp.pitch(ind5(1):ind5(end));
        aquadopp_burst.roll=aquadopp.roll(ind5(1):ind5(end));
        aquadopp_burst.status_code=aquadopp.status_code(ind5(1):ind5(end));
        for jj=1:ncells
            
            
            % Note for Ben
            % the depth indices go from jj+2 here as the first two columns
            % are the number of the data point in the burst and something
            % else (I've forgotton what). So you can just get rid of these.
            % Just double check the (size,2) of vel files are ncells+2 in
            % case this has changed!
            
            
            aquadopp_burst.vel_b1(:,jj)=vel_b1(ind5(1):ind5(end),jj+2);
            aquadopp_burst.vel_b2(:,jj)=vel_b2(ind5(1):ind5(end),jj+2);
            aquadopp_burst.vel_b3(:,jj)=vel_b3(ind5(1):ind5(end),jj+2);
        end
        burst_filename=['../Processed/HRaquadopps/HR' num2str(adnumber) '/' fnmonth  ...
            '/HR' num2str(adnumber)  'download' num2str(dlnumber)   '_burst' num2str(ii)];
        
        save (burst_filename,'aquadopp_burst','blankdist','burstinterval','cellsize',...
            'ncells','orientation','serialno' ,'pressure_corrected')

    end


    clear vel_b1 vel_b2 vel_b3


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Repeat for backscatter (splitting into burst files)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     [beam1]= textread([dirname filename '.a1'],'','delimiter',' /:','headerlines',0);
    beam1=beam1(ind3:ind4,:);   %%%THIS IS JUST GETTING RID OF STUFF OUT OF WATER
    [beam2]= textread([dirname filename '.a2'],'','delimiter',' /:','headerlines',0);
    beam2=beam2(ind3:ind4,:);
    [beam3]= textread([dirname filename '.a3'],'','delimiter',' /:','headerlines',0);
    beam3=beam3(ind3:ind4,:);


    for ii=aquadopp.burstcounter(1):numbersections
        clear ind5 aquadopp_burst burst_filename
  burst_filename=['../Processed/HRaquadopps/HR' num2str(adnumber) '/' fnmonth  ...
            '/HR' num2str(adnumber)  'download' num2str(dlnumber)   '_burst' num2str(ii)];
               load(burst_filename)
        ind5=find(aquadopp.burstcounter==ii);

        for jj=1:ncells
            aquadopp_burst.beam1(:,jj)=beam1(ind5(1):ind5(end),jj+2);
            aquadopp_burst.beam2(:,jj)=beam2(ind5(1):ind5(end),jj+2);
            aquadopp_burst.beam3(:,jj)=beam3(ind5(1):ind5(end),jj+2);
        end

        save (burst_filename,'aquadopp_burst','blankdist','burstinterval','cellsize',...
            'ncells','orientation','serialno')

    end

    clear beam1 beam2 beam3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Repeat for correlation (splitting into burst files)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    [cor1]= textread([dirname filename '.c1'],'','delimiter',' /:','headerlines',0);
    cor1=cor1(ind3:ind4,:);
    [cor2]= textread([dirname filename '.c2'],'','delimiter',' /:','headerlines',0);
    cor2=cor2(ind3:ind4,:);
    [cor3]= textread([dirname filename '.c3'],'','delimiter',' /:','headerlines',0);
    cor3=cor3(ind3:ind4,:);

    for ii=aquadopp.burstcounter(1):numbersections
        clear ind5 aquadopp_burst burst_filename
        
        burst_filename=['../Processed/HRaquadopps/HR' num2str(adnumber) '/' fnmonth  ...
            '/HR' num2str(adnumber)  'download' num2str(dlnumber)   '_burst' num2str(ii)];
          load(burst_filename)

        ind5=find(aquadopp.burstcounter==ii);

        for jj=1:ncells
            aquadopp_burst.cor1(:,jj)=cor1(ind5(1):ind5(end),jj+2);
            aquadopp_burst.cor2(:,jj)=cor2(ind5(1):ind5(end),jj+2);
            aquadopp_burst.cor3(:,jj)=cor3(ind5(1):ind5(end),jj+2);
        end


        save (burst_filename,'aquadopp_burst','blankdist','burstinterval','cellsize',...
            'ncells','orientation','serialno')

    end
    






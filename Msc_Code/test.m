clear
close all
doplot=1
 
    close all
    burst_filename='/Users/julia/JULIA/Research/FieldExperiments/MekongDelta/SeptemberFieldExperiment/mfiles/HR4_26092014_psub.mat';
    load(burst_filename)

    %if mean(aquadopp_burst.pressure)<=0.3 && (mean(aquadopp_burst.beam1(:,4))>60) && (mean(aquadopp_burst.cor1(:,4))>40)

    aquadopp_burst=aqdp;
    
    if doplot==1
    figure(1)
    a(1) = subplot(311);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.cor1')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
%    title(['aquadopp burst ' num2str(adnumber) ': correlations beam1'])
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    a(2) = subplot(312);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins,aquadopp_burst.cor2')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    c = colorbar;
    axis(c);
    title('beam 2')
    ylabel('range (m)') 
    a(3) = subplot(313);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.cor3')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    title('beam 3')

    xlabel('yearday 2008')
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    set(a,'YDir','Normal','YLim',[0 0.8])

    figure(2)
    a(1) = subplot(311);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.beam1')
    hold on
    plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
   % title(['aquadopp burst ' num2str(adnumber) ': backscatter beam1'])
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    a(2) = subplot(312);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins,aquadopp_burst.beam2')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    c = colorbar;
    axis(c);
    title('beam 2')
    ylabel('range (m)')
    a(3) = subplot(313);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.beam3')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    title('beam 3')

    xlabel('yearday 2008')
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    set(a,'YDir','Normal','YLim',[0 0.8])

    figure(3)
    a(1) = subplot(311);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.vel_b1')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
 %  title(['aquadopp burst ' ': VEL beam1'])
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    a(2) = subplot(312);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins,aquadopp_burst.vel_b2')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    c = colorbar;
    axis(c);
    title('beam 2')
    ylabel('range (m)')
    a(3) = subplot(313);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.vel_b3')
    hold on
    %plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    title('beam 3')

    xlabel('yearday 2008')
    c = colorbar;
    axis(c);
    ylabel('range (m)')
    set(a,'YDir','Normal','YLim',[0 0.8])

    end
    
    ccrit=70;
    
     %%%%%% Add in check on percentage of good correlations
    number_goodbeam1=find(aquadopp_burst.cor1>=ccrit);
    number_goodbeam2=find(aquadopp_burst.cor2>=ccrit);
    number_goodbeam3=find(aquadopp_burst.cor3>=ccrit);
    pc_good=sum([length(number_goodbeam1),length(number_goodbeam2),length(number_goodbeam3)])/(aqdp.nCells*4096*3);
    
    
    vwrap=(max(aquadopp_burst.vel_b1(:))-min(aquadopp_burst.vel_b1(:)))*0.5;
    aquadopp_burst=interpolate_bad_correlations(aquadopp_burst,ccrit); 
    
    %%%%%FLIP HERE FOR UPWARD LOOKING
    
    aquadopp_burst.vel_b1=fliplr(aquadopp_burst.vel_b1);
    aquadopp_burst.vel_b2=fliplr(aquadopp_burst.vel_b2);    
    aquadopp_burst.vel_b3=fliplr(aquadopp_burst.vel_b3);
  

%%% UNCOMMENT TO FLIP BACK FOR UPWARD LOOKING
    w1=unwrap_w_prof1(aquadopp_burst.vel_b1,vwrap);
   w1=fliplr(w1);
    w2=unwrap_w_prof1(aquadopp_burst.vel_b2,vwrap);
   w2=fliplr(w2);
    w3=unwrap_w_prof1(aquadopp_burst.vel_b3,vwrap);
   w3=fliplr(w3);
    
    
    aquadopp_burst.vel_b1_uw=w1; 
    aquadopp_burst.vel_b2_uw=w2; 
    aquadopp_burst.vel_b3_uw=w3;
    
    if doplot==1
    figure(5)
    a(1) = subplot(311);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.vel_b1_uw')
    hold on
    plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
   title(['aquadopp burst '  ': VELS UNWRAPPED beam1'])
    c = colorbar;
    caxis([-vwrap vwrap])
    axis(c);
    ylabel('range (m)')
    a(2) = subplot(312);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins,aquadopp_burst.vel_b2_uw')
    hold on
    %lot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    c = colorbar;
    caxis([-vwrap vwrap])
    axis(c);
    title('beam 2')
    ylabel('range (m)')
    a(3) = subplot(313);
    imagesc(aquadopp_burst.yearday,aquadopp_burst.rangebins, aquadopp_burst.vel_b3_uw')
    hold on
    plot(aquadopp_burst.yearday,aquadopp_burst.pressure,'k')
    title('beam 3')

    xlabel('yearday 2008')
    c = colorbar;
    caxis([-vwrap vwrap])
    axis(c);
    ylabel('range (m)')
    set(a,'YDir','Normal','YLim',[0 0.8])

    end
    tilefigs([3 5])
    pause;
    
    unwrapOK='n';%input('Is unwrap satisfactory (y or n)?: ','s');
%     
   if unwrapOK=='y'
         % [aquadopp_burst]=rotate_unwrapped_pipe(aquadopp_burst,'u',heading);
        
       
    
   
    burst_filename=(['../../Processed/HRaquadopps/HR' num2str(adnumber) '/' month '/HR'...
        num2str(adnumber) 'download' num2str(dlnumber) '_burst' num2str(bn) '_unwrapped']);
     save(burst_filename,'aquadopp_burst','ncells','serialno','adnumber','orientation',...                 
                'blankdist','burstinterval','cellsize','vwrap')
      end
   
  %  end

 

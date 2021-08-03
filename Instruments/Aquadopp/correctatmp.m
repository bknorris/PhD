function finalpres = correctatmp(met_time,met_pres,inst_datetime,inst_pres,inst_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax: 
%   
%      p = correctatmp(met_time,met_pres,aqdp_datetime,aqdp_pres,aqdp_temp)
%
%  Inputs: met_time: time from weather station in datenum format
%          met_pressure: pressure from a weather station in dBar
%          inst_datetime: datetime in datenum format
%          inst_pres: pressure from instrument in dBar
%          inst_temp: temperature data from instrument
%          Note: inputs must be vectors and be clipped to the same
%          start/end times.
%
%  Interpolates (with spline) met station data to instrument times and subtracts
%  atmospheric pressure variations. Corrects Pressure signal due to temperature 
%  variance and P-T relationship.
%
%  Developed by Dr. Julia C Mullarney, University of Waikato, New Zealand 
%  c. 2013
%
%  Additions made by Benjamin K Norris, 2014
%  Edits: Cleaned up and organized original code fragments
%         removed extraneous bits of code
%         Added routine to average the user selected 'out of water'
%         pressure data and correct pressure data with this offset
%         Added routine to prompt user for data processing progress plots
%         Removed bit of code that zeros all the out of water pressure data
%         in an attempt to correct offset seen when plotting pressure and
%         backscatter data together.
%         Added check for statistical fit of P-T relationship, set limit to
%         85%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    help(mfilename)
    return
end

disp('Correcting for barometric pressure')
pause(1)
yd = met_time;
airpressure = spline(yd,met_pres,inst_datetime);
newpressure = inst_pres+mean(airpressure)-airpressure;
disp('Pressure signal corrected')
pause(1)
rsq = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P-T varation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while rsq < 0.5 %set a limit for the polynomial fit r-squared value
    
    disp('USER select out of water start and end points')
    
    figure(6)
    plot(inst_datetime,newpressure)
    datetick('x','keepticks')
    ylabel('dBar')
    title('Barometric Corrected Pressure Signal')
    hold on
    nsec = 2;  %Number of selections: start and finish of Out of Water period
    clear q
    q = ginput(2);
    for ii = 1:nsec
        [~, startout(ii)] = min(abs(q(1,1)-inst_datetime));
        [~, endout(ii)] = min(abs(q(2,1)-inst_datetime));
    end
    for ii = 1:nsec   %Combines all data into single variable to plot it.
        if ii == 1
            p_all = newpressure(startout(ii):endout(ii))';
            temp_all = inst_temp(startout(ii):endout(ii))';
        else
            p_all = [p_all newpressure(startout(ii):endout(ii))'];
            temp_all = [temp_all inst_temp(startout(ii):endout(ii))'];
        end
    end
    p = polyfit(temp_all,p_all,1);
    
    %calculate r-squared value of the fit
    stat = fitlm(temp_all,p_all);
    rsq = stat.Rsquared.Ordinary;
    if rsq < 0.5
        disp(['The r-squared value of temperature and pressure is ' num2str(rsq)]);
        disp('WARNING: This value is below the recommended statistical fit of 50%')
    else
        disp(['The r-squared value of temperature and pressure is ' num2str(rsq)]);
    end
end

newpressure2 = newpressure-p(1)*inst_temp; %Removes T-variation
disp('Pressure-temperature dependent variance corrected')
pause(1)
disp('Applying offset')
pause(1)
close(6)
%use the user selected q coordinates to designate the offset values
ind = (inst_datetime >= q(1,1) & inst_datetime <= q(2,1));
pres = newpressure2(ind);pave = sum(pres)/numel(pres); %take the average of out of water values
finalpres = newpressure2 - pave; %apply the offset
finalpres(finalpres < 0 )=0; %negate any negative pressures
        
prompt = 'Plot figure [y/n]? ';
result = input(prompt,'s');
if strcmp(result,'y') || strcmp(result,'yes');
    figure(7)
    plot(inst_datetime,inst_pres,'k')
    hold on
    plot(inst_datetime,newpressure,'r')
    hold on
    plot(inst_datetime,finalpres,'m')
    datetick('x','keepticks')
    title('Pressure Adjustment')
    legend('Uncorrected','Pres Corrected','Temp & Pres Corrected')
    ylabel('dBar')
    pause(4)
    close(7)
end

if strcmp(result,'n') || strcmp(result,'no');
end


end
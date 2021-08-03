function day=date2day(month,dayte,hour,mn,sc,yr);

% DATE2DAY(MONTH,DAYTE,HR,MN,SC,YR) returns time in yearday corresponding 
% to specified time in month, date, hour etc 
%
%	MONTH=month of year.
%	DAYTE=day of month.
%	HOUR=hour of day.
%	MN=minute of hour.
%	SC=second of minute.
%   YR=year (default=non-leap-year).  
%
% Alternatively, DATE2DAY('MMDDHHNN',yr) returns time in yearday for a string,
% with MM=Month, DD=day, HH=hour (24 hr clock), NN=minute.  yr=year, for
% dealing with leap years, default=non-leap year, in not specified.
%
%    Examples:
%	day=date2day(9,20,15,23,2.4);
%	  ...sets the variable 'day' equal to 263.641.  
%	OR can omit latter parameters such as hour, min & sc
%	day=date2day(9,20);
%	  ...sets the variable 'day' equal to 263
%   OR can specify string,
%   day=date2day('09201523')
%     ...sets the variable 'day' equal to 263.64097...
%
%    See also day2date.  

if isstr(month)
    if exist('dayte')
        yr=dayte;
    end
    sc=0;
    mn=str2num(month(:,7:8));
    hour=str2num(month(:,5:6));
    dayte=str2num(month(:,3:4));
    month=str2num(month(:,1:2));
end
if ~exist('sc')
  sc=0;
end
if ~exist('mn')
  mn=0;
end
if ~exist('hour')
  hour=0;
end
if ~exist('dayte')
  dayte=0;
end
if ~exist('yr')
    yr=1997;
end
day=datenum(yr,month,dayte,hour,mn,sc)-datenum(yr-1,12,31);

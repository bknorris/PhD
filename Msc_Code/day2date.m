function date=day2date(day,yr);

% DAY2DATESHORT(DAY) converts time in yearday, specified by DAY, to 
% mm/dd/hh/nn/ss.s format, where 
% mm=month, dd=day, hh=hour, nn=minute, ss.s=second.  
% Pre:  DAY=time in yearday;
%       YR=year (important only for determining whether leap year.  If not
%       specified, defaults to non-leap-year).
% Post: date=string, mm/dd/hh/nn/ss.s
%
%    Example:
%	date=day2date(263.641)
%	  ...sets the variable 'date' equal to '09/20/15/23/02.4'.
%
%    See also date2day.

if ~exist('yr')
    yr=1997;
end
mmdd=datestr(day+datenum(yr-1,12,31),6);
time_of_day=day-floor(day);
hh=floor(time_of_day*24);
thh=hh/24;
nn=floor(60*24*(time_of_day-thh));
tnn=thh+nn/(24*60);
ss=floor(60*60*24*(time_of_day-tnn));
tss=tnn+ss/(24*60*60);
d=round(10*60*60*24*(time_of_day-tss));
if d==10
  ss=ss+1;
  if ss==60
    nn=nn+1;
    ss=0;
    if nn==60
      hh=hh+1;
      nn=0;
      if hh==24
        mmdd(4:5)=num2str(str2num(mmdd(4:5))+1);
        hh=0;
      end
    end
  end
  d=0;
end
if hh<10
  hh=['0' num2str(hh)];
else
  hh=num2str(hh);
end
if nn<10
  nn=['0' num2str(nn)];
else
  nn=num2str(nn);
end
if ss<10
  ss=['0' num2str(ss)];
else
  ss=num2str(ss);
end
d=num2str(d);
date=[mmdd '/' hh '/' nn '/' ss '.' d];

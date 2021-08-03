%function: input = (year,month,day); RPM - Mar/04
function [yearday] = yearday(year,month,day)

date = datenum(year,month,day);
jan1 = datenum(year,1,1) -1;

yearday = date - jan1;
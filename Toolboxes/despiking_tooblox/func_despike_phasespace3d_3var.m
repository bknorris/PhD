function [vxc,vyc,vzc,ip] = func_despike_phasespace3d_3var(vx,vy,vz,i_opt)
%======================================================================
%
% Version 1.00
%
%======================================================================
%
% This function excludes spike noise from Acoustic Doppler Velocimetry
% (ADV) data using phase-space method by Goring and Nikora (2002).
% This function calls func_despike_phasespace.m.
%
% Input
%   vx    : input x-direction velocity component
%   vy    : input y-direction velocity component
%   vz    : input z-direction velocity component
%   i_opt : = 1 or noinput ; return spike noise as NaN
%           = 2            ; remove spike noise and variable becomes 
%                            shorter than input
%           = 3            ; interpolate NaN using cubic polynomial
%           = negative     ; plot results for check
%
% Output
%   vxc   : corrected x-direction velocity component
%   vyc   : corrected x-direction velocity component
%   vzc   : corrected x-direction velocity component
%   ip    : excluded array element number of vx, vy vz
%
% Example: 
%   [V,SNR,COR,time] = winadv_read('test.txt');
%    vx      = V(:,1);
%    vy      = V(:,2);
%    vz      = V(:,3);
%    [vxc, vyc, vzc, ip] = func_despike_phasespace3d_3var( vx, vy, vz, 2 );
%   
%
%======================================================================
% Terms:
%
%       Distributed under the terms of the terms of the BSD License
%
% Copyright:
%
%       Nobuhito Mori
%           Disaster Prevention Research Institue
%           Kyoto University
%           mori@oceanwave.jp
%
%========================================================================
%
% Update:
%       1.01    2009/06/09 modified func_despike_phasespace3d.m
%       1.00    2004/09/01 Nobuhito Mori
%
%========================================================================

if nargin<=3
  i_opt = 1;
end
disp('Despiking...')
%
% --- despike for each velocity component
%
disp('Beam 1')
[~, ipx] = func_despike_phasespace3d(vx);
disp('Beam 2')
[~, ipy] = func_despike_phasespace3d(vy);
disp('Beam 3')
[~, ipz] = func_despike_phasespace3d(vz);

%
% --- surmmarize excluding data 
%

uxc = vx;
uyc = vy;
uzc = vz;

uxc(ipx) = NaN;
uxc(ipy) = NaN;
uxc(ipz) = NaN;

uyc(ipx) = NaN;
uyc(ipy) = NaN;
uyc(ipz) = NaN;

uzc(ipx) = NaN;
uzc(ipy) = NaN;
uzc(ipz) = NaN;

ip = find(isnan(uxc));

%
% --- interpolation or shorten NaN data
%

% remove NaN from data
if abs(i_opt) >= 2
  inan = find(~isnan(uxc));
  vxc = uxc(inan);
  vyc = uyc(inan);
  vzc = uzc(inan);

  % interpolate NaN data
  if abs(i_opt) == 3
    x   = find(~isnan(uxc));
    xi  = 1:max(size(vx));
    vxc = interp1(x,vxc,xi,'pchip')';
    vyc = interp1(x,vyc,xi,'pchip')';
    vzc = interp1(x,vzc,xi,'pchip')';
  end
else
  vxc = uxc;
  vyc = uyc;
  vzc = uzc;
end

%
% --- for check
%

if sign(i_opt) == -1


figure('PaperOrientation','portrait',...
    'position',[400 400   1000   550]);
ax1 = subplot(3,1,1);
plot(vx,'b-');
hold on
plot(vxc,'r-');
title('Despiking Routine comparison')
legend({'Original';'New'},-1)
ylabel('m/s')
xlabel('n Samples')
hold off

ax2 = subplot(3,1,2);
plot(vy,'b-');
hold on
plot(vyc,'r-');
ylabel('m/s')
xlabel('n Samples')
hold off

ax3 = subplot(3,1,3);
plot(vz,'b-');
hold on
plot(vzc,'r-');
ylabel('m/s')
xlabel('n Samples')
hold off
ax=get(ax1,'Position');
ax(2) = ax(2)-0.3;
set(ax2,'Position',ax);
ax(2) = ax(2)-0.3;
set(ax3,'Position',ax);
linkaxes([ax1 ax2 ax3],'x')
end

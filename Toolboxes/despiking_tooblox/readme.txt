%======================================================================
Outline
%======================================================================
This subroutine excludes spike noise from Acoustic Doppler Velocimetry (ADV) data using phasce-space method by Modified Goring and Nikora (2002) method by Nobuhito Mori (2005).

%======================================================================
Programs
%======================================================================

(1) Despiking a single variable 'fi'
- func_despike_phasespace3d.m
	This program calls
		func_excludeoutlier_ellipsoid3d
		func_excludeoutlier_cor

(1) Despiking three variables 'xi', 'yi' and 'zi' (such as fliud velocity components) 
- func_despike_phasespace3d_3var
	This program calls
		func_despike_phasespace3d_3var.m
		func_excludeoutlier_ellipsoid3d
		func_excludeoutlier_cor

%======================================================================
Usage
%======================================================================

(1) func_despike_phasespace3d.m
	[vxc, vyc, vzc, ip] = func_despike_phasespace3d_3var( vx, vy, vz, i_opt );
	Input
		fi     	input data with dimension (n,1)
		i_plot 	=9 	plot results (optional)
			=other	noplot
		i_opt :	= 0 or not specified  ; return spike noise as NaN
			= 1            ; remove spike noise and variable becomes shorter than input length
			= 2            ; interpolate NaN using cubic polynomial
	Output
		fo	output (filterd) data
		ip	excluded array element number in xi and yi

(2) func_despike_phasespace3d_3var
	[fo, ip] = func_despike_phasespace3d( fi, i_plot, i_opt );
	Input
		vx,vy,vz,	input data with dimension (n,1)
		i_opt :	= 0 or not specified  ; return spike noise as NaN
			= 1            ; remove spike noise and variable becomes shorter than input length
			= 2            ; interpolate NaN using cubic polynomial
	Output
		vxc,vyc,vzc
			output (filterd) data
		ip	excluded array element number in xi and yi

%======================================================================
Example
%======================================================================

(1)  func_despike_phasespace3d.m
》 x=0:0.1:10*pi;
》 fi=sin(x);
》 fi(100)=100;
》 [fo, ip] = func_despike_phasespace3d( fi, 9 );		<- plot process and length(fo)=length(fi)
							fo contains NaN
》 [fo, ip] = func_despike_phasespace3d( fi, 9, 1 );	<- plot process and length(fo)<length(fi)
							NaN is excluded from fo
》 [fo, ip] = func_despike_phasespace3d( fi, 9, 2 );	<- plot process and length(fo)=length(fi)
							NaN in fo is interpolated
>> Number of outlier   = 1 : Number of iteration = 1

(2) func_despike_phasespace3d_3var
Similar to (1) but no input for i_plot.

%======================================================================
References
%======================================================================

- Mori, N., T. Suzuki and S. Kakuno (2007) Noise of acoustic Doppler velocimeter data in bubbly flow, Journal of Engineering Mechanics, American Society of Civil Engineers, Vol.133, Issue 1, pp.122-125. (doi:10.1061/(ASCE)0733-9399(2007)133:1(122))
- http://www.oceanwave.jp/softwares/mace/index.php?MACE%20Softwares

%======================================================================
Update history
%======================================================================

2014/03/19	func_despike_phasespace3d.m has updated by Joseph Z. Ulanowski
2009/06/09	func_despike_phasespace3d has been modified to fill interpolated data
2009/06/09	readme.txt revised

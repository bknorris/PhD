Updates on ADV Processing script:

This run script has a series of pre-processing to provide quality control for acoustic doppler velocimeter data.
We follow the procedures originally set out by K.J. Rosenberger, USGS for post processing, and have added in 
additional parameters for signal noise reduction based on the findings of B. Gunawan & V. Neary (2011). This 
document attempts to explain the processing in the order in which it is applied to the raw data:

Steps:

1. Define metadata; metadata fields for the instrument deployment are defined. Additional fields for processing
   scripts are also supplied, such as whether or not to zero pitch and/or roll, apply a beam swap to align beams 
   accordinng to the orientation of the deployed instrument, and specifying cutoff amplitudes for final QC.
2. Data are loaded in using a script 'readADVdata.m' which additonally calls 'readVEChdr.m'. These scripts load 
   the text files created from Nortek's software and populate metadata fields from the instrument header file,
   appropriate and load sensor data into a structure, and appropriate and load velocity data into another 
   structure. Lastly, these data are combined into a final structure, ADV, that contains all the metadata,
   sensor data, and velocity data. 
3. Pressure signal adjustment; The pressure signal is adjusted firstly by removing atmospheric pressure. This is
   completed by splining external pressure sensor data (in our case, a weather station) and is subtracted from the
   ADV pressure signal. Then, the script prompts the user to select points where the instrument was out of the 
   water to designate pressure samples that would be affected by solar heating. These samples are defined, and the
   script computes the r-squared value of pressure versus temperature at these sample points. If the statistical
   fit is deemed sufficient (>90%), a pressure-temperature adjustment is applied to the pressure data to subtract
   the slope of the fit from out-of-water pressure values. 
4. A speed-of-sound correction is applied to the velocity data based on A.B. Coppens (1981).
5. Conversion to beam coordinates and Despiking; based on the findings of Gunawan & V.Neary (2011), velocity data
   are converted to beam coordinates. Then, a series of filtering software is run to replace the data that are 
   defined as outliers. Here, despiking is done on all of the velocities, including the times when the instrument
   was out of the water as this data will be thrown out eventually anyway.
6. Coordinate transformation; velocities in beam coordinates are transformed to Earth coordinates using the method
   outlined by Rusello (2009). 
7. Amplitude check; any velocity samples corresponding to amplitudes below a certain (user-defined) threshold are
   given NaN values. This will remove all out-of-water velocity data from the final dataset. During this time, 
   statistics (stdev, min and max) are calcuated for the velocity data (U,V,W) and the pressure data. 
8. File is saved.


Citations:

Coppens, Alan B. "Simple equations for the speed of sound in Neptunian waters." The Journal of the Acoustical 
Society of America 69.3 (1981): 862-863.

Gunawan, Budi, Vincent S. Neary, and James R. McNutt. "ORNL ADV post-processing guide and MATLAB algorithms 
for MHK site flow and turbulence analysis." ORNL/TML-2011/338, September (2011).

Rusello, P. J. "A practical primer for pulse coherent instruments." Nortek technical note No.: TN-027 (2009).

# Faraday cube - Readme file
# Outlier detection in 3 dimensional data with constraints

In radio astronomy it is pursued to detect so-called polarised sources which will help to map magnetic fields in galaxies 
(for the interested there is an additional file describing what this is). 
The astronomers are obtaining their data from radio telescopes which are stored in FITS-files 
(Files used by the large mass of astronomers, these files can be managed in Matlab and Python and probably in more languages). 
These FITS files contains a lot of information which makes it impractical to obtain information manually. 
The observational data contains observation from radio telescopes for several hours of selected parts of the sky 
in an interval of radio frequencies with different resolution. This data is saved in a "frequency cube" 
(the facing plane is the observed part of the sky while the third dimension represent the frequency).
To detect polarised radio sources one makes a so called Fourier-Faraday transformed from the frequency domain to the polarisation domain Φ (phi).
In these (new) cubes, called Faraday cubes, we search for polarised sources which are characterised by a large SNR (signal-to-noise ratio, with values of 5 or 6 and more) and having Φ-values not in interval -3 to 1 (Assumption: High SNR for Φ in interval -3 to 1 are not physical detections of polarised sources, they come from the telescope and instruments which are called instrumental polarisation and should thus be ignored). 
The idea here is to find a time-efficient and memory saving code for detection of polarised sources in the Faraday cubes. The interested will only work on a code detecting outliers in the Faraday cube and not be concerned with the transformation.


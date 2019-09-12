# Detection-of-polarised-sources: Outlier detection in 3 dimensional data with constraints

In radio astronomy it is pursued to detect so-called polarised sources which will help to map magnetic fields in galaxies (for the interested there is an additional file describing what this is). The astronomers are obtaining their data from radio telescopes which are stored in FITS-files (Files used by the large mass of astronomers, these files can be managed in Matlab and Python and probably in more languages). These FITS files contains a lot of information which makes it impractical to study by hand. The observational data contains observation for several hours of selected parts of the sky in an interval of radio frequencies with different resolution. This data is a "frequency cube".
To detect polarised radio sources they are so called Fourier-Faraday transformed from frequency domain to a different domain Φ (phi).
In these cubes (Faraday cubes) we search for polarised sources which are characterised by a large SNR (signal-to-noise ratio, 5 or 6) and having Φ-values not in interval -3 to 1 (Assumption: High SNR for Φ in interval -3 to 1 are not physical detections of polarised sources, they come from the telescope and instruments which are called instrumental polarisation and should thus be ignored). 
The idea here is to find a time-efficient and memory saving code for detection of polarised sources in the Faraday cubes. 

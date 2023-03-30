# assess-coordinate-transformation
ERT Assessment Challenge -- coordinate transformation

Use "make" to compile the code. 
Use "make run" to run the code using the file "GISradartest.exe."

Initial and final latitude and longitude values are stored in the file "GIS2Rinputs.txt."

This code uses the Vincenty algorithm[1] to transform between radar coordinates (bearing, range, initial latitude and longitude) and GIS coordinates (latitude, longitude for two points). The header declaration is as suggested in the problem statement. Vincenty's algorithm, which is iterative, may fail to converge for nearly-antipodal points, but converges very rapidly for the sample points given. When inverted, it gives the initial coordinates to very high precision.

When assessing this problem, I considered three possible solutions. In ascending order of time commitment and accuracy, they are:
  1. Spherical-Earth approximation.
  2. Vincenty's algorithm.
  3. Karney's algorithms for geodesics [2].
  
The spherical-Earth approximation presented as too prone to error (up to 0.5% for certain inputs), and Karney's algorithm would have involved a significant time commitment for proper implementation from scratch. Thus, Vincenty's algorithm was chosen.

[1] https://www.tandfonline.com/doi/abs/10.1179/sre.1975.23.176.88 , also available at https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

[2] https://arxiv.org/pdf/1109.4448.pdf

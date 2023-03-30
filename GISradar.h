#ifndef GISradar
#define GISradar

/*
 * RtoG takes radar coordinates (range, bearing) and latitude/longitude of the 
 * initial point (glonInit, glatInit) and returns the GIS coordinates of the 
 * final point (glonFinal, glatFinal).
 * 
 * GIS2Radar converts from GIS quantities (latitude and longitude) to radar coordinates 
 * (range and bearing) from the first point to the second.
 * 
 * Radar bearing is decimal degrees, clockwise from north. 
 * 
 * Latitude and longitude must be in decimal degrees.
 * 
*/


void GIS2Radar(double *range, double *bearing, double glonInit, double glatInit, double glonFinal, double glatFinal);
void RtoG(double range,  double bearing, double  glonInit, double glatInit, double *glonFinal, double *glatFinal);

#endif

#include <stdio.h>
#include <math.h>
#include "GISradar.h"

// int GIS2Radar(double *range, double *bearing, double glonInit, double glatInit, double glonFinal, double glatFinal);
// int RtoG (double range,  double bearing, double  glonInit, double glatInit, double *glonFinal, double *glatFinal);
void read_inverse_inputs(double* glonInit, double* glatInit, double* glonFinal, double* glatFinal);

int main(void){
    
    double glonInit;
    double glatInit;
    double glonFinal;
    double glatFinal;
    double range;
    double bearing;
    
    
    
    read_inverse_inputs(&glonInit, &glatInit, &glonFinal, &glatFinal);
    
//     printf("Initial/final lat/long: %lf %lf, %lf %lf\n", glatInit, glonInit, glatFinal, glonFinal);
    
    GIS2Radar(&range, &bearing, glonInit, glatInit, glonFinal, glatFinal);
    
    printf("range = %lf km, bearing = %lf deg\n", range, bearing);
    
    RtoG(range, bearing, glonInit, glatInit, &glonFinal, &glatFinal);
    
    printf("Test: reversing range and bearing. Should get original lat and long \nglatFinal = %lf deg, glonFinal = %lf deg\n", glatFinal, glonFinal);
    
    
}

void read_inverse_inputs(double* glonInit, double* glatInit, double* glonFinal, double* glatFinal){
    
    
    FILE* inputs = fopen("GIS2Rinputs.txt", "r");
    if(inputs == NULL){perror("inputs.txt does not exist");}
    
    char inputstring[200];
    
    if( fgets (inputstring, 200, inputs)!=NULL ) 
    {
      sscanf(inputstring, "%lf %lf %lf %lf", glatInit, glonInit, glatFinal, glonFinal);
      
      printf("initial lat and lon = %lf, %lf\n final lat and long = %lf, %lf\n", *glatInit, *glonInit, *glatFinal, *glonFinal);
    }
    
}

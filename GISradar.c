// This code implements Vincenty's method for solving the geodesic on an ellipsoid of rotation, from Survey Review 23, 176, 1975.


#include <stdio.h>
#include <math.h>


const double a = 6378.1370; //a, f, and b are all constants related to the shape of the Earth (WGS84).
const double f = 1/298.257223563;
const double b = 6356.752314245;
const double degtoRads = M_PI/180.0; //degree to radians conversion
const double eps = 10E-12; //The desired precision to which to converge.


int GIS2Radar(double *range, double *bearing, double glonInit, double glatInit, double glonFinal, double glatFinal){
    
    
    double phi1 = glatInit*degtoRads;
    double phi2 = glatFinal*degtoRads;
    double lambda, lambda_new, sinsigma, cossigma, cos2alpha, sin2alpha, cos2sigmam, U1, U2, u2, sigma, sinalpha;
    double L = (glonFinal - glonInit)*degtoRads;
    double xx, yy, C;
    
//     int tracker = 0;
    
    lambda = L;
    
    if(lambda >= M_PI){
     
        printf("Lambda > pi, will not converge");
        
        return 1;
        
    }
    
    lambda_new = 0;
    U1 = atan((1.0 - f)*tan(phi1));
    U2 = atan((1.0 - f)*tan(phi2));
    
    while(fabs(lambda_new - lambda) > eps){
        
        lambda = lambda_new;
        
        xx = pow(cos(U2)*sin(lambda),2.0);
        yy = pow(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda), 2.0);
        

        
        sinsigma = pow(xx + yy, 0.5);
        cossigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
        
        sigma = atan2(sinsigma, cossigma);
        
        sinalpha = cos(U1)*cos(U2)*sin(lambda)/sinsigma;
        
        sin2alpha = pow(sinalpha, 2.0);
        cos2alpha = 1.0 - sin2alpha;
        
        cos2sigmam = cossigma - 2.0*sin(U1)*sin(U2)/(cos2alpha);
        
        C = (f/16.0)*cos2alpha*(4.0 + f*(4.0 - 3.0*cos2alpha));
        
        lambda_new = L + (1.0 - C)*f*sinalpha*(sigma + C*sinsigma*(cos2sigmam + C*cossigma*(-1.0+ 2.0*pow(cos2sigmam, 2.0))));
        
//         ++tracker;
        
//         if(tracker > 100){printf("GIS2radar %d %lf\n", tracker, lambda_new - lambda);}
        
    }
        
    lambda = lambda_new;
    u2 = (1.0 - sin2alpha)*(a*a - b*b)/(b*b);
//     double A = 1.0 + (u2/16384.0)*(4096.0 + u2*(-768.0 + u2*(320.0 - 175.0*u2)));
//     double B = (u2/1024.0)*(256.0 + u2*(-128.0 + u2*(74.0 - 47.0*u2)));
    
    double A = 1.0 + u2/256.0*(64.0 + u2*(12.0+5*u2));
    double B = u2/512.0*(128.0 + u2*(-64.0 + 37.0*u2));
    
    sigma = atan2(sinsigma, cossigma);

    double temp1 = cossigma*(-1.0+ 2.0*pow(cos2sigmam, 2.0));
    double temp2 = -(1.0/6.0)*B*cos2sigmam*(-3.0 + 4.0*pow(sinsigma, 2.0))*(-3.0 + 4.0*pow(cos2sigmam, 2.0));
        
    double dsigma = B*sinsigma*(cos2sigmam + 0.25*B*(temp1 + temp2));
    
//     double dsigma = B*sinsigma*(cos2sigmam + 0.25*B*cossigma*(-1.0 + 2.0*pow(cos2sigmam, 2.0)));
    
    *range = b*A*(sigma - dsigma);
    
    yy = cos(U2)*sin(lambda);
    xx = cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda);
    
    double alpha1 = atan2(yy, xx);
    
    *bearing = alpha1/degtoRads;
    
    return 0;
    
//     printf("iterations %d\n", tracker);
}

int RtoG (double range,  double bearing, double  glonInit, double glatInit, double *glonFinal, double *glatFinal){
    
    double alpha1 = bearing*degtoRads;
    double phi1 = glatInit*degtoRads;
    double alpha, u2, sigma, sigma1, sigma_m, sigma_new, sinalpha, U1, U, A, B, dsigma, phi2, lambda2, cossigma, sinsigma, cosalpha, C;
    double s = range;
    double temp1, temp2, cos2sigmam;
    
    int tracker = 0;
    
    U1 = atan((1.0 - f)*tan(phi1));
    sigma1 = atan2((1.0 - f)*tan(phi1), cos(alpha1));
    sinalpha = cos(U1)*sin(alpha1);
       
    u2 = (1.0 - sinalpha*sinalpha)*(a*a - b*b)/(b*b);
    
//     A = 1 + (u2/16384.0)*(4096.0 + u2*(-768.0 + u2*(320.0 - 175.0*u2)));
//     B = (u2/1024.0)*(256.0 + u2*(-128.0 + u2*(74.0 - 47.0*u2)));
    A = 1.0 + u2/256.0*(64.0 + u2*(12.0+5*u2));
    B = u2/512.0*(128.0 + u2*(-64.0 + 37.0*u2));
    
    sigma_new = s /(b * A);
    
    while(fabs(sigma_new - sigma) > eps){
     
        sigma = sigma_new;
        cossigma = cos(sigma);
        sinsigma = sin(sigma);
        sigma_m = sigma1 + 0.5*(sigma); //eqn 5
        
        cos2sigmam = cos(2.0*sigma_m);
        
        temp1 = cos(sigma)*(-1.0 + 2.0*pow(cos2sigmam, 2.0));
        temp2 = -(1.0/6.0)*B*cos2sigmam*(-3.0 + 4.0*pow(sin(sigma), 2))*(-3.0 + 4.0*pow(cos2sigmam, 2));
        
        dsigma = B*sin(sigma)*(cos2sigmam + 0.25*B*(temp1 + temp2));
        
//         dsigma = B*sinsigma*(cos2sigmam + 0.25*B*cossigma*(-1.0 + 2.0*pow(cos2sigmam, 2.0)));
        
        sigma_new = s/(b * A) + dsigma;
        
        ++tracker;
        
//         if(tracker > 100){printf("RtoG %d\n", tracker);}
    }
    
//     sigma = sigma_new;
    
    double y = sin(U1)*cos(sigma) + cos(U1)*sin(sigma)*cos(alpha1);
    double x = (1.0 - f)*pow((sinalpha*sinalpha + pow((sin(U1)*sin(sigma) - cos(U1)*cos(sigma)*cos(alpha1)), 2.0)), 0.50);
    
    phi2 = atan2(y, x);
    
    y = sin(sigma)*sin(alpha1);
    x = cos(U1)*cos(sigma) - sin(U1)*sin(sigma)*cos(alpha1);
    
    cosalpha = cos(alpha1);
    C = (f/16.0)*cosalpha*cosalpha*(4.0 + f*(4.0 - 3.0*cosalpha*cosalpha));
    
    lambda2 = atan2(y,x) - (1.0 - C)*f*sinalpha*(sigma + C*sinsigma*(cos2sigmam + C*cossigma*(-1.0+ 2.0*pow(cos2sigmam, 2.0))));
    
    *glatFinal = phi2/degtoRads;
    *glonFinal = (glonInit) + (lambda2)/degtoRads;
    
//     printf("iterations %d\n", tracker);
    
    return 0;
    
}

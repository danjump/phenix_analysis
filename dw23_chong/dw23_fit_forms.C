#include "dw23_common.h"


// x[0] <- dw23
// this assumes a normalized, coaxial, double gaussian distribution of dw23
double double_gaus_1d_coax(double * x, double * par) {
  double offset = par[0];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double pi; pi = 3.14159265358979323846;
  double val = 1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-offset)/sigma1,2)) + factor*exp(-0.5*pow((x[0]-offset)/sigma2,2)) );
  return val;
}

// x[0] <- dw23
// this assumes a normalized, double gaussian distribution of dw23
double double_gaus_1d(double * x, double * par) {
  double offset = par[0];
  double mean_offset = par[4];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double pi; pi = 3.14159265358979323846;
  double val = 1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-(offset))/sigma1,2)) + factor*exp(-0.5*pow((x[0]-(offset+mean_offset))/sigma2,2)) );
  return val;
}

double single_gaus_slice(double *x, double *par) {
  double offset = par[0];
  double sigma = par[1];
  double scaling = par[2];
  double poly_scale = par[3];
  double pi; pi = 3.14159265358979323846;
  double val = scaling*poly_scale*1/(sigma*sqrt(2*pi)) * exp(-0.5*pow((x[0]-offset)/sigma,2));
  return val;
}

// x[0] <- dw23
// this assumes coaxial, double gaussian distribution of dw23
double double_gaus_slice_coax(double * x, double * par) {
  double offset = par[0];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double poly_scale = par[4];
  double pi; pi = 3.14159265358979323846;
  double val = poly_scale*1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-offset)/sigma1,2)) + factor*exp(-0.5*pow((x[0]-offset)/sigma2,2)) );
  return val;
}

// x[0] <- dw23
// this assumes a normalized, double gaussian distribution of dw23
double double_gaus_slice(double * x, double * par) {
  double offset = par[0];
  double mean_offset = par[5];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double poly_scale = par[4];
  double pi; pi = 3.14159265358979323846;
  double val = poly_scale*1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-(offset))/sigma1,2)) + factor*exp(-0.5*pow((x[0]-(offset+mean_offset))/sigma2,2)) );
  return val;
}

// x[0] <- dw23
// this assumes a normalized, tripple gaussian distribution of dw23
double tripple_gaus_slice(double * x, double * par) {
  double mean1 = par[0];
  double offset2 = par[5];
  double offset3 = par[6];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double sigma3 = par[7]; //should be a small bump
  double factor1 = par[3];
  double factor2 = par[8];
  double poly_scale = par[4];
  double pi; pi = 3.14159265358979323846;
  double val = poly_scale*1/( sigma1*sqrt(2*pi)+ factor1*sigma2*sqrt(2*pi) + factor2*sigma3*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-(mean1))/sigma1,2)) + factor1*exp(-0.5*pow((x[0]-(mean1+offset2))/sigma2,2)) + factor2*exp(-0.5*pow((x[0]-(mean1+offset3))/sigma3,2)) );
  return val;
}

// x[0] <- Wness
// x[1] <- dw23
// this assumes a normalized, coaxial, double gaussian distribution of dw23, whos parameters varry linearly with Wness
double double_gaus_2d_coax(double * x, double * par) {
  double gaus_mean = par[0] + par[1]*x[0];
  double gaus_sigma1 = par[2] + par[3]*x[0]; //should be wider
  double gaus_sigma2 = par[4] + par[5]*x[0]; //should be narrower
  double gaus_factor = par[6] + par[7]*x[0];
  double pi; pi = 3.14159265358979323846;

  // polynomial describes the normalized 1D wness distribution. 
  // used as a weight factor for the otherwise normalized dw23 distribution and different regions of wness
  double polynomial_factor = par[8] + par[9]*x[0] + par[10]*pow(x[0],2) + par[11]*pow(x[0],3) + par[12]*pow(x[0],4);
  
  double val = polynomial_factor*( 1/(gaus_sigma1*sqrt(2*pi)+ gaus_factor*gaus_sigma2*sqrt(2*pi)) ) * 
    ( exp(-0.5*pow((x[1]-gaus_mean)/gaus_sigma1,2)) + gaus_factor*exp(-0.5*pow((x[1]-gaus_mean)/gaus_sigma2,2)) );
  
  return val;
}
// x[0] <- Wness
// x[1] <- dw23
// this assumes a normalized, double gaussian distribution of dw23, whos parameters varry linearly with Wness
double double_gaus_2d(double * x, double * par) {
  double gaus_mean = par[0] + par[1]*x[0];
  double mean_offset = par[13] + par[14]*x[0]; //offset between each of the double gaussian parts
  double gaus_sigma1 = par[2] + par[3]*x[0]; //should be wider
  double gaus_sigma2 = par[4] + par[5]*x[0]; //should be narrower
  double gaus_factor = par[6] + par[7]*x[0];
  double pi; pi = 3.14159265358979323846;

  // polynomial describes the normalized 1D wness distribution. 
  // used as a weight factor for the otherwise normalized dw23 distribution and different regions of wness
  double polynomial_factor = par[8] + par[9]*x[0] + par[10]*pow(x[0],2) + par[11]*pow(x[0],3) + par[12]*pow(x[0],4);
  
  double val = polynomial_factor*( 1/(gaus_sigma1*sqrt(2*pi)+ gaus_factor*gaus_sigma2*sqrt(2*pi)) ) * 
    ( exp(-0.5*pow((x[1]-(gaus_mean))/gaus_sigma1,2)) + gaus_factor*exp(-0.5*pow((x[1]-(gaus_mean+mean_offset))/gaus_sigma2,2)) );
  
  return val;
}

// x[0] <- Wness
// x[1] <- dw23
// this assumes a normalized, double gaussian distribution of dw23, whos parameters varry linearly with Wness
double tripple_gaus_2d(double * x, double * par) {
  double mean1 = par[0] + par[1]*x[0];
  double offset2 =par[13] + par[14]*x[0];
  double offset3 = par[19] + par[20]*x[0];
  double sigma1 = par[2] + par[3]*x[0]; //should be wider
  double sigma2 = par[4] + par[5]*x[0]; //should be narrower
  double sigma3 = par[15] + par[16]*x[0]; //should be a small bump
  double factor1 = par[6] + par[7]*x[0];
  double factor2 = par[17] + par[18]*x[0];

  double pi; pi = 3.14159265358979323846;

  // polynomial describes the normalized 1D wness distribution. 
  // used as a weight factor for the otherwise normalized dw23 distribution and different regions of wness
  double polynomial_factor = par[8] + par[9]*x[0] + par[10]*pow(x[0],2) + par[11]*pow(x[0],3) + par[12]*pow(x[0],4);
  
  double val = polynomial_factor*1/( sigma1*sqrt(2*pi)+ factor1*sigma2*sqrt(2*pi) + factor2*sigma3*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-(mean1))/sigma1,2)) + factor1*exp(-0.5*pow((x[0]-(mean1+offset2))/sigma2,2)) + factor2*exp(-0.5*pow((x[0]-(mean1+offset3))/sigma3,2)) );
  
  return val;
}


void dw23_fit_forms(
    TH2F * h2_dw23_vs_wness_data[][2],
    TF1 * f_double_gaus_decomp[][2][2] 
    ) {
    
  


}



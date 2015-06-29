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


void dw23_2d_fits(bool do_plots,
    TH2F * h2_dw23_vs_wness[][2],
    TH1F * h_wness_slices[][2][10],
    TH1F * h_dw23_slices[][2][10],
    TF1 * f_double_gaus_coax[][2][10], 
    TF1 * f_double_gaus_offset[][2][10], 
    TF1 * f_tripple_gaus[][2][10],
    TF2 * f2_double_gaus_coax[][2], 
    TF2 * f2_double_gaus_offset[][2], 
    TF2 * f2_tripple_gaus[][2],
    float target_wness[][2],
    TF1 * f_double_gaus_coax_target[][2], 
    TF1 * f_double_gaus_offset_target[][2], 
    TF1 * f_tripple_gaus_target[][2]
    ) {
  
  char name[200];
  
  int narm=2, ncharge=2;
  const int num_methods = 3;
  //methods key:
  //0=double gaus coaxial
  //1=double gaus non-coax
  //2=tripple gaus
  int n_pars[num_methods] = {13,15,21};
  int n_pars_1d[num_methods] = {5,6,9};
  
  double fit_chi2ndf[num_methods][narm][ncharge];
  double mean_linear_pars[num_methods][narm][ncharge][2];
  double mean_offset_linear_pars[num_methods][narm][narm][2];
  double sigma1_linear_pars[num_methods][narm][ncharge][2];
  double sigma2_linear_pars[num_methods][narm][ncharge][2];
  double factor_linear_pars[num_methods][narm][ncharge][2];
  double sigma3_linear_pars[num_methods][narm][ncharge][2];
  double sigma3_factor_linear_pars[num_methods][narm][ncharge][2];
  double sigma3_offset_linear_pars[num_methods][narm][ncharge][2];
  double poly_pars[num_methods][narm][ncharge][5];

  double mean_linear_pars_err[num_methods][narm][ncharge][2];
  double mean_offset_linear_pars_err[num_methods][narm][ncharge][2];
  double sigma1_linear_pars_err[num_methods][narm][ncharge][2];
  double sigma2_linear_pars_err[num_methods][narm][ncharge][2];
  double factor_linear_pars_err[num_methods][narm][ncharge][2];
  double sigma3_linear_pars_err[num_methods][narm][ncharge][2];
  double sigma3_factor_linear_pars_err[num_methods][narm][ncharge][2];
  double sigma3_offset_linear_pars_err[num_methods][narm][ncharge][2];
  double poly_pars_err[num_methods][narm][ncharge][5];

  double dw23_vs_wness_integral;
  TH1F *h_wness[2][2];
  double wness_mean[2][2][10];
  TF1 *wness_pol[2][2];
  double integral[2][2];
 
  //start main loop of fitting process
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      dw23_vs_wness_integral = h2_dw23_vs_wness[arm][charge]->Integral("width");
      if(dw23_vs_wness_integral!=0)
        h2_dw23_vs_wness[arm][charge]->Scale(1.0/dw23_vs_wness_integral);

      //fit wness distribution with polynomial to get parameter seeds
      /*if(do_extrap_plots) {
        c_wness_fits->cd(arm+2*charge+1);
        gPad->SetLogy();
        }*/

      int firstbin, lastbin;
      firstbin = 1;
      lastbin = nbins_dw23;
      h_wness[arm][charge] = (TH1F*) h2_dw23_vs_wness[arm][charge]->ProjectionX(Form("h_wness_a%d_c%d",arm,charge),firstbin,lastbin);
      //h_wness[arm][charge]->Sumw2();
      integral[arm][charge] = h_wness[arm][charge]->Integral("width");
      if(integral != 0)
        h_wness[arm][charge]->Scale(1/integral[arm][charge]);
      /*if(do_extrap_plots)
        h_wness[arm][charge]->Draw("eP");
      if(save_fit_hists) {
        h_wness[arm][charge]->Write();
      }*/
      char name[200];
      sprintf(name,"wness_pol_0%d_charge%d_arm%d",0,arm,charge);
      wness_pol[arm][charge] = new TF1(name,"pol4(0)",0.1,0.9);
      h_wness[arm][charge]->Fit(name,"QLRN");
      TF1 *wness_pol_tmp = new TF1("wness_pol_tmp","pol4(0)",0,1);
      for(int i=0; i<5; i++)
        wness_pol_tmp->SetParameter(i,wness_pol[arm][charge]->GetParameter(i));
      wness_pol_tmp->SetLineColor(1);
      wness_pol_tmp->SetLineStyle(2);
      /*if(do_extrap_plots) {
        wness_pol_tmp->DrawCopy("same");
        if(arm==1 && charge==1)
          c_wness_fits->SaveAs("/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/dw23_extrapolation/temp/wness_fits.png");
      }*/
      TF2 *f2_temp[num_methods];
      for(int method=0; method < num_methods; method++) {//start loop for 2d dw23 fit

        if(method==0) {
          f2_temp[method] = new TF2("f2_temp_coax",double_gaus_2d_coax,0.1,0.9,-.3,.3,n_pars[method]);
        } else if(method==1) {
          f2_temp[method] = new TF2("f2_temp_offset",double_gaus_2d,0.1,0.9,-.3,.3,n_pars[method]);
        } else if(method==2){
          f2_temp[method] = new TF2("f2_temp_tripple",tripple_gaus_2d,0.1,0.9,-.3,.3,n_pars[method]);
        }

        //linear gaussian parameter seeds (from hide thesis mostly)
        /*f2_temp[method]->SetParameter(0,0.01);        // mean(offset) constant
        f2_temp[method]->SetParameter(1,0.01); // mean(offset) slope
        
        f2_temp[method]->SetParameter(2,0.1);               // sigma1 constant
        f2_temp[method]->SetParameter(3,-0.1);             // sigma1 slope
        
        f2_temp[method]->SetParameter(4,0.02);              // sigma2 constant
        f2_temp[method]->SetParameter(5,-0.1);             // sigma2 slope
        
        f2_temp[method]->SetParameter(6,1);                 // constant(factor) constant
        f2_temp[method]->SetParameter(7,-1);                 // constant(factor) slope*/

        double dw23_average = (h_dw23_slices[arm][charge][1]->GetMean()+
            h_dw23_slices[arm][charge][2]->GetMean()+
            h_dw23_slices[arm][charge][3]->GetMean())/3.0;

        f2_temp[method]->SetParameter(0,dw23_average-.005);        // mean(offset) constant
        f2_temp[method]->SetParameter(1,(charge?-1:1)*0.005); // mean(offset) slope
        
        
        f2_temp[method]->SetParameter(2,0.1);               // sigma1 constant
        f2_temp[method]->SetParameter(3,-0.06);             // sigma1 slope
        
        f2_temp[method]->SetParameter(4,0.04);              // sigma2 constant
        f2_temp[method]->SetParameter(5,-0.01);             // sigma2 slope
        
        f2_temp[method]->SetParameter(6,1);                 // constant(factor) constant
        f2_temp[method]->SetParameter(7,1);                 // constant(factor) slope

        if(method==1) {
          //the seed values were chosen by eye
          f2_temp[method]->SetParameter(13,(charge?1:-1)*.025);        // mean(offset) constant
          f2_temp[method]->SetParameter(14,(charge?-1:1)*.015); // mean(offset) slope
        } else if(method==2) {
          //the seed values were chosen by eye
          for(int i=0; i<15; i++) {
            f2_temp[method]->SetParameter(13,f2_temp[1]->GetParameter(i));
          }
          //f2_temp[method]->SetParameter(13,(charge?1:-1)*.025);        // mean(offset) constant
          //f2_temp[method]->SetParameter(14,(charge?-1:1)*.015); // mean(offset) slope
          f2_temp[method]->FixParameter(15,.01);                // sigma3
          f2_temp[method]->FixParameter(16,0);                  // sigma3 slope
          f2_temp[method]->FixParameter(17,.025);        // sigma3 factor
          f2_temp[method]->FixParameter(18,0); // sigma3 factor slope
          f2_temp[method]->FixParameter(19,(charge?-1:1)*.05);        // sigma3 offset
          f2_temp[method]->FixParameter(20,0); // sigma3 offset slope
        }
        
        
        //parameter limits also estimated from hides thesis for an initial fit
        /*f2_temp[method]->SetParLimits(0,dw23_average*0.5,dw23_average*2);        // mean(offset) constant
        f2_temp[method]->SetParLimits(1,(charge?-1:0)*0.01,(charge?0:1)*0.01);   // mean(offset) slope

        f2_temp[method]->SetParLimits(2,0.05,0.8);               // sigma1 constant
        f2_temp[method]->SetParLimits(3,-0.5,0);                 // sigma1 slope
        
        f2_temp[method]->SetParLimits(4,0.008,0.15);             // sigma2 constant
        f2_temp[method]->SetParLimits(5,-0.001,0);               // sigma2 slope
        
        f2_temp[method]->SetParLimits(6,-30,30);                 // constant(factor) constant
        f2_temp[method]->SetParLimits(7,-30,30);*/                 // constant(factor) slope
        
        //set wness parameters according to the previous 1d wness fit
        f2_temp[method]->FixParameter(8,wness_pol[arm][charge]->GetParameter(0));
        f2_temp[method]->FixParameter(9,wness_pol[arm][charge]->GetParameter(1));
        f2_temp[method]->FixParameter(10,wness_pol[arm][charge]->GetParameter(2));
        f2_temp[method]->FixParameter(11,wness_pol[arm][charge]->GetParameter(3));
        f2_temp[method]->FixParameter(12,wness_pol[arm][charge]->GetParameter(4));
        
        f2_temp[method]->SetParLimits(2,.09,1);
        f2_temp[method]->SetParLimits(4,0,.09);
        
        if(method==0) {
          h2_dw23_vs_wness[arm][charge]->Fit("f2_temp_coax","QLRN");
          f2_double_gaus_coax[arm][charge] = (TF2*) h2_dw23_vs_wness[arm][charge]->Clone(Form("f2_double_gaus_coax_a%d_c%d",arm,charge)); 
        } else if(method==1) {
          h2_dw23_vs_wness[arm][charge]->Fit("f2_temp_offset","QLRN");
          f2_double_gaus_offset[arm][charge] = (TF2*) h2_dw23_vs_wness[arm][charge]->Clone(Form("f2_double_gaus_offset_a%d_c%d",arm,charge)); 
        } else if(method==2) {
          h2_dw23_vs_wness[arm][charge]->Fit("f2_temp_tripple","QLRN");
          f2_tripple_gaus[arm][charge] = (TF2*) h2_dw23_vs_wness[arm][charge]->Clone(Form("f2_tripple_gaus_a%d_c%d",arm,charge)); 
        }

        fit_chi2ndf[method][arm][charge] = f2_temp[method]->GetChisquare() / f2_temp[method]->GetNDF();
        
        mean_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(0);
        mean_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(1);
        if(method==1) {
          mean_offset_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(13);
          mean_offset_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(14);
        } else if(method==2) {
          mean_offset_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(13);
          mean_offset_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(14);
          sigma3_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(15);
          sigma3_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(16);
          sigma3_factor_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(17);
          sigma3_factor_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(18);
          sigma3_offset_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(19);
          sigma3_offset_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(20);
        }
        sigma1_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(2);
        sigma1_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(3);
        sigma2_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(4);
        sigma2_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(5);
        factor_linear_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(6);
        factor_linear_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(7);
        poly_pars[method][arm][charge][0] = f2_temp[method]->GetParameter(8);
        poly_pars[method][arm][charge][1] = f2_temp[method]->GetParameter(9);
        poly_pars[method][arm][charge][2] = f2_temp[method]->GetParameter(10);
        poly_pars[method][arm][charge][3] = f2_temp[method]->GetParameter(11);
        poly_pars[method][arm][charge][4] = f2_temp[method]->GetParameter(12);
        
        mean_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(0);
        mean_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(1);
        if(method==1) {
          mean_offset_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(13);
          mean_offset_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(14);
        } else if(method==2) {
          mean_offset_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(13);
          mean_offset_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(14);
          sigma3_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(15);
          sigma3_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(16);
          sigma3_factor_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(17);
          sigma3_factor_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(18);
          sigma3_offset_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(19);
          sigma3_offset_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(20);
        }
        sigma1_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(2);
        sigma1_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(3);
        sigma2_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(4);
        sigma2_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(5);
        factor_linear_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(6);
        factor_linear_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(7);
        poly_pars_err[method][arm][charge][0] = f2_temp[method]->GetParError(8);
        poly_pars_err[method][arm][charge][1] = f2_temp[method]->GetParError(9);
        poly_pars_err[method][arm][charge][2] = f2_temp[method]->GetParError(10);
        poly_pars_err[method][arm][charge][3] = f2_temp[method]->GetParError(11);
        poly_pars_err[method][arm][charge][4] = f2_temp[method]->GetParError(12);
        
      }
    }
  }//end dw23 2d fitting loop
  
  // now project to the target area:
  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      for(int m=0; m<3; m++) {
        double pars[n_pars_1d[m]];

        pars[0] = mean_linear_pars[m][a][c][0] + target_wness[a][c]*mean_linear_pars[m][a][c][1];
        pars[1] = sigma1_linear_pars[m][a][c][0] + target_wness[a][c]*sigma1_linear_pars[m][a][c][1];
        pars[2] = sigma2_linear_pars[m][a][c][0] + target_wness[a][c]*sigma2_linear_pars[m][a][c][1];
        pars[3] = factor_linear_pars[m][a][c][0] + target_wness[a][c]*factor_linear_pars[m][a][c][1];
        pars[4] = poly_pars[m][a][c][0]
         + target_wness[a][c]*poly_pars[m][a][c][1]
         + pow(target_wness[a][c],2)*poly_pars[m][a][c][2]
         + pow(target_wness[a][c],3)*poly_pars[m][a][c][3]
         + pow(target_wness[a][c],4)*poly_pars[m][a][c][4];
        //pars[4]*=integral[a][c];
       
        if(m==0) {
          f_double_gaus_coax_target[a][c] = new TF1(Form("f_double_gaus_coax_target_a%d_c%d",a,c),double_gaus_slice_coax,-.3,.3,n_pars_1d[m]); 
          f_double_gaus_coax_target[a][c]->SetParameters(pars);
        } else if(m==1) {
          pars[5] = mean_offset_linear_pars[m][a][c][0] + target_wness[a][c]*mean_offset_linear_pars[m][a][c][1];
          f_double_gaus_offset_target[a][c] = new TF1(Form("f_double_gaus_offset_target_a%d_c%d",a,c),double_gaus_slice,-.3,.3,n_pars_1d[m]);
          f_double_gaus_offset_target[a][c]->SetParameters(pars);
        } else {
          pars[5] = mean_offset_linear_pars[m][a][c][0] + target_wness[a][c]*mean_offset_linear_pars[m][a][c][1];
          pars[6] = sigma3_offset_linear_pars[m][a][c][0] + target_wness[a][c]*sigma3_offset_linear_pars[m][a][c][1];
          pars[7] = sigma3_linear_pars[m][a][c][0]  + target_wness[a][c]*sigma3_linear_pars[m][a][c][1];
          pars[8] = sigma3_factor_linear_pars[m][a][c][0] + target_wness[a][c]*sigma3_factor_linear_pars[m][a][c][1];
          f_tripple_gaus_target[a][c] = new TF1(Form("f_tripple_gaus_target_a%d_c%d",a,c),tripple_gaus_slice,-.3,.3,n_pars_1d[m]);
          f_tripple_gaus_target[a][c]->SetParameters(pars);
        }
        
        
      }
    }
  }

  if(do_plots) {
    TCanvas *c_2d_params[2][2];
    TF1 *f_param_linear[2][2][num_methods][3];
    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        sprintf(name,"c_2d_params_a%d_c%d",a,c);
        c_2d_params[a][c] = new TCanvas(name,name,900,300);
        c_2d_params[a][c]->Divide(3);

        for(int m=0; m<num_methods; m++) {
          float mean_const = mean_linear_pars_err[m][a][c][0] + (m==0)?mean_offset_linear_pars_err[m][a][c][0]:0;
          float mean_slope = mean_linear_pars_err[m][a][c][1] + (m==0)?mean_offset_linear_pars_err[m][a][c][1]:0;
          sprintf(name,"f_mean2_param_linear_a%d_c%d_method%d",a,c,m);
          f_param_linear[a][c][m][0] = new TF1(name,"pol1(0)",0,1);
          f_param_linear[a][c][m][0]->SetParameter(0,mean_const);
          f_param_linear[a][c][m][0]->SetParameter(1,mean_slope);
          f_param_linear[a][c][m][0]->SetLineColor(m+1);
          c_2d_params[a][c]->cd(1);
          f_param_linear[a][c][m][0]->Draw(m==0?"":"same");
          
          sprintf(name,"f_sigma1_param_linear_a%d_c%d_method%d",a,c,m);
          f_param_linear[a][c][m][1] = new TF1(name,"pol1(0)",0,1);
          f_param_linear[a][c][m][1]->SetParameter(0,sigma1_linear_pars[m][a][c][0]);
          f_param_linear[a][c][m][1]->SetParameter(1,sigma1_linear_pars[m][a][c][1]);
          f_param_linear[a][c][m][1]->SetLineColor(m+1);
          c_2d_params[a][c]->cd(2);
          f_param_linear[a][c][m][1]->Draw(m==0?"":"same");
          
          sprintf(name,"f_sigma2_param_linear_a%d_c%d_method%d",a,c,m);
          f_param_linear[a][c][m][2] = new TF1(name,"pol1(0)",0,1);
          f_param_linear[a][c][m][2]->SetParameter(0,sigma2_linear_pars[m][a][c][0]);
          f_param_linear[a][c][m][2]->SetParameter(1,sigma2_linear_pars[m][a][c][1]);
          f_param_linear[a][c][m][2]->SetLineColor(m+1);
          c_2d_params[a][c]->cd(3);
          f_param_linear[a][c][m][2]->Draw(m==0?"":"same");
          
        }

      }
    }
  }//end do_plots conditional


}



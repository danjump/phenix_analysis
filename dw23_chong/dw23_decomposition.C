#include <cmath>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "dw23_common.h"

using namespace std;

//global constants for this file
int decomp_nthr=3;

//Function to separate central and sideband of dw23
void decompose_dw23(int fileindex, float wness_cut_low, float wness_cut_high, float threshold_factor, 
    TH2F *h2_dw23_vs_wness[][2], TH1F *h_dw23[][2], TH1F *h_dw23_sideband[][2], TH1F *h_dw23_central[][2], 
    TF1 * f_dw23[][2], TF1 *f_dw23_sideband[][2], TF1 *f_dw23_central[][2]);
//plotting function
void plot_sideband_extrap(TCanvas *c_dw23, TH1F *h_dw23[][2], TH1F *h_dw23_sideband[][2], TH1F *h_dw23_central[][2],
    TF1 *f_dw23[][2], TF1 *f_dw23_sideband[][2], TF1 *f_dw23_central[][2]);
//cint if statment necessary to allow 4+ dimension objects to be passes
#ifndef __CINT__
//feed in an array of already fitted functions to produce trend vs wness graphs.
void generate_trend_graphs(TGraphErrors *g_param_trends_sideband[][2][2], TGraphErrors *g_param_trends_central[][2][2], 
    TF1 *f_dw23_sideband[][2][2], TF1 *f_dw23_central[][2][2]);
//plotting function
void plot_trend_graphs(TCanvas * c_trends[][2], TGraphErrors * g_param_trends_sideband[][3][2][2], TGraphErrors * g_param_trends_central[][3][2][2]);
void generate_chisquare_hist(TH2F *h2_fulldw23chisqr_thr_vs_wness[][2],TF1 *f_dw23_data[][10][2][2],TH1F *h_dw23_data[][10][2][2]);

#endif

// x[0] <- dw23
double double_gaus_slice_decomp(double * x, double * par) {
  double mean1 = par[0];
  double sigma1 = par[1]; //should be wider
  double mean2 = par[2];
  double sigma2 = par[3]; //should be narrower
  double relative_factor = par[4];
  double pi; pi = 3.14159265358979323846;
  double val = 1/( sigma1*sqrt(2*pi)+ relative_factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-(mean1))/sigma1,2)) + relative_factor*exp(-0.5*pow((x[0]-(mean2))/sigma2,2)) );
  return val;
}


//MAIN FUNCTION -----------
void dw23_decomposition(
    bool do_plots,
    TH2F * h2_dw23_vs_wness_data[][2],
    TH1F * h_wness_slices[][2][10],
    TF1 * f_double_gaus_decomp[][2][10],
    TF1 * f_decomp_params[][6][2][2],
    float target_wness[][2],
    TF1 * f_double_gaus_decomp_target[][2][3]
    ) {
  

  //plot 2d dw23 vs wness
  if(do_plots && 0) {//plots
    TCanvas * c_dw23_vs_wness_data = new TCanvas("c_dw23_vs_wness_data","Data dw23 vs wness",1000,700);
    c_dw23_vs_wness_data->Divide(2,2);
    
    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        c_dw23_vs_wness_data->cd(a*2+c+1);
        h2_dw23_vs_wness_data[a][c]->Draw();
      }
    }
  }//plots
  
  TH1F *h_dw23_data[decomp_nthr][10][2][2];
  TF1 *f_dw23_data[decomp_nthr][10][2][2];
  TH1F *h_dw23_sideband_data[decomp_nthr][10][2][2];
  TF1 *f_dw23_sideband_data[decomp_nthr][10][2][2];
  TH1F *h_dw23_central_data[decomp_nthr][10][2][2];
  TF1 *f_dw23_central_data[decomp_nthr][10][2][2];

  TCanvas * c_dw23_data[decomp_nthr][10];
  
  for(int i=0; i<decomp_nthr; i++) {
    float thr = .80+.2*i;
    //separate dw23 sidebands in to one histogram, and central dw23 into another
    decompose_dw23(0, .1, .9, thr, h2_dw23_vs_wness_data, h_dw23_data[i][0], h_dw23_sideband_data[i][0], 
        h_dw23_central_data[i][0], f_dw23_data[i][0], f_dw23_sideband_data[i][0], f_dw23_central_data[i][0]);
    //decompose_dw23(1, .92, 1, 1.35, h2_dw23_vs_wness_w, h_dw23_w, h_dw23_sideband_w, h_dw23_central_w, f_dw23_sideband_w, f_dw23_central_w);

    for(int wness_bin=1; wness_bin<10; wness_bin++) {
      float wness_thresh_low = (wness_bin==9)?.92:(float)wness_bin/10.;
      float wness_thresh_hi = (float)wness_bin/10.+.1;
      decompose_dw23(0, wness_thresh_low, wness_thresh_hi, thr, 
          h2_dw23_vs_wness_data, h_dw23_data[i][wness_bin], h_dw23_sideband_data[i][wness_bin], h_dw23_central_data[i][wness_bin], 
          f_dw23_data[i][wness_bin], f_dw23_sideband_data[i][wness_bin], f_dw23_central_data[i][wness_bin]);
    }

    //plot 1D dw23 distributions for sideband extrapolation
    if(do_plots && 0) {//plots

      c_dw23_data[i][0] = new TCanvas(Form("c_dw23_data_full_thr%4.2f",thr),Form("thr%4.2f_full data dw23 extrap",thr),1000,700);
      plot_sideband_extrap(c_dw23_data[i][0], h_dw23_data[i][0], h_dw23_sideband_data[i][0], h_dw23_central_data[i][0],
          f_dw23_data[i][0], f_dw23_sideband_data[i][0], f_dw23_central_data[i][0] );

      for(int wness_bin=1; wness_bin<10; wness_bin++) {
        c_dw23_data[i][wness_bin] = new TCanvas(Form("c_dw23_data_w_thr%4.2f_%1.1f-%1.1f",thr,(double)wness_bin/10,(double)wness_bin/10+.1),
            Form("%4.2f %1.1f-%1.1f data dw23 extrap",thr,(double)wness_bin/10,(double)wness_bin/10+.1),1000,700);
        plot_sideband_extrap(c_dw23_data[i][wness_bin], h_dw23_data[i][wness_bin], h_dw23_sideband_data[i][wness_bin], h_dw23_central_data[i][wness_bin],
            f_dw23_data[i][wness_bin], f_dw23_sideband_data[i][wness_bin], f_dw23_central_data[i][wness_bin] );
      }
    }//plots
  }
  
  //graph trend in parameters
  TGraphErrors *g_param_trends_sideband[decomp_nthr][3][2][2];
  TGraphErrors *g_param_trends_central[decomp_nthr][3][2][2];
  TF1 * f_param_trends[decomp_nthr][6][2][2];
  float target_parameters[decomp_nthr][2][2][6];

  char * param_names[6] = {(char*)"nar_scale",(char*)"nar_mean",(char*)"nar_sigma",(char*)"wid_scale",(char*)"wid_mean",(char*)"wid_sigma"};
  char name[200];
  for(int i=0; i<decomp_nthr; i++) {
    generate_trend_graphs(g_param_trends_sideband[i], g_param_trends_central[i], f_dw23_sideband_data[i], f_dw23_central_data[i]);
    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        for(int p=0; p<6; p++) {
          if(p%3==0) {
            sprintf(name,"f_decomp_%s_a%d_c%d_thr%d",param_names[p],a,c,i);
            f_decomp_params[i][p][a][c] = new TF1(name,"pol14(0)",0,1);
            sprintf(name,"f_temp_decomp_%s_a%d_c%d_thr%d",param_names[p],a,c,i);
            f_param_trends[i][p][a][c] = new TF1(name,"pol4(0)",.1,1);
          } else {
            sprintf(name,"f_decomp_%s_a%d_c%d_thr%d",param_names[p],a,c,i);
            f_decomp_params[i][p][a][c] = new TF1(name,"pol1(0)",0,1);
            sprintf(name,"f_temp_decomp_%s_a%d_c%d_thr%d",param_names[p],a,c,i);
            f_param_trends[i][p][a][c] = new TF1(name,"pol1(0)",.1,.9);
          }
          if(p<3) {
            g_param_trends_central[i][p][a][c]->Fit(name,"QR");
          } else {
            g_param_trends_sideband[i][p-3][a][c]->Fit(name,"QR");
          }
          f_decomp_params[i][p][a][c]->SetParameters(f_param_trends[i][p][a][c]->GetParameters());
          
          float range_mean = h_wness_slices[a][c][9]->GetMean();
          target_parameters[i][a][c][p] = f_param_trends[i][p][a][c]->Eval(range_mean); 
          //printf("i%d a%dc%d p%d=%f\n",i,a,c,p,target_parameters[i][a][c][p]);
        }
      }
    }
  }
  
  if(do_plots) {
    TCanvas *c_trends[2][2];

    plot_trend_graphs(c_trends,g_param_trends_sideband,g_param_trends_central);
  }

  // assign selected TF1 to input TF1 object to pass out of this function
  for(int w=0; w<10; w++) {
    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        f_double_gaus_decomp[a][c][w] = (TF1*) f_dw23_data[2][w][a][c]->Clone(Form("f_double_gaus_decomp_wbin%d_a%d_c%d",w,a,c));
        if(w==9) {
          f_double_gaus_decomp_target[a][c][0] = (TF1*) new TF1(Form("f_double_gaus_decomp_target_a%d_c%d_narthr",a,c),double_gaus_slice_decomp,-.3,.3,5);
          f_double_gaus_decomp_target[a][c][0]->SetParameter(0,target_parameters[0][a][c][1]);
          f_double_gaus_decomp_target[a][c][0]->SetParameter(1,target_parameters[0][a][c][2]);
          f_double_gaus_decomp_target[a][c][0]->SetParameter(2,target_parameters[0][a][c][4]);
          f_double_gaus_decomp_target[a][c][0]->SetParameter(3,target_parameters[0][a][c][5]);
          f_double_gaus_decomp_target[a][c][0]->SetParameter(4,target_parameters[0][a][c][3]/target_parameters[0][a][c][0]);
          f_double_gaus_decomp_target[a][c][1] = (TF1*) new TF1(Form("f_double_gaus_decomp_target_a%d_c%d_midthr",a,c),double_gaus_slice_decomp,-.3,.3,5);
          f_double_gaus_decomp_target[a][c][1]->SetParameter(0,target_parameters[1][a][c][1]);
          f_double_gaus_decomp_target[a][c][1]->SetParameter(1,target_parameters[1][a][c][2]);
          f_double_gaus_decomp_target[a][c][1]->SetParameter(2,target_parameters[1][a][c][4]);
          f_double_gaus_decomp_target[a][c][1]->SetParameter(3,target_parameters[1][a][c][5]);
          f_double_gaus_decomp_target[a][c][1]->SetParameter(4,target_parameters[1][a][c][3]/target_parameters[1][a][c][0]);
          f_double_gaus_decomp_target[a][c][2] = (TF1*) new TF1(Form("f_double_gaus_decomp_target_a%d_c%d_widthr",a,c),double_gaus_slice_decomp,-.3,.3,5);
          f_double_gaus_decomp_target[a][c][2]->SetParameter(0,target_parameters[2][a][c][1]);
          f_double_gaus_decomp_target[a][c][2]->SetParameter(1,target_parameters[2][a][c][2]);
          f_double_gaus_decomp_target[a][c][2]->SetParameter(2,target_parameters[2][a][c][4]);
          f_double_gaus_decomp_target[a][c][2]->SetParameter(3,target_parameters[2][a][c][5]);
          f_double_gaus_decomp_target[a][c][2]->SetParameter(4,target_parameters[2][a][c][3]/target_parameters[2][a][c][0]);
        }
      }
    }
  }

  /*TH2F * h2_fulldw23chisqr_thr_vs_wness[2][2];
  generate_chisquare_hist(h2_fulldw23chisqr_thr_vs_wness,f_dw23_data,h_dw23_data);
  
  if(do_plots) {
    TCanvas * c_chi2 = new TCanvas("c_chi2","c_chi2",1000,700);
    
    c_chi2->Divide(2,2);

    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        c_chi2->cd(a*2+c+1);
        h2_fulldw23chisqr_thr_vs_wness[a][c]->Draw("COLZ");
      }
    }
  }*/
  
}// END MAIN FUNCTION -----------------






//Function to separate central and sideband of dw23.
//first characterize sidebands then subtract that from the whole
//histogram, then characterize the central part.
void decompose_dw23(int fileindex, float wness_cut_low, float wness_cut_high, float threshold_factor, 
    TH2F *h2_dw23_vs_wness[][2], TH1F *h_dw23[][2], TH1F *h_dw23_sideband[][2], TH1F *h_dw23_central[][2], 
    TF1 * f_dw23[][2], TF1 *f_dw23_sideband[][2], TF1 *f_dw23_central[][2]) {

  
  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      char * name = Form("h_dw23_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10));
      h_dw23[a][c] = (TH1F*)h2_dw23_vs_wness[a][c]->ProjectionY(name,
          h2_dw23_vs_wness[a][c]->GetXaxis()->FindBin(wness_cut_low+.0001),
          h2_dw23_vs_wness[a][c]->GetXaxis()->FindBin(wness_cut_high-.0001));
      h_dw23[a][c]->SetTitle(name);
      //h_dw23[a][c]->Sumw2();

      //temporarily copy the whole histogram to be modified later
      name = Form("h_dw23_central_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10));
      h_dw23_central[a][c] = (TH1F*)h_dw23[a][c]->Clone(name);
      
      //fit whole dw23 with a simple gaussian as a basis for chosing separation range
      f_dw23[a][c] = new TF1(Form("f_dw23_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10)),"gaus(0)",-.3,.3);
      f_dw23[a][c]->SetParameter(0,h_dw23[a][c]->GetMaximum());
      f_dw23[a][c]->SetParameter(1,h_dw23[a][c]->GetMean());
      f_dw23[a][c]->SetParameter(2,h_dw23[a][c]->GetRMS());
      h_dw23[a][c]->Fit(Form("f_dw23_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10)),"NQ");

      //create sidebands only histogram then fit with gaussian
      float lower_boundary = f_dw23[a][c]->GetParameter(1) - f_dw23[a][c]->GetParameter(2)*threshold_factor;
      float upper_boundary = f_dw23[a][c]->GetParameter(1) + f_dw23[a][c]->GetParameter(2)*threshold_factor;

      h_dw23_sideband[a][c] = new TH1F(
          Form("h_dw23_sideband_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10)),
          Form("h_dw23_sideband_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10)),
          nbins_dw23,-.3,.3);
      for(int i=1; i<=nbins_dw23; i++) {
        if( i <= h_dw23[a][c]->FindBin(lower_boundary) || i >= h_dw23[a][c]->FindBin(upper_boundary) )
          h_dw23_sideband[a][c]->SetBinContent(i,h_dw23[a][c]->GetBinContent(i));
        else
          h_dw23_sideband[a][c]->SetBinContent(i,0);
      }
      //h_dw23_sideband[a][c]->Sumw2();

      //set central histogram to full minus sideband fit
      char * sidename = Form("f_dw23_sideband_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10));
      f_dw23_sideband[a][c] = (TF1*)f_dw23[a][c]->Clone(sidename);
      h_dw23_sideband[a][c]->Fit(sidename,"MENQ");
      
      h_dw23_central[a][c]->Add(f_dw23_sideband[a][c],-1);
      
      char * centralname = Form("f_dw23_central_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10));
      f_dw23_central[a][c] = (TF1*)f_dw23[a][c]->Clone(centralname);
      //set the range to only the central area for the fit, then set it back to the full range
      //f_dw23_central[a][c]->SetRange(lower_boundary,upper_boundary);
      //f_dw23_central[a][c]->SetParLimits(2,0,.3);
      f_dw23_central[a][c]->SetParameter(0,h_dw23_central[a][c]->GetMaximum()/3);
      f_dw23_central[a][c]->SetParameter(2,f_dw23_central[a][c]->GetParameter(2)/2);
      h_dw23_central[a][c]->Fit(centralname,"MENRQ");
      f_dw23_central[a][c]->SetRange(-.3,.3);
    
      //redefine whole function as sum of the sideband and central gaussian functions
      f_dw23[a][c] = new TF1(Form("f_dw23_%s_thr%4.2f_a%d_c%d_%d_%d",processname[fileindex],threshold_factor,a,c,(int)floor(wness_cut_low*10),(int)floor(wness_cut_high*10)) ,
          Form("%s + %s", sidename, centralname), -.3, .3 );
    }
  }
  

}

void plot_sideband_extrap(TCanvas *c_dw23, TH1F *h_dw23[][2], TH1F *h_dw23_sideband[][2], TH1F *h_dw23_central[][2],
    TF1 * f_dw23[][2], TF1 *f_dw23_sideband[][2], TF1 *f_dw23_central[][2]) {
    
    c_dw23->Divide(2,2);
    
    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        c_dw23->cd(a*2+c+1);
        h_dw23[a][c]->SetLineColor(1);
        h_dw23[a][c]->Draw();
        f_dw23[a][c]->SetLineColor(1);
        f_dw23[a][c]->SetLineStyle(2);
        f_dw23[a][c]->Draw("same");
        h_dw23_sideband[a][c]->SetLineColor(3);
        h_dw23_sideband[a][c]->Draw("same");
        f_dw23_sideband[a][c]->SetLineColor(3);
        f_dw23_sideband[a][c]->SetLineStyle(2);
        f_dw23_sideband[a][c]->Draw("same");
        h_dw23_central[a][c]->SetLineColor(4);
        h_dw23_central[a][c]->Draw("same");
        f_dw23_central[a][c]->SetLineColor(4);
        f_dw23_central[a][c]->SetLineStyle(2);
        f_dw23_central[a][c]->Draw("sameL");
      }
    }
}

//cint if statment necessary to allow 4+ dimension objects to be passed
#ifndef __CINT__
void generate_trend_graphs(TGraphErrors *g_param_trends_sideband[][2][2], TGraphErrors *g_param_trends_central[][2][2], 
    TF1 *f_dw23_sideband[][2][2], TF1 *f_dw23_central[][2][2]) {
  const int n = 9;

  double x[n],y[n],xe[n],ye[n];

  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      for(int p=0; p<3; p++) {
        
        for(int w=1; w<=n; w++) {
          y[w-1] = f_dw23_sideband[w][a][c]->GetParameter(p);
          ye[w-1] = f_dw23_sideband[w][a][c]->GetParError(p);
          x[w-1] = (double)w/10+.05;
          xe[w-1] = .05;
        }
        g_param_trends_sideband[p][a][c] = new TGraphErrors(n,x,y,xe,ye);
        
        for(int w=1; w<=n; w++) {
          y[w-1] = f_dw23_central[w][a][c]->GetParameter(p);
          ye[w-1] = f_dw23_central[w][a][c]->GetParError(p);
        }
        g_param_trends_central[p][a][c] = new TGraphErrors(n,x,y,xe,ye);
      }
    }
  }
}

void plot_trend_graphs(TCanvas * c_trends[][2], TGraphErrors * g_param_trends_sideband[][3][2][2], TGraphErrors * g_param_trends_central[][3][2][2]) {
    char * param[3] = {(char*)"Gaus Amplitude",(char*)"Gaus Mean",(char*)"Gaus Sigma"};
    char * region[2] = {(char*)"Sideband ",(char*)"Central "};

    for(int a=0; a<2; a++) {
      for(int c=0; c<2; c++) {
        c_trends[a][c] = new TCanvas(Form("c_trends_a%d_c%d",a,c),Form("%d%dTrends in fit params",a,c),1000,700);
        c_trends[a][c]->Divide(3,2);
        
        for(int p=0; p<3; p++) {
          for(int t=0; t<decomp_nthr; t++) {
            c_trends[a][c]->cd(1+p);
            if(t==0) g_param_trends_sideband[t][p][a][c]->SetTitle(Form("%s%s",region[0],param[p]));
            g_param_trends_sideband[t][p][a][c]->SetLineColor(floor(t*15)+51);
            if(p>0 && t==0) {
              double ymin = g_param_trends_sideband[t][p][a][c]->GetYaxis()->GetXmin();
              double ymax = g_param_trends_sideband[t][p][a][c]->GetYaxis()->GetXmax();
              g_param_trends_sideband[t][p][a][c]->GetYaxis()->SetRangeUser(ymin, (ymax-ymin)*1.3+ymin);
            }
            g_param_trends_sideband[t][p][a][c]->Draw((t==0)?"ALP":"LP");
            c_trends[a][c]->cd(4+p);
            if(t==0) g_param_trends_central[t][p][a][c]->SetTitle(Form("%s%s",region[1],param[p]));
            g_param_trends_central[t][p][a][c]->SetLineColor(floor(t*15)+51);
            if(p>0 && t==0) {
              double ymin = g_param_trends_central[t][p][a][c]->GetYaxis()->GetXmin();
              double ymax = g_param_trends_central[t][p][a][c]->GetYaxis()->GetXmax();
              g_param_trends_central[t][p][a][c]->GetYaxis()->SetRangeUser(ymin, (ymax-ymin)*1.6+ymin);
            }
            g_param_trends_central[t][p][a][c]->Draw((t==0)?"ALP":"LP");
          }
        }
      }
    }
}

void generate_chisquare_hist(TH2F *h2_fulldw23chisqr_thr_vs_wness[][2],TF1 *f_dw23_data[][10][2][2],TH1F *h_dw23_data[][10][2][2]) {
  for(int a=0; a<2; a++) {//loops
    for(int c=0; c<2; c++) {
      h2_fulldw23chisqr_thr_vs_wness[a][c] = new TH2F(Form("h_fulldw23chisqr_thr_vs_wness_a%d_c%d",a,c),
          Form("h_fulldw23chisqr_thr_vs_wness_a%d_c%d",a,c),10,0,1.0,decomp_nthr,.85,.85+.05*decomp_nthr);
      for(int w=1; w<10; w++) {
        for(int t=0; t<decomp_nthr; t++) {
          double x[nbins_dw23],y[nbins_dw23],xe[nbins_dw23],ye[nbins_dw23];
          for(int b=1; b<=nbins_dw23; b++) {
            y[b] = h_dw23_data[t][w][a][c]->GetBinContent(b);
            ye[b] = h_dw23_data[t][w][a][c]->GetBinError(b);
            x[b] = -.3 + (b-1+0.5)*(.6/nbins_dw23);
            xe[b] = (.6/nbins_dw23)*.5;
          }
          TGraphErrors *g_temp = new TGraphErrors(nbins_dw23,x,y,xe,ye);
          double chi2 = g_temp->Chisquare(f_dw23_data[t][w][a][c]);
          h2_fulldw23chisqr_thr_vs_wness[a][c]->SetBinContent(w+1,t+1,chi2);
        }
      }
    }
  }//loops
}

#endif



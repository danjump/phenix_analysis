#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>
#include <TF2.h>
#include <time.h>
#include <TLeaf.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include <string>
#include <vector>
#include <RooFit.h>
#include "RooStats/SPlot.h"
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "sbg_fit.h"
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"
//#include "phys_hist_definitions.h"
#include "define_sbg_fit_hist.h"

using namespace RooFit ;
using namespace RooStats ;

// x[0] <- dw23
// this assumes a normalized, coaxial, double gaussian distribution of dw23
double double_gaus_1d(double * x, double * par) {
  double offset = par[0];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double pi; pi = 3.14159265358979323846;
  double val = 1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-offset)/sigma1,2)) + factor*exp(-0.5*pow((x[0]-offset)/sigma2,2)) );
  return val;
}

// x[0] <- dw23
// this assumes a normalized, coaxial, double gaussian distribution of dw23
double double_gaus_slice(double * x, double * par) {
  double offset = par[0];
  double sigma1 = par[1]; //should be wider
  double sigma2 = par[2]; //should be narrower
  double factor = par[3];
  double poly_scale = par[4];
  double pi; pi = 3.14159265358979323846;
  double val = poly_scale*1/( sigma1*sqrt(2*pi)+ factor*sigma2*sqrt(2*pi) ) * ( exp(-0.5*pow((x[0]-offset)/sigma1,2)) + factor*exp(-0.5*pow((x[0]-offset)/sigma2,2)) );
  return val;
}

// x[0] <- Wness
// x[1] <- dw23
// this assumes a normalized, coaxial, double gaussian distribution of dw23, whos parameters varry linearly with Wness
double double_gaus_2d(double * x, double * par) {
  double gaus_offset = par[0] + par[1]*x[0];
  double gaus_sigma1 = par[2] + par[3]*x[0]; //should be wider
  double gaus_sigma2 = par[4] + par[5]*x[0]; //should be narrower
  double gaus_factor = par[6] + par[7]*x[0];
  double pi; pi = 3.14159265358979323846;

  // polynomial describes the normalized 1D wness distribution. 
  // used as a weight factor for the otherwise normalized dw23 distribution and different regions of wness
  double polynomial_factor = par[8] + par[9]*x[0] + par[10]*pow(x[0],2) + par[11]*pow(x[0],3) + par[12]*pow(x[0],4);
  
  double val = polynomial_factor*( 1/(gaus_sigma1*sqrt(2*pi)+ gaus_factor*gaus_sigma2*sqrt(2*pi)) ) * 
    ( exp(-0.5*pow((x[1]-gaus_offset)/gaus_sigma1,2)) + gaus_factor*exp(-0.5*pow((x[1]-gaus_offset)/gaus_sigma2,2)) );
  
  return val;
}

double get_eta_trigeff(int eta_range, int arm, int charge) {
  FILE * trig_eff_file = fopen("/direct/phenix+u/rseidl/FORWARD/roofit/db/totaltrigefficiencies2_new13.dat","r");

  if(trig_eff_file==NULL) {
    cout << "bad file\n";
    fclose(trig_eff_file);
    return -999;
  }

  char text_line[300];
  float unused;

  int line_count=0;
  int num_eta_ranges = 15;
  int which_line = arm*num_eta_ranges + eta_range;

  if(which_line>29) {
    cout << "bad line num\n";
    fclose(trig_eff_file);
    return -999;
  }

  float trig_eff[2];

  while(fgets(text_line, 300, trig_eff_file) != NULL) {
    sscanf(text_line,"%[^\n]\n",text_line);
    sscanf(text_line,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",&unused,&unused,&trig_eff[0],&unused,&unused,&trig_eff[1],&unused,&unused);

    if(line_count==which_line) {
      fclose(trig_eff_file);
      return trig_eff[charge];
    }

    line_count++;
  }

  cout << "effor pulling trigger efficiency line\n";
  fclose(trig_eff_file);
  return -999;
}

int bad_eta_count = 0;
double get_eta_trigeff(float eta, float trig_eff_factors[15]) {

  int eta_range=-1;
  if(1.1 <= eta && eta <= 2.6) {
    eta_range = floor((eta-1.1)*10);
  } else {
    //printf("bad eta %f\n",eta);
    bad_eta_count++;
    return 0;
  }
  if(eta_range==15) eta_range=14;

  return trig_eff_factors[eta_range];
}

void sbg_fit( const std::string bkg_had_infilename,
    const std::string sig_w_infilename,
    const std::string bkg_dy_infilename,
    const std::string bkg_light_infilename,
    const std::string bkg_onium_infilename,
    const std::string bkg_onlyz_infilename,
    const std::string bkg_openbottom_infilename,
    const std::string bkg_opencharm_infilename,
    const std::string bkg_whad_infilename,
    const std::string bkg_wtau_infilename,
    const std::string bkg_zsum_infilename) { 

  //bool save_fit_hists=false;
  //bool read_fit_hists=false;
  //bool read_dw23_hist=true;
  //bool do_extrap_plots=false;
  bool save_sbg_hists=false;
  int sbg_rpcdca_choice = 3;
  int sbg_fvtx_choice = 1;

  char name[300]; sprintf(name,"prevent unused variable error");
  
  Int_t
  Run_Number,  Evt_Number,  triggerbit,
  clockcross,  fvtx_cone;
  
  Float_t 
  Evt_bbcZ,    Wness,
  tree_charge,      pT,          pz,
  phi,         eta,         DG0,
  DDG0,        DG4,         chi2,
  DCA_z,       DCA_r,       dphi23,
  dw23,        Rpc1dca,     Rpc3dca,
  fvtx_dphi,   fvtx_dr,     fvtx_dtheta, fvtx_dr_dtheta;
  
  TFile *infile[11];
  std::string temp_str[11];

  temp_str[0] = bkg_had_infilename;
  temp_str[1] = sig_w_infilename;
  temp_str[2] = bkg_dy_infilename;
  temp_str[3] = bkg_light_infilename;
  temp_str[4] = bkg_onium_infilename;
  temp_str[5] = bkg_onlyz_infilename;
  temp_str[6] = bkg_openbottom_infilename;
  temp_str[7] = bkg_opencharm_infilename;
  temp_str[8] = bkg_whad_infilename;
  temp_str[9] = bkg_wtau_infilename;
  temp_str[10] = bkg_zsum_infilename;
  
  infile[0]= new TFile(bkg_had_infilename.c_str());
  infile[1]= new TFile(sig_w_infilename.c_str());
  infile[2] = new TFile(bkg_dy_infilename.c_str());
  infile[3] = new TFile(bkg_light_infilename.c_str());
  infile[4] = new TFile(bkg_onium_infilename.c_str());
  infile[5] = new TFile(bkg_onlyz_infilename.c_str());
  infile[6] = new TFile(bkg_openbottom_infilename.c_str());
  infile[7] = new TFile(bkg_opencharm_infilename.c_str());
  infile[8] = new TFile(bkg_whad_infilename.c_str());
  infile[9] = new TFile(bkg_wtau_infilename.c_str());
  infile[10] = new TFile(bkg_zsum_infilename.c_str());

  //load trigger efficiencies 
  float trig_eff_factor[2][2][15]; // 15=eta range bins
  for(int ebin=0; ebin < 15; ebin++) {
    for(int arm=0; arm < 2; arm++) {
      for(int charge=0; charge < 2; charge++) {
        trig_eff_factor[arm][charge][ebin] = get_eta_trigeff(ebin,arm,charge); 
        printf("e%d a%d c%d ef%f\n",ebin,arm,charge,trig_eff_factor[arm][charge][ebin]);
      }
    }
  }

  define_dw23_vs_eta_hists();
  
  TH2F * dw23_vs_eta = new TH2F("dw23_vs_eta","dw23_vs_eta",100,distmin[5],distmax[5],100,distmin[14],distmax[14]);
  
  int in_count[11][2][2],over_count[11][2][2],under_count[11][2][2],raw_target_wness_count[11][2][2];
  float trigeff_target_wness_count[11][2][2];


  //loop over the different input files and all events within each file to get distributions
  for(int file_index=0; file_index<11; file_index++) {//read loop over files
    for(int arm=0; arm<2; arm++) {
      for(int charge_index=0; charge_index<2; charge_index++) {
        in_count[file_index][arm][charge_index] = 0;
        over_count[file_index][arm][charge_index] = 0;
        under_count[file_index][arm][charge_index] = 0;
        raw_target_wness_count[file_index][arm][charge_index] = 0;
        trigeff_target_wness_count[file_index][arm][charge_index] = 0;
      }
    }

    printf("test2\n");
    
    printf("Reading file:\n%s\n",temp_str[file_index].c_str());
    define_wness_tree(infile[file_index]);

    wness_tree->SetBranchAddress("Run_Number",&Run_Number);
    wness_tree->SetBranchAddress("Evt_Number",&Evt_Number);
    wness_tree->SetBranchAddress("triggerbit",&triggerbit);
    wness_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
    wness_tree->SetBranchAddress("clockcross",&clockcross);
    wness_tree->SetBranchAddress("Wness",&Wness);
    wness_tree->SetBranchAddress("charge",&tree_charge);
    wness_tree->SetBranchAddress("pT",&pT);
    wness_tree->SetBranchAddress("pz",&pz);
    wness_tree->SetBranchAddress("phi",&phi);
    wness_tree->SetBranchAddress("eta",&eta);
    wness_tree->SetBranchAddress("DG0",&DG0);
    wness_tree->SetBranchAddress("DDG0",&DDG0);
    wness_tree->SetBranchAddress("DG4",&DG4);
    wness_tree->SetBranchAddress("chi2",&chi2);
    wness_tree->SetBranchAddress("DCA_z",&DCA_z);
    wness_tree->SetBranchAddress("DCA_r",&DCA_r);
    wness_tree->SetBranchAddress("dphi23",&dphi23);
    wness_tree->SetBranchAddress("dw23",&dw23);
    wness_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
    wness_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
    wness_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
    wness_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
    wness_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
    wness_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
    wness_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);


    int entries = wness_tree->GetEntries();

    int percent_done=0;
    int percent_incriment=20;
    int percent_done_previous;

    time_t rawtime;

    cout << "\nNumber of events:  " << entries << endl;
    cout << "Starting Main Event Loop...\n\n";
    time(&rawtime);
    printf("Start time:  %s",ctime(&rawtime));

    printf("test3\n");

    //Get events from file, fill necessary histograms
    for(int i=0; i<entries; i++) {//Events loop
      //loop progress command line output
      percent_done_previous=percent_done;
      percent_done=(int)floor((float)(i+1)/(float)entries*(float)100);
      if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
        printf("%3i%% done",percent_done);
        if(percent_done==100) {
          cout << "^_^";
          time( &rawtime );
          printf(" %s",ctime(&rawtime));
        } else {
          cout << "...";
          time( &rawtime );
          printf(" %s",ctime(&rawtime));
        }
      }

      wness_tree->GetEntry(i);

      // Set arm (0 = South, 1 = North)
      int arm = 0;
      if		(pz < 0) arm = 0;
      else if (pz > 0) arm = 1;


      int charge_index = 0;
      if		(tree_charge < 0) charge_index = 0;
      else if (tree_charge > 0) charge_index = 1;

      int rpcdca_condition = -1;
      if((Rpc1dca < 100) && (Rpc3dca < 100)) {
        rpcdca_condition = 2;
      } else if(Rpc1dca < 100) {
        rpcdca_condition = 0;
        //continue;
      } else if(Rpc3dca < 100) {
        rpcdca_condition = 1;
        //continue;
      }

      int goodfvtx = 0;
      if(-1.5 < fvtx_dphi && fvtx_dphi < 1.5 &&
          0 < fvtx_dtheta && fvtx_dtheta < 1.5  &&
          0 < fvtx_dr && fvtx_dr < 100) {
        bool fvtx_outofbounds = false;
        if(-1.5 > fvtx_dphi || fvtx_dphi > distmax[19]) {
          printf("fvtx_dphi: %f  ",fvtx_dphi);
          fvtx_outofbounds = true;
        }
        if(distmin[20] > fvtx_dr || fvtx_dr > distmax[20]) {
          printf("fvtx_dr: %f  ",fvtx_dr);
          fvtx_outofbounds = true;
        }
        if(distmin[21] > fvtx_dtheta || fvtx_dtheta > distmax[21]) {
          printf("fvtx_dtheta: %f  ",fvtx_dtheta);
          fvtx_outofbounds = true;
        }
        //if(distmin[17] > fvtx_dr_dtheta || fvtx_dr_dtheta > distmax[17]) {
          //printf("fvtx_dr_dtheta: %f",fvtx_dr_dtheta);
          //fvtx_outofbounds = true;
       // }

        if(fvtx_outofbounds) {
          printf("\n");
          continue;
        }

        goodfvtx = 1;

      } else {
        goodfvtx = 0;
      }

      int wness_section;
      if(Wness==1.0) wness_section = 9;
      else wness_section = floor(Wness*10);

      float trig_scaling_factor;
      if(file_index==0)
        trig_scaling_factor = 1;
      else
        trig_scaling_factor = get_eta_trigeff(eta,trig_eff_factor[arm][charge_index]);
      

      if(Wness >= 0.92 && file_index != 0) {
        h_sbg_fit_eta[file_index][arm][charge_index][rpcdca_condition][goodfvtx]->Fill(eta,trig_scaling_factor);
        h_sbg_fit_dw23[file_index][arm][charge_index][rpcdca_condition][goodfvtx]->Fill(dw23,trig_scaling_factor);
        h_target_wness[file_index][arm][charge_index]->Fill(Wness,trig_scaling_factor);
      } else if(Wness >= 0.1 && Wness <= 0.7 && file_index==0) {
        h_sbg_fit_eta[file_index][arm][charge_index][rpcdca_condition][goodfvtx]->Fill(eta);
      }
      
      // rpcdca=3 goodfvtx=1 is a shortcut to represent the combination of all rpcdca/fvtx conditions
      if(Wness >= 0.92 && file_index != 0) {
        raw_target_wness_count[file_index][arm][charge_index]++;
        trigeff_target_wness_count[file_index][arm][charge_index] += trig_scaling_factor;
        h_sbg_fit_eta[file_index][arm][charge_index][3][1]->Fill(eta,trig_scaling_factor);
        h_sbg_fit_dw23[file_index][arm][charge_index][3][1]->Fill(dw23,trig_scaling_factor);
      } else if(Wness >= 0.1 && Wness <= 0.7 && file_index==0) {
        raw_target_wness_count[file_index][arm][charge_index]++;
        trigeff_target_wness_count[file_index][arm][charge_index] += trig_scaling_factor;
        h_sbg_fit_eta[file_index][arm][charge_index][3][1]->Fill(eta);
      }
      
      // rpcdca=3 goodfvtx=0 is a shortcut to represent the combination of both fvtx conditions, but no rpc1only
      if(Wness >= 0.92 && file_index != 0 && Rpc3dca<100) {
        h_sbg_fit_eta[file_index][arm][charge_index][3][0]->Fill(eta,trig_scaling_factor);
        h_sbg_fit_dw23[file_index][arm][charge_index][3][0]->Fill(dw23,trig_scaling_factor);
      } else if(Wness >= 0.1 && Wness <= 0.7 && file_index==0) {
        h_sbg_fit_eta[file_index][arm][charge_index][3][0]->Fill(eta);
      }

      
      
      if(wness_section==9) dw23_vs_eta->Fill(eta,dw23,trig_scaling_factor);
      h_dw23[file_index][arm][charge_index][wness_section]->Fill(dw23,trig_scaling_factor);
      h_wness_sections[file_index][arm][charge_index][wness_section]->Fill(Wness,trig_scaling_factor);
      h_wness[file_index][arm][charge_index]->Fill(Wness,trig_scaling_factor);
      h_eta[file_index][arm][charge_index][wness_section]->Fill(eta,trig_scaling_factor);
      h2_dw23_vs_eta[file_index][arm][charge_index][wness_section]->Fill(eta,dw23,trig_scaling_factor);
      h2_dw23_vs_wness[file_index][arm][charge_index]->Fill(Wness,dw23,trig_scaling_factor);
      
      if(Wness>0 && Wness<1)
        in_count[file_index][arm][charge_index]++;
      else if(Wness>=1) {
        over_count[file_index][arm][charge_index]++;
      } else if(Wness <= 0) {
        under_count[file_index][arm][charge_index]++;
      }


    } //end main event loop
    printf("test1\n");
    printf("Num bad etas file %d: %d out of %d\n", file_index,bad_eta_count,entries); bad_eta_count = 0;
  }//end read loop over files
  
  

  //Signal/Background Ratio dw23 vs. eta fit
  TFile * mu_bg_hists_file = new TFile("temp_mu_bg_hists_file.root","RECREATE");
  TFile * sbg_hists_file;
  if(save_sbg_hists) {
    sbg_hists_file = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit/output/sbg_fit_pdfs.root","RECREATE");
  }
  
  define_wness_tree(infile[0]);
  TTree *temp_wness_tree;
  double sig_w_count[2][2],bkg_mu_count[2][2],bkg_had_count[2][2],data_count[2][2];
  double sig_w_err_lo[2][2],sig_w_err_hi[2][2],bkg_mu_err_lo[2][2],bkg_mu_err_hi[2][2],bkg_had_err_lo[2][2],bkg_had_err_hi[2][2];
  double mu_bg_yields[2][2][9];
//----
  //TFile * dw23_infile = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_extrapolation/");
  TFile * dw23_infile = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_chong/output.root");
  TH1F * h_temp[4][2][2];
  TCanvas * c_dw23_forms = new TCanvas("c_dw23_forms","c_dw23_forms",800,600);
  c_dw23_forms->Divide(2,2);

  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      c_dw23_forms->cd(1+arm*2+charge);
      TF1 * f_temp;
      double integral;

      sprintf(name,"f_double_gaus_decomp_target_a%d_c%d_narthr",arm,charge);
      f_temp = (TF1*)dw23_infile->Get(name);
      h_temp[2][arm][charge] = new TH1F(Form("temp_2_a%d_c%d",arm,charge),Form("temp_2_a%d_c%d",arm,charge),10*nhistbins_sbg_dw23,-.1,.1);
      h_temp[2][arm][charge]->FillRandom(name,10000000);
      integral = h_temp[2][arm][charge]->Integral("width");
      h_temp[2][arm][charge]->Scale(1.0/integral);
      h_temp[2][arm][charge]->SetLineColor(3);
      h_temp[2][arm][charge]->SetLineStyle(2);
      h_temp[2][arm][charge]->Draw();
      
      sprintf(name,"f_double_gaus_decomp_target_a%d_c%d_widthr",arm,charge);
      f_temp = (TF1*)dw23_infile->Get(name);
      h_temp[3][arm][charge] = new TH1F(Form("temp_3_a%d_c%d",arm,charge),Form("temp_3_a%d_c%d",arm,charge),10*nhistbins_sbg_dw23,-.1,.1);
      h_temp[3][arm][charge]->FillRandom(name,10000000);
      integral = h_temp[3][arm][charge]->Integral("width");
      h_temp[3][arm][charge]->Scale(1.0/integral);
      h_temp[3][arm][charge]->SetLineColor(3);
      h_temp[3][arm][charge]->SetLineStyle(2);
      h_temp[3][arm][charge]->Draw("same");
      
      sprintf(name,"f_double_gaus_coax_target_a%d_c%d",arm,charge);
      f_temp = (TF1*)dw23_infile->Get(name);
      h_temp[0][arm][charge] = new TH1F(Form("temp_0_a%d_c%d",arm,charge),Form("temp_0_a%d_c%d",arm,charge),10*nhistbins_sbg_dw23,-.1,.1);
      h_temp[0][arm][charge]->FillRandom(name,10000000);
      integral = h_temp[0][arm][charge]->Integral("width");
      h_temp[0][arm][charge]->Scale(1.0/integral);
      h_temp[0][arm][charge]->SetLineColor(1);
      h_temp[0][arm][charge]->Draw("same");
      
      sprintf(name,"f_double_gaus_offset_target_a%d_c%d",arm,charge);
      f_temp = (TF1*)dw23_infile->Get(name);
      h_temp[1][arm][charge] = new TH1F(Form("temp_1_a%d_c%d",arm,charge),Form("temp_1_a%d_c%d",arm,charge),10*nhistbins_sbg_dw23,-.1,.1);
      h_temp[1][arm][charge]->FillRandom(name,10000000);
      integral = h_temp[1][arm][charge]->Integral("width");
      h_temp[1][arm][charge]->Scale(1.0/integral);
      h_temp[1][arm][charge]->SetLineColor(2);
      h_temp[1][arm][charge]->Draw("same");
      
     
      sprintf(name,"h_sbg_fit_dw23_source0_a%d_c%d",arm,charge);
      h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]= (TH1F*)h_temp[1][arm][charge]->Clone(name); 
    }
  }

  
  
  //cross section in picobarns
  double mu_bg_cross_section[9] = {
    pow(10,7)*5.3200, //dy
    pow(10,10)*5.94000, //light
    pow(10,8)*1.3500, //onium
    //pow(10,2)*(-1.33), //onlyz old value
    pow(10,2)*(-3.37), //onlyz
    pow(10,6)*7.3000, //openbottom
    pow(10,8)*5.7100, //opencharm
    pow(10,3)*1.6600, //whad
    pow(10,3)*1.6600, //wtau
    pow(10,4)*1.5900}; //z_sum
  double mu_bg_generated_events[9] = { // all values are in k events
    /*6400000, //dy
    193600, //light
    55470000, //onium
    106500, //onlyz
    4003000, //openbottom
    134000000, //opencharm
    81000, //whad
    82000, //wtau
    245200}; //z_sum// run 12 sims */
    58420000, //dy
    7351500, //light
    150420000, //onium
    173100, //onlyz
    7363000, //openbottom
    584940000, //opencharm
    342000, //whad
    343000, //wtau
    292900}; //z_sum
  double mu_bg_dimuon_factor[9] = {
    1.4, //dy
    1.4, //light
    1.3, //onium
    //1, //onlyz old value
    2.2, //onlyz
    .77, //openbottom
    3.4, //opencharm
    1, //whad
    1, //wtau
    2.2}; //z_sum

  double mu_bg_scaling_factor[9] = {
    mu_bg_dimuon_factor[0]*mu_bg_cross_section[0]/(1000*mu_bg_generated_events[0])*228/1.4, //dy
    mu_bg_dimuon_factor[1]*mu_bg_cross_section[1]/(1000*mu_bg_generated_events[1])*228/1.4*0, //light
    mu_bg_dimuon_factor[2]*mu_bg_cross_section[2]/(1000*mu_bg_generated_events[2])*228/1.4, //onium
    mu_bg_dimuon_factor[3]*mu_bg_cross_section[3]/(1000*mu_bg_generated_events[3])*228/1.4, //onlyz
    mu_bg_dimuon_factor[4]*mu_bg_cross_section[4]/(1000*mu_bg_generated_events[4])*228/1.4, //openbottom
    mu_bg_dimuon_factor[5]*mu_bg_cross_section[5]/(1000*mu_bg_generated_events[5])*228/1.4, //opencharm
    mu_bg_dimuon_factor[6]*mu_bg_cross_section[6]/(1000*mu_bg_generated_events[6])*228/1.4, //whad
    mu_bg_dimuon_factor[7]*mu_bg_cross_section[7]/(1000*mu_bg_generated_events[7])*228/1.4, //wtau
    mu_bg_dimuon_factor[8]*mu_bg_cross_section[8]/(1000*mu_bg_generated_events[8])*228/1.4}; //z_sum

  TH1F * h_dw23_data[2][2];
  TH1F * h_eta_data[2][2];
  //TH1F *h_combined_mu_sig_dw23[2][2];
  //TH1F *h_combined_mu_sig_eta[2][2];
  TH1F *h_combined_mu_bg_dw23[2][2];
  TH1F *h_combined_mu_bg_eta[2][2];
  TH1F *h_sbg_pdf_sig_mu_eta[2][2];
  TH1F *h_sbg_pdf_sig_mu_dw23[2][2];
  TH1F *h_sbg_pdf_bkg_mu_eta[2][2];
  TH1F *h_sbg_pdf_bkg_mu_dw23[2][2];
  TH1F *h_sbg_pdf_bkg_had_eta[2][2];
  TH1F *h_sbg_pdf_bkg_had_dw23[2][2];
  //TH1F *h_combined_had_bg_dw23[2][2];
  //TH1F *h_combined_had_bg_eta[2][2];
  Float_t temp_dw23, temp_eta;
  
  TCanvas * c_sbg_fits_dw23 = new TCanvas("c_sbg_fits_dw23","c_sbg_fits_dw23",1500,1500);
  c_sbg_fits_dw23->Divide(2,2);
  TCanvas * c_sbg_fits_eta = new TCanvas("c_sbg_fits_eta","c_sbg_fits_eta",1500,1500);
  c_sbg_fits_eta->Divide(2,2);
  RooPlot *frame_dw23[2][2];
  RooPlot *frame_eta[2][2];
      
  RooRealVar * v_eta[2][2];
  RooRealVar * v_dw23[2][2];
  RooRealVar * v_wness[2][2];
  RooDataSet * data[2][2];

  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      v_eta[arm][charge] = new RooRealVar("eta","eta",1.1,2.6);
      v_dw23[arm][charge] = new RooRealVar("dw23","dw23",-.1,.1);
      v_wness[arm][charge] = new RooRealVar("Wness","Wness",0,1);

      char conditions[500];
      sprintf(conditions,"Wness>.92 && 16<pT && pT<60");
      if(sbg_rpcdca_choice==3) {
        
      } else {
        if(sbg_rpcdca_choice!=1)
          sprintf(conditions,"%s && Rpc1dca<100", conditions);
        if(sbg_rpcdca_choice!=0)
          sprintf(conditions,"%s && Rpc3dca<100", conditions);
        
        if(sbg_fvtx_choice==0)
          sprintf(conditions,"%s && fvtx_dr<-1", conditions);
        else
          sprintf(conditions,"%s && fvtx_dr>-1", conditions);
      }
      if(arm==0)
        sprintf(conditions,"%s && pz<0", conditions);
      else
        sprintf(conditions,"%s && pz>0", conditions);

      if(charge==0)
        sprintf(conditions,"%s && charge<0", conditions);
      else
        sprintf(conditions,"%s && charge>0", conditions);

      temp_wness_tree = (TTree*)wness_tree->CopyTree(conditions);
      data_count[arm][charge] = 0.0+temp_wness_tree->GetEntries();
      data[arm][charge] = new RooDataSet("data","data",temp_wness_tree,RooArgSet(*v_eta[arm][charge],*v_dw23[arm][charge],*v_wness[arm][charge]));
      
      temp_wness_tree->SetBranchAddress("dw23",&temp_dw23);
      temp_wness_tree->SetBranchAddress("eta",&temp_eta);
      
      sprintf(name,"h_dw23_data_arm%d_charge%d",arm,charge);
      h_dw23_data[arm][charge] = new TH1F(name,name,40,distmin[14],distmax[14]);
      h_dw23_data[arm][charge]->Sumw2();
      sprintf(name,"h_eta_data_arm%d_charge%d",arm,charge);
      h_eta_data[arm][charge] = new TH1F(name,name,40,distmin[5],distmax[5]);
      h_eta_data[arm][charge]->Sumw2();
      
      for(int event=0; event<temp_wness_tree->GetEntries(); event++) {
        temp_wness_tree->GetEntry(event);
        
        h_dw23_data[arm][charge]->Fill(temp_dw23);
        h_eta_data[arm][charge]->Fill(temp_eta);
        
      }

      //compose hadronic bkg dw23
      /*TF1 * f_temp;
      TH1F * h_temp;
      sprintf(name,"f_double_gaus_decomp_target_a%d_c%d",arm,charge);
      f_temp = (TF1*)dw23_infile->Get(name);
      h_temp = new TH1F(Form("temp_a%d_c%d",arm,charge),Form("temp_a%d_c%d",arm,charge),nhistbins_sbg_dw23,-.1,.1);
      h_temp->FillRandom(name,10000);
      double integral = h_temp->Integral("width");
      h_temp->Scale(1.0/integral);
      sprintf(name,"h_sbg_fit_dw23_source0_a%d_c%d",arm,charge);
      h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]= (TH1F*)h_temp->Clone(name); */


      //compose combined muonic background distributions
      //compose combined histograms
      sprintf(name,"h_combined_mu_bg_dw23_arm%d_charge%d",arm,charge);
      h_combined_mu_bg_dw23[arm][charge] = new TH1F(name,name,nhistbins_sbg_dw23,-0.1,0.1);
      sprintf(name,"h_combined_mu_bg_eta_arm%d_charge%d",arm,charge);
      h_combined_mu_bg_eta[arm][charge] = new TH1F(name,name,nhistbins_sbg_eta,1.1,2.6);
      h_combined_mu_bg_eta[arm][charge]->Sumw2();
      double mu_bg_yield=0;

      printf("Muon Bkg Yields:\n");
      for(int process=0; process<9; process++) {
        h_combined_mu_bg_eta[arm][charge]->Add(h_sbg_fit_eta[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice],mu_bg_scaling_factor[process]);
        h_combined_mu_bg_dw23[arm][charge]->Add(h_sbg_fit_dw23[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice],mu_bg_scaling_factor[process]);
        mu_bg_yield +=h_sbg_fit_dw23[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Integral()*mu_bg_scaling_factor[process];
        mu_bg_yields[arm][charge][process] = h_sbg_fit_dw23[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Integral()*mu_bg_scaling_factor[process];
        
        mu_bg_hists_file->cd();
        h_sbg_fit_eta[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Write();
        h_sbg_fit_dw23[2+process][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Write();
      }
      
      //add only z component to w signal distributions (according to scale factor from luminosity difference from 
      h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Add(h_sbg_fit_eta[3][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice],0.235);
      h_sbg_fit_dw23[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Add(h_sbg_fit_dw23[3][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice],0.235);

      cout << "arm" << arm << " charge" << charge << " bgeta pre:" << h_combined_mu_bg_eta[arm][charge]->Integral();
      h_combined_mu_bg_eta[arm][charge]->Smooth(6);
      cout << " post:" << h_combined_mu_bg_eta[arm][charge]->Integral() << endl;
      
      cout << "arm" << arm << " charge" << charge << " seta pre:" << h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Integral();
      h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Smooth(6);
      cout << " post:" << h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Integral() << endl;

      cout << endl;

      //first clone input hists to save to file
      
      
      //set up s/bg extraction pdfs
      RooDataHist * h_bkg_had_eta = new RooDataHist("h_bkg_had_eta","h_bkg_had_eta",RooArgSet(*v_eta[arm][charge]),h_sbg_fit_eta[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]);
      RooDataHist * h_bkg_had_dw23 = new RooDataHist("h_bkg_had_dw23","h_bkg_had_eta",RooArgSet(*v_dw23[arm][charge]),h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]);
      RooDataHist * h_bkg_mu_eta = new RooDataHist("h_bkg_mu_eta","h_bkg_mu_eta",RooArgSet(*v_eta[arm][charge]),h_combined_mu_bg_eta[arm][charge]);
      RooDataHist * h_bkg_mu_dw23 = new RooDataHist("h_bkg_mu_dw23","h_bkg_mu_dw23",RooArgSet(*v_dw23[arm][charge]),h_combined_mu_bg_dw23[arm][charge]);
      RooDataHist * h_sig_mu_eta = new RooDataHist("h_sig_mu_eta","h_sig_mu_eta",RooArgSet(*v_eta[arm][charge]),h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]);
      RooDataHist * h_sig_mu_dw23 = new RooDataHist("h_sig_mu_dw23","h_sig_mu_dw23",RooArgSet(*v_dw23[arm][charge]),h_sbg_fit_dw23[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]);

      RooHistPdf * pdf_bkg_had_eta = new RooHistPdf("pdf_bkg_had_eta","pdf_bkg_had_eta",RooArgSet(*v_eta[arm][charge]),*h_bkg_had_eta);
      RooHistPdf * pdf_bkg_had_dw23 = new RooHistPdf("pdf_bkg_had_dw23","pdf_bkg_had_eta",RooArgSet(*v_dw23[arm][charge]),*h_bkg_had_dw23);
      RooHistPdf * pdf_bkg_mu_eta = new RooHistPdf("pdf_bkg_mu_eta","pdf_bkg_mu_eta",RooArgSet(*v_eta[arm][charge]),*h_bkg_mu_eta);
      RooHistPdf * pdf_bkg_mu_dw23 = new RooHistPdf("pdf_bkg_mu_dw23","pdf_bkg_mu_dw23",RooArgSet(*v_dw23[arm][charge]),*h_bkg_mu_dw23);
      RooHistPdf * pdf_sig_mu_eta = new RooHistPdf("pdf_sig_mu_eta","pdf_sig_mu_eta",RooArgSet(*v_eta[arm][charge]),*h_sig_mu_eta);
      RooHistPdf * pdf_sig_mu_dw23 = new RooHistPdf("pdf_sig_mu_dw23","pdf_sig_mu_dw23",RooArgSet(*v_dw23[arm][charge]),*h_sig_mu_dw23);

      if(save_sbg_hists) {
        sbg_hists_file->cd();
        h_combined_mu_bg_eta[arm][charge]->Write();
        h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Write();
        sprintf(name,"h_sig_mu_eta_a%d_q%d",arm,charge);
        //h_sbg_pdf_sig_mu_eta[arm][charge] = (TH1F*)h_sbg_fit_eta[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Clone(name);
        h_sbg_pdf_sig_mu_eta[arm][charge] = new TH1F(name,name,nhistbins_sbg_eta,1.1,2.6);
        pdf_sig_mu_eta->fillHistogram(h_sbg_pdf_sig_mu_eta[arm][charge],RooArgList(*v_eta[arm][charge]));
        h_sbg_pdf_sig_mu_eta[arm][charge]->Write();
        sprintf(name,"h_sig_mu_dw23_a%d_q%d",arm,charge);
        //h_sbg_pdf_sig_mu_dw23[arm][charge] = (TH1F*)h_sbg_fit_dw23[1][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Clone(name);
        h_sbg_pdf_sig_mu_dw23[arm][charge] = new TH1F(name,name,nhistbins_sbg_dw23,-.1,.1);
        pdf_sig_mu_dw23->fillHistogram(h_sbg_pdf_sig_mu_dw23[arm][charge],RooArgList(*v_dw23[arm][charge]));
        h_sbg_pdf_sig_mu_dw23[arm][charge]->Write();
        sprintf(name,"h_bkg_had_eta_a%d_q%d",arm,charge);
        //h_sbg_pdf_bkg_had_eta[arm][charge] = (TH1F*)h_sbg_fit_eta[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Clone(name);
        h_sbg_pdf_bkg_had_eta[arm][charge] = new TH1F(name,name,nhistbins_sbg_eta,1.1,2.6);
        pdf_bkg_had_eta->fillHistogram(h_sbg_pdf_bkg_had_eta[arm][charge],RooArgList(*v_eta[arm][charge]));
        h_sbg_pdf_bkg_had_eta[arm][charge]->Write();
        sprintf(name,"h_bkg_had_dw23_a%d_q%d",arm,charge);
        //h_sbg_pdf_bkg_had_dw23[arm][charge] = (TH1F*)h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Clone(name);
        h_sbg_pdf_bkg_had_dw23[arm][charge] = new TH1F(name,name,nhistbins_sbg_dw23,-.1,.1);
        pdf_bkg_had_dw23->fillHistogram(h_sbg_pdf_bkg_had_dw23[arm][charge],RooArgList(*v_dw23[arm][charge]));
        h_sbg_pdf_bkg_had_dw23[arm][charge]->Write();
        sprintf(name,"h_bkg_mu_eta_a%d_q%d",arm,charge);
        //h_sbg_pdf_bkg_mu_eta[arm][charge] = (TH1F*)h_combined_mu_bg_eta[arm][charge]->Clone(name);
        h_sbg_pdf_bkg_mu_eta[arm][charge] = new TH1F(name,name,nhistbins_sbg_eta,1.1,2.6);
        pdf_bkg_mu_eta->fillHistogram(h_sbg_pdf_bkg_mu_eta[arm][charge],RooArgList(*v_eta[arm][charge]));
        h_sbg_pdf_bkg_mu_eta[arm][charge]->Write();
        sprintf(name,"h_bkg_mu_dw23_a%d_q%d",arm,charge);
        //h_sbg_pdf_bkg_mu_dw23[arm][charge] = (TH1F*)h_combined_mu_bg_dw23[arm][charge]->Clone(name);
        h_sbg_pdf_bkg_mu_dw23[arm][charge] = new TH1F(name,name,nhistbins_sbg_dw23,-.1,.1);
        pdf_bkg_mu_dw23->fillHistogram(h_sbg_pdf_bkg_mu_dw23[arm][charge],RooArgList(*v_dw23[arm][charge]));
        h_sbg_pdf_bkg_mu_dw23[arm][charge]->Write();
      }

      RooProdPdf * pdf_bkg_had = new RooProdPdf("pdf_bkg_had","pdf_bkg_had",RooArgList(*pdf_bkg_had_eta,*pdf_bkg_had_dw23),0);
      RooProdPdf * pdf_bkg_mu = new RooProdPdf("pdf_bkg_mu","pdf_bkg_mu",RooArgList(*pdf_bkg_mu_eta,*pdf_bkg_mu_dw23),0);
      RooProdPdf * pdf_sig_mu = new RooProdPdf("pdf_sig_mu","pdf_sig_mu",RooArgList(*pdf_sig_mu_eta,*pdf_sig_mu_dw23),0);

      RooRealVar * v_nw = new RooRealVar("v_nw", "W number", 300, 1.0, 3000);
      RooRealVar * v_nbg = new RooRealVar("v_nbg", "BG number", 1500, 1.0, 3000);
      RooRealVar * v_nhad = new RooRealVar("v_nhad", "fake number", 1, 1.0, 3000);

      RooAddPdf *component_model;
      component_model= new RooAddPdf("component_model", "model", RooArgList(*pdf_sig_mu, *pdf_bkg_mu, *pdf_bkg_had), RooArgList(*v_nw, *v_nbg, *v_nhad));
      *v_nbg = mu_bg_yield;
      //*v_nbg = 30;

      printf("\n\nmuonic bg yields: %f\n\n",mu_bg_yield);
      v_nbg->setConstant(true);

      component_model->fitTo(*data[arm][charge],
          Extended(kTRUE),
          Minos(kTRUE),
          NumCPU(2),
          SumW2Error(kTRUE),
          PrintLevel(-1));

      sig_w_count[arm][charge] = v_nw->getVal();
      bkg_mu_count[arm][charge] =  v_nbg->getVal();
      bkg_had_count[arm][charge] = v_nhad->getVal();

      sig_w_err_lo[arm][charge] = v_nw->getErrorLo();
      sig_w_err_hi[arm][charge] = v_nw->getErrorHi();
      bkg_mu_err_lo[arm][charge] =  v_nbg->getErrorLo();
      bkg_mu_err_hi[arm][charge] =  v_nbg->getErrorHi();
      bkg_had_err_lo[arm][charge] = v_nhad->getErrorLo();
      bkg_had_err_hi[arm][charge] = v_nhad->getErrorHi();

      //do range fraction
      v_eta[arm][charge]->setRange("testrange",distmin[5],distmax[5]);
      v_dw23[arm][charge]->setRange("testrange",distmin[14],distmax[14]);//(-0.05+charge*0.04),(0.01+charge*0.04));
      RooAbsReal* sig_mu_range_fraction = pdf_bkg_had->createIntegral(RooArgSet(*v_eta[arm][charge],*v_dw23[arm][charge]),NormSet(RooArgSet(*v_eta[arm][charge],*v_dw23[arm][charge])),Range("testrange"));
      std::cout << "\n\ntest: " << " " << sig_mu_range_fraction->getVal() << std::endl;

      
      c_sbg_fits_dw23->cd(charge+arm*2+1);
      frame_dw23[arm][charge] = v_dw23[arm][charge]->frame();
      data[arm][charge]->plotOn(frame_dw23[arm][charge], Binning(nhistbins_sbg_dw23), DataError(RooAbsData::SumW2) );
      component_model->plotOn(frame_dw23[arm][charge],LineWidth(1),LineColor(kBlack));
      component_model->plotOn(frame_dw23[arm][charge],Components(*pdf_sig_mu),LineWidth(1),LineStyle(kDashed),LineColor(kRed));
      component_model->plotOn(frame_dw23[arm][charge],Components(*pdf_bkg_mu),LineWidth(1),LineStyle(kDashed),LineColor(kGreen));
      component_model->plotOn(frame_dw23[arm][charge],Components(*pdf_bkg_had),LineWidth(1),LineStyle(kDashed),LineColor(kBlue));
      sprintf(name,"S/BG Fit  dw23  %c%c",arm?'N':'S',charge?'+':'-');
      frame_dw23[arm][charge]->SetTitle(name);
      frame_dw23[arm][charge]->Draw();
      
      c_sbg_fits_eta->cd(charge+arm*2+1);
      frame_eta[arm][charge] = v_eta[arm][charge]->frame(1.1, 2.6);
      data[arm][charge]->plotOn(frame_eta[arm][charge], Binning(nhistbins_sbg_eta), DataError(RooAbsData::SumW2) );
      component_model->plotOn(frame_eta[arm][charge],LineWidth(1),LineColor(kBlack));
      component_model->plotOn(frame_eta[arm][charge],Components(*pdf_sig_mu),LineWidth(1),LineStyle(kDashed),LineColor(kRed));
      component_model->plotOn(frame_eta[arm][charge],Components(*pdf_bkg_mu),LineWidth(1),LineStyle(kDashed),LineColor(kGreen));
      component_model->plotOn(frame_eta[arm][charge],Components(*pdf_bkg_had),LineWidth(1),LineStyle(kDashed),LineColor(kBlue));
      sprintf(name,"S/BG Fit  eta  %c%c",arm?'N':'S',charge?'+':'-');
      frame_eta[arm][charge]->SetTitle(name);
      frame_eta[arm][charge]->Draw(); 
    }
  }

  if(save_sbg_hists)
    sbg_hists_file->Close();

  printf("\nMuon Bkg Scaling Factors (including eff):");
  for(int i=0; i<9; i++)
    printf("\n%f",mu_bg_scaling_factor[i]);

  printf("\n\n");
  
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      printf("\n-- Mu Bg Yields Arm%d Charge %d --\n",arm,charge);
      for(int process=0; process<9; process++) {
        printf("%f\n",mu_bg_yields[arm][charge][process]);
        //printf("%s %f\n",mu_bg_labels[process],mu_bg_yields[arm][charge][process]);
      }
      printf("\n");
      for(int process=0; process<9; process++) {
        printf("%f\n",mu_bg_yields[arm][charge][process]/mu_bg_scaling_factor[process]);
        //printf("%s %f\n",mu_bg_labels[process],mu_bg_yields[arm][charge][process]);
      }
    }
  }

  ofstream smout;
  smout.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit/output/sig_mu_yields.txt");
  smout << "raw and trigger efficiency scaled yields for signal muon simulation events in the target wness region (wness>.92)" << endl;
  smout << "         arm      charge      rawyld  trigscayld" << endl;
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      smout << std::setw(12) << arm << std::setw(12) << charge << std::setw(12) << raw_target_wness_count[1][arm][charge] << std::setw(12) << trigeff_target_wness_count[1][arm][charge] << endl;
    }
  }
  smout.close();
  
  ofstream bhout;
  bhout.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit/output/bkg_had_yields.txt");
  bhout << "raw and trigger efficiency scaled yields for background hadron events from data in the background wness region (0.1<wness<0.7)" << endl;
  bhout << "         arm      charge      rawyld  trigscayld" << endl;
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      bhout << std::setw(12) << arm << std::setw(12) << charge << std::setw(12) << raw_target_wness_count[0][arm][charge] << std::setw(12) << trigeff_target_wness_count[0][arm][charge] << endl;
    }
  }
  bhout.close();

  ofstream bmout;
  bmout.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit/output/bkg_mu_scaling_yields.txt");
  bmout << "raw,trigeff scaled, and total scaled yields (and scale factors) for bkg muon events from simulation in the target wness region (wness>.92)" << endl;
  bmout << "     process         arm      charge      rawyld  trigscayld    #genevts    xsection   dimuonfac   totsclfac   totscayld" << endl;
  for(int process=0; process<9; process++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        bmout << std::setw(12) << mu_bg_labels[process] << std::setw(12) << arm << std::setw(12) << charge << std::setw(12) << raw_target_wness_count[process+2][arm][charge] << std::setw(12) << mu_bg_yields[arm][charge][process]/mu_bg_scaling_factor[process] << std::setw(12) << mu_bg_generated_events[process] << std::setw(12) << mu_bg_cross_section[process] << std::setw(12) << mu_bg_dimuon_factor[process] << std::setw(12) << mu_bg_scaling_factor[process] << std::setw(12) << mu_bg_yields[arm][charge][process] << endl;
      }
    }
  }
  bmout.close();

  ofstream sbg_outfile;
  sbg_outfile.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit/output/sbg_values.txt");
  sbg_outfile << "arm0charg0 - arm0charge1 - arm1charge0 - arm1charge1\n";
  sbg_outfile << "sbg_ratio\twsig_yld\twsig_errhi\twsig_errlo\thad_yld\thad_errhi\thad_errlo\tmubkg_yld\n";


  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      printf("\n-- Arm%d Charge %d --\n",arm,charge);
      printf("Data entries:       %.2f\n", data_count[arm][charge]);
      printf("Summed fit entries: %.2f\n", sig_w_count[arm][charge]+bkg_mu_count[arm][charge]+bkg_had_count[arm][charge]);
      printf("Signal fit:         %.2f (+%.2f %.2f)\n", sig_w_count[arm][charge],sig_w_err_hi[arm][charge],sig_w_err_lo[arm][charge]);
      printf("Had bkg fit:        %.2f (+%.2f %.2f)\n", bkg_had_count[arm][charge],bkg_had_err_hi[arm][charge],bkg_had_err_lo[arm][charge]);
      printf("Fixed mu bkg:       %.2f\n", bkg_mu_count[arm][charge]);
      printf("s/bg ratio:         %.3f\n", sig_w_count[arm][charge]/(bkg_mu_count[arm][charge]+bkg_had_count[arm][charge]));

      sbg_outfile << sig_w_count[arm][charge]/(bkg_mu_count[arm][charge]+bkg_had_count[arm][charge]) << "\t" << 
        sig_w_count[arm][charge] << "\t" << sig_w_err_hi[arm][charge] << "\t" << sig_w_err_lo[arm][charge] << "\t" << 
        bkg_had_count[arm][charge] << "\t" << bkg_had_err_hi[arm][charge] << "\t" << bkg_had_err_lo[arm][charge]  << "\t" << 
        bkg_mu_count[arm][charge] << endl;

      
      /*double temp_integral;

      c_sbg_fit_dw23->cd(charge+2*arm+1);
      temp_integral = h_dw23_data[arm][charge]->Integral("width");
      if(temp_integral > 0)
        h_dw23_data[arm][charge]->Scale(h_dw23_data[arm][charge]->GetEntries()/temp_integral);
      h_dw23_data[arm][charge]->Draw("L");
      h_dw23_data[arm][charge]->SetLineColor(1);
      
      temp_integral = h_sbg_fit_dw23[0][arm][charge]->Integral("width");
      if(temp_integral > 0)
        h_sbg_fit_dw23[0][arm][charge]->Scale(bkg_had_count[arm][charge]/temp_integral);
      h_sbg_fit_dw23[0][arm][charge]->Draw("same");
      h_sbg_fit_dw23[0][arm][charge]->SetLineColor(2);
      
      temp_integral = h_sbg_fit_dw23[1][arm][charge]->Integral("width");
      if(temp_integral > 0)
        h_sbg_fit_dw23[1][arm][charge]->Scale(sig_w_count[arm][charge]/temp_integral);
      h_sbg_fit_dw23[1][arm][charge]->Draw("same");
      h_sbg_fit_dw23[1][arm][charge]->SetLineColor(3);
      
      temp_integral = h_combined_mu_bg_dw23[arm][charge]->Integral("width");
      if(temp_integral > 0)
        h_combined_mu_bg_dw23[arm][charge]->Scale(bkg_mu_count[arm][charge]/temp_integral);
      h_combined_mu_bg_dw23[arm][charge]->Draw("same");
      h_combined_mu_bg_dw23[arm][charge]->SetLineColor(4);
      
      c_sbg_fit_eta->cd(charge+2*arm+1);
      temp_integral = h_eta_data[arm][charge]->Integral();
      if(temp_integral > 0)
        h_eta_data[arm][charge]->Scale(h_eta_data[arm][charge]->GetEntries()/temp_integral);
      h_eta_data[arm][charge]->Draw("L");
      h_eta_data[arm][charge]->SetLineColor(1);
      
      temp_integral = h_sbg_fit_eta[0][arm][charge]->Integral();
      if(temp_integral > 0)
        h_sbg_fit_eta[0][arm][charge]->Scale(bkg_had_count[arm][charge]/temp_integral);
      h_sbg_fit_eta[0][arm][charge]->Draw("same");
      h_sbg_fit_eta[0][arm][charge]->SetLineColor(2);
      
      temp_integral = h_sbg_fit_eta[1][arm][charge]->Integral();
      if(temp_integral > 0)
        h_sbg_fit_eta[1][arm][charge]->Scale(sig_w_count[arm][charge]/temp_integral);
      h_sbg_fit_eta[1][arm][charge]->Draw("same");
      h_sbg_fit_eta[1][arm][charge]->SetLineColor(3);
      
      temp_integral = h_combined_mu_bg_eta[arm][charge]->Integral();
      if(temp_integral > 0)
        h_combined_mu_bg_eta[arm][charge]->Scale(bkg_mu_count[arm][charge]/temp_integral);
      h_combined_mu_bg_eta[arm][charge]->Draw("same");
      h_combined_mu_bg_eta[arm][charge]->SetLineColor(4);
      
      */
    }
  }
  
  //mu_bg_hists_file->Close();
  sbg_outfile.close();
  printf("testend\n");  
}



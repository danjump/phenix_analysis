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

  bool save_fit_hists=false;
  bool read_fit_hists=false;
  bool do_extrap_plots=false;
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
  
  TFile * hist_infile;
  if(read_fit_hists) {
    hist_infile = new TFile("/direct/phenix+spin2/beaumim/WAnalysisRun13/AnalysisCode/macros/rootTrees/Run13_Data_Plots.root");
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        sprintf(name,"dw23VsWness_Arm%d_Charge%d_data",arm,charge);
        h2_dw23_vs_wness[0][arm][charge] = (TH2F*)hist_infile->Get(name);
        sprintf(name,"Wness_Fit_Seed_Arm%d_Charge%d_data",arm,charge);
        h_wness[0][arm][charge] = (TH1F*)hist_infile->Get(name);
      }
    }
  }
  //DW23 EXTRAPOLATION FITTING

  //2d dw23 vs wness fitting
  //double count_region_2d[2][2][6];
  //double count_region_1d[2][2][6];
  int lower_w_bin, upper_w_bin;
  
  /*
  for(int source=0; source < 1; source++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        for(int i=0; i<4; i++) {
          lower_w_bin = (2*i+1)*nhistbins/10+1;
          upper_w_bin = lower_w_bin+2*nhistbins/10-1;
          count_region_2d[arm][charge][i+1] = h2_dw23_vs_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin,1,nhistbins); 
          count_region_1d[arm][charge][i+1] = h_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin); 
        }
        
        lower_w_bin = (2*4+1)*nhistbins/10+1;
        upper_w_bin = lower_w_bin+1*nhistbins/10-1;
        count_region_2d[arm][charge][5] = h2_dw23_vs_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin,1,nhistbins); 
        count_region_1d[arm][charge][5] = h_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin); 

        lower_w_bin = 1;
        upper_w_bin = lower_w_bin+1*nhistbins/10-1;
        count_region_2d[arm][charge][0] = h2_dw23_vs_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin,1,nhistbins);
        count_region_1d[arm][charge][0] = h_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin); 
        
        double sum = count_region_2d[arm][charge][0]+count_region_2d[arm][charge][1]+count_region_2d[arm][charge][2]+
                     count_region_2d[arm][charge][3]+count_region_2d[arm][charge][4]+count_region_2d[arm][charge][5];
        printf("2D STATS: total:%6.0f .0-.1:%6.0f .1-.3:%4.0f .3-.5:%4.0f .5-.7:%4.0f .7-.9:%4.0f .9-1:%4.0f\n",
            sum,count_region_2d[arm][charge][0],count_region_2d[arm][charge][1],count_region_2d[arm][charge][2],
            count_region_2d[arm][charge][3],count_region_2d[arm][charge][4],count_region_2d[arm][charge][5]);
        
        sum = count_region_1d[arm][charge][0]+count_region_1d[arm][charge][1]+count_region_1d[arm][charge][2]+
                     count_region_1d[arm][charge][3]+count_region_1d[arm][charge][4]+count_region_1d[arm][charge][5];
        printf("1D STATS: total:%6.0f .0-.1:%6.0f .1-.3:%4.0f .3-.5:%4.0f .5-.7:%4.0f .7-.9:%4.0f .9-1:%4.0f\n",
            sum,count_region_1d[arm][charge][0],count_region_1d[arm][charge][1],count_region_1d[arm][charge][2],
            count_region_1d[arm][charge][3],count_region_1d[arm][charge][4],count_region_1d[arm][charge][5]);
        printf("under:%d in:%d over:%d\n",under_count[source][arm][charge]++,in_count[source][arm][charge]++,over_count[source][arm][charge]++);
      }
    }
  }*/
 
  int n_pars_2d = 13;
  
  double fit_chi2ndf[3][2][2];
  double mean_linear_pars[3][2][2][2];
  double sigma1_linear_pars[3][2][2][2];
  double sigma2_linear_pars[3][2][2][2];
  double factor_linear_pars[3][2][2][2];
  double poly_pars[3][2][2][5];

  double mean_linear_pars_err[3][2][2][2];
  double sigma1_linear_pars_err[3][2][2][2];
  double sigma2_linear_pars_err[3][2][2][2];
  double factor_linear_pars_err[3][2][2][2];
  double poly_pars_err[3][2][2][5];
  
  double dw23_vs_wness_integral[3][2][2];

  TCanvas * c_wness_fits; 
  if(do_extrap_plots) {
    c_wness_fits = new TCanvas("c_wness_fits","c_wness_fits",1000,1000);
    c_wness_fits->Divide(2,2);
  }
  TF1 *wness_pol[3][2][2];
  

  TFile * f_fit_hists;
  if(save_fit_hists) {
    f_fit_hists = new TFile("temp_fit_hists.root","RECREATE");
  }

  for(int source=0; source < 1; source++) {//start loop for 2d dw 23 extrapolation fit
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        h2_dw23_vs_wness[source][arm][charge]->Sumw2();
        dw23_vs_wness_integral[source][arm][charge] = h2_dw23_vs_wness[source][arm][charge]->Integral("width");
        if(dw23_vs_wness_integral[source][arm][charge]!=0)
          h2_dw23_vs_wness[source][arm][charge]->Scale(1.0/dw23_vs_wness_integral[source][arm][charge]);
        
        if(save_fit_hists) {
          f_fit_hists->cd();
          h2_dw23_vs_wness[source][arm][charge]->Write();
        }
        TF2 * f2_double_gaus_2d = new TF2("f2_double_gaus_2d",double_gaus_2d,0.1,0.9,-.3,.3,n_pars_2d);
        
        //linear gaussian parameter seeds (from hide thesis mostly)
        /*f2_double_gaus_2d->SetParameter(0,0.01);        // mean(offset) constant
        f2_double_gaus_2d->SetParameter(1,0.01); // mean(offset) slope
        
        f2_double_gaus_2d->SetParameter(2,0.1);               // sigma1 constant
        f2_double_gaus_2d->SetParameter(3,-0.1);             // sigma1 slope
        
        f2_double_gaus_2d->SetParameter(4,0.02);              // sigma2 constant
        f2_double_gaus_2d->SetParameter(5,-0.1);             // sigma2 slope
        
        f2_double_gaus_2d->SetParameter(6,1);                 // constant(factor) constant
        f2_double_gaus_2d->SetParameter(7,-1);                 // constant(factor) slope*/

        double dw23_average = (h_dw23[source][arm][charge][1]->GetMean()+
            h_dw23[source][arm][charge][2]->GetMean()+
            h_dw23[source][arm][charge][3]->GetMean())/3.0;

        f2_double_gaus_2d->SetParameter(0,dw23_average);        // mean(offset) constant
        f2_double_gaus_2d->SetParameter(1,(charge?-1:1)*0.005); // mean(offset) slope
        
        f2_double_gaus_2d->SetParameter(2,0.1);               // sigma1 constant
        f2_double_gaus_2d->SetParameter(3,-0.06);             // sigma1 slope
        
        f2_double_gaus_2d->SetParameter(4,0.04);              // sigma2 constant
        f2_double_gaus_2d->SetParameter(5,-0.01);             // sigma2 slope
        
        f2_double_gaus_2d->SetParameter(6,1);                 // constant(factor) constant
        f2_double_gaus_2d->SetParameter(7,1);                 // constant(factor) slope

        
        //parameter limits also estimated from hides thesis for an initial fit
        /*f2_double_gaus_2d->SetParLimits(0,dw23_average*0.5,dw23_average*2);        // mean(offset) constant
        f2_double_gaus_2d->SetParLimits(1,(charge?-1:0)*0.01,(charge?0:1)*0.01);   // mean(offset) slope

        f2_double_gaus_2d->SetParLimits(2,0.05,0.8);               // sigma1 constant
        f2_double_gaus_2d->SetParLimits(3,-0.5,0);                 // sigma1 slope
        
        f2_double_gaus_2d->SetParLimits(4,0.008,0.15);             // sigma2 constant
        f2_double_gaus_2d->SetParLimits(5,-0.001,0);               // sigma2 slope
        
        f2_double_gaus_2d->SetParLimits(6,-30,30);                 // constant(factor) constant
        f2_double_gaus_2d->SetParLimits(7,-30,30);*/                 // constant(factor) slope
        
        //fit wness distribution with polynomial to get parameter seeds
        if(do_extrap_plots) {
          c_wness_fits->cd(arm+2*charge+1);
          gPad->SetLogy();
        }
        h_wness[source][arm][charge]->Sumw2();
        double integral = h_wness[source][arm][charge]->Integral("width");
        if(integral != 0)
          h_wness[source][arm][charge]->Scale(1/integral);
        if(do_extrap_plots)
          h_wness[source][arm][charge]->Draw("eP");
        if(save_fit_hists) {
          h_wness[source][arm][charge]->Write();
        }
        sprintf(name,"wness_pol_source%d_charge%d_arm%d",source,arm,charge);
        wness_pol[source][arm][charge] = new TF1(name,"pol4(0)",0.1,0.9);
        h_wness[source][arm][charge]->Fit(name,"QLR");
        TF1 *wness_pol_tmp = new TF1("wness_pol_tmp","pol4(0)",0,1);
        for(int i=0; i<5; i++)
          wness_pol_tmp->SetParameter(i,wness_pol[source][arm][charge]->GetParameter(i));
        wness_pol_tmp->SetLineColor(1);
        wness_pol_tmp->SetLineStyle(2);
        if(do_extrap_plots) {
          wness_pol_tmp->DrawCopy("same");
          if(arm==1 && charge==1)
            c_wness_fits->SaveAs("/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/sbg_fit/temp/wness_fits.png");
        }
        
        f2_double_gaus_2d->FixParameter(8,wness_pol[source][arm][charge]->GetParameter(0));
        f2_double_gaus_2d->FixParameter(9,wness_pol[source][arm][charge]->GetParameter(1));
        f2_double_gaus_2d->FixParameter(10,wness_pol[source][arm][charge]->GetParameter(2));
        f2_double_gaus_2d->FixParameter(11,wness_pol[source][arm][charge]->GetParameter(3));
        f2_double_gaus_2d->FixParameter(12,wness_pol[source][arm][charge]->GetParameter(4));
        
        f2_double_gaus_2d->SetParLimits(2,.09,1);
        f2_double_gaus_2d->SetParLimits(4,0,.09);
        
        h2_dw23_vs_wness[source][arm][charge]->Fit("f2_double_gaus_2d","QLRN");

        fit_chi2ndf[source][arm][charge] = f2_double_gaus_2d->GetChisquare() / f2_double_gaus_2d->GetNDF();
        
        mean_linear_pars[source][arm][charge][0] = f2_double_gaus_2d->GetParameter(0);
        mean_linear_pars[source][arm][charge][1] = f2_double_gaus_2d->GetParameter(1);
        sigma1_linear_pars[source][arm][charge][0] = f2_double_gaus_2d->GetParameter(2);
        sigma1_linear_pars[source][arm][charge][1] = f2_double_gaus_2d->GetParameter(3);
        sigma2_linear_pars[source][arm][charge][0] = f2_double_gaus_2d->GetParameter(4);
        sigma2_linear_pars[source][arm][charge][1] = f2_double_gaus_2d->GetParameter(5);
        factor_linear_pars[source][arm][charge][0] = f2_double_gaus_2d->GetParameter(6);
        factor_linear_pars[source][arm][charge][1] = f2_double_gaus_2d->GetParameter(7);
        poly_pars[source][arm][charge][0] = f2_double_gaus_2d->GetParameter(8);
        poly_pars[source][arm][charge][1] = f2_double_gaus_2d->GetParameter(9);
        poly_pars[source][arm][charge][2] = f2_double_gaus_2d->GetParameter(10);
        poly_pars[source][arm][charge][3] = f2_double_gaus_2d->GetParameter(11);
        poly_pars[source][arm][charge][4] = f2_double_gaus_2d->GetParameter(12);
        
        mean_linear_pars_err[source][arm][charge][0] = f2_double_gaus_2d->GetParError(0);
        mean_linear_pars_err[source][arm][charge][1] = f2_double_gaus_2d->GetParError(1);
        sigma1_linear_pars_err[source][arm][charge][0] = f2_double_gaus_2d->GetParError(2);
        sigma1_linear_pars_err[source][arm][charge][1] = f2_double_gaus_2d->GetParError(3);
        sigma2_linear_pars_err[source][arm][charge][0] = f2_double_gaus_2d->GetParError(4);
        sigma2_linear_pars_err[source][arm][charge][1] = f2_double_gaus_2d->GetParError(5);
        factor_linear_pars_err[source][arm][charge][0] = f2_double_gaus_2d->GetParError(6);
        factor_linear_pars_err[source][arm][charge][1] = f2_double_gaus_2d->GetParError(7);
        poly_pars_err[source][arm][charge][0] = f2_double_gaus_2d->GetParError(8);
        poly_pars_err[source][arm][charge][1] = f2_double_gaus_2d->GetParError(9);
        poly_pars_err[source][arm][charge][2] = f2_double_gaus_2d->GetParError(10);
        poly_pars_err[source][arm][charge][3] = f2_double_gaus_2d->GetParError(11);
        poly_pars_err[source][arm][charge][4] = f2_double_gaus_2d->GetParError(12);
        
      }
    }
  }

  printf("test2\n");

  //combine and scale dw23 histograms for display
  TH1F *h_dw23_sums[3][2][2][6];
  TH1F *h_wness_sums[3][2][2][5];
  double integral_1d, integral_region_2d[2][2][6];
  
  for(int source=0; source < 1; source++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        for(int i=0; i<4; i++) {
          if(read_fit_hists) {
            sprintf(name,"dw23Wness_WnessBin_%d_Arm%d_Charge%d_data",i+1,arm,charge);
            h_dw23_sums[0][arm][charge][i] = (TH1F*)hist_infile->Get(name);
          } else {

            sprintf(name,"h_dw23_sums_source%d_arm%d_charge%d_range%d",source,arm,charge,i);
            h_dw23_sums[source][arm][charge][i] = new TH1F(name,name,nhistbins_2ddw23,distmin[14],distmax[14]);

            h_dw23_sums[source][arm][charge][i]->Add(h_dw23[source][arm][charge][(2*i)+1],h_dw23[source][arm][charge][(2*i)+2]);
          }
          sprintf(name,"h_wness_sums_source%d_arm%d_charge%d_range%d",source,arm,charge,i);
          h_wness_sums[source][arm][charge][i] = new TH1F(name,name,nhistbins_2dwness,0,1);
          
          h_wness_sums[source][arm][charge][i]->Add(h_wness_sections[source][arm][charge][(2*i)+1],h_wness_sections[source][arm][charge][(2*i)+2]);
          
          //scale the histogram to match the function slice that will be drawn
          integral_1d = h_dw23_sums[source][arm][charge][i]->Integral("width");
          lower_w_bin = (2*i+1)*nhistbins_2dwness/10+1;
          upper_w_bin = lower_w_bin+2*nhistbins_2dwness/10-1;
          integral_region_2d[arm][charge][i+1] = h2_dw23_vs_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin,1,nhistbins_2ddw23,"width"); 
          h_dw23_sums[source][arm][charge][i]->Sumw2();
          if(integral_1d != 0) {
            //set the integral of the 1d curve equal to the integral of the 2d section/wness width
            h_dw23_sums[source][arm][charge][i]->Scale((integral_region_2d[arm][charge][i+1]/0.2)/integral_1d); 
          }
          if(save_fit_hists) {
            h_dw23_sums[source][arm][charge][i]->Write();
          }
        }
        if(read_fit_hists) {
          sprintf(name,"dw23Wness_WnessBin_%d_Arm%d_Charge%d_data",4+1,arm,charge);
          h_dw23_sums[0][arm][charge][4] = (TH1F*)hist_infile->Get(name);
        } else {
          sprintf(name,"h_dw23_sums_source%d_arm%d_charge%d_range%d",source,arm,charge,4);
          h_dw23_sums[source][arm][charge][4] = (TH1F*)h_dw23[source][arm][charge][9]->Clone(name);
        }
        
        sprintf(name,"h_wness_sums_source%d_arm%d_charge%d_range%d",source,arm,charge,4);
        h_wness_sums[source][arm][charge][4] = (TH1F*)h_wness_sections[source][arm][charge][9]->Clone(name);

        //scale the histogram to match the function slice that will be drawn
        integral_1d = h_dw23_sums[source][arm][charge][4]->Integral("width");
        lower_w_bin = (2*4+1)*nhistbins_2dwness/10+1;
        upper_w_bin = lower_w_bin+1*nhistbins_2dwness/10-1;
        integral_region_2d[arm][charge][5] = h2_dw23_vs_wness[source][arm][charge]->Integral(lower_w_bin,upper_w_bin,1,nhistbins_2ddw23,"width"); 
        h_dw23_sums[source][arm][charge][4]->Sumw2();
        if(integral_1d != 0) {
          //set the integral of the 1d curve equal to the integral of the 2d section/wness width
          h_dw23_sums[source][arm][charge][4]->Scale((integral_region_2d[arm][charge][5]/0.1)/integral_1d); 
        }

        if(save_fit_hists) {
          h_dw23_sums[source][arm][charge][4]->Write();
        }
      }
    }
  }

  
  if(do_extrap_plots) {
    TCanvas * c_h2_dw23_vs_wness[3][2][2];
    TCanvas * c_dw23_fits_summary;
    c_dw23_fits_summary = new TCanvas("c_dw23_fits_summary","c_dw23_fits_summary",1600,1600);
    c_dw23_fits_summary->Divide(4,5);

    char * arm_label[2] = {(char*)"South",(char*)"North"};
    char * charge_label[2] = {(char*)"#mu^{-}",(char*)"#mu^{+}"};
    for(int source=0; source < 1; source++) {//display loop for 2d fitting
      for(int arm=0; arm<2; arm++) {
        for(int charge=0; charge<2; charge++) {
          sprintf(name,"c_h2_dw23_vs_wness_source%d_arm%d_charge%d",source,arm,charge);
          c_h2_dw23_vs_wness[source][arm][charge] = new TCanvas(name,name,1800,1000);
          c_h2_dw23_vs_wness[source][arm][charge]->Divide(3,3);

          /*printf("\n\n--- Arm%d Charge%d ---\n",arm,charge);
            double sum = count_region_2d[arm][charge][0]+count_region_2d[arm][charge][1]+
            count_region_2d[arm][charge][2]+count_region_2d[arm][charge][3]+count_region_2d[arm][charge][4]+count_region_2d[arm][charge][5];
            printf("Statistics: total:%.0f .0-.1:%.0f .1-.3:%.0f .3-.5:%.0f .5-.7:%.0f .7-.9:%.0f .9-1:%.0f\n",sum,count_region_2d[arm][charge][0],count_region_2d[arm][charge][1],
            count_region_2d[arm][charge][2],count_region_2d[arm][charge][3],count_region_2d[arm][charge][4],count_region_2d[arm][charge][5]);
            printf("offset: const=%9f+-%9f  slope=%9f+-%9f\n",
            mean_linear_pars[source][arm][charge][0],mean_linear_pars_err[source][arm][charge][0],
            mean_linear_pars[source][arm][charge][1],mean_linear_pars_err[source][arm][charge][1]);
            printf("Sigma1: const=%9f+-%9f  slope=%9f+-%9f\n",
            sigma1_linear_pars[source][arm][charge][0],sigma1_linear_pars_err[source][arm][charge][0],
            sigma1_linear_pars[source][arm][charge][1],sigma1_linear_pars_err[source][arm][charge][1]);
            printf("Sigma2: const=%9f+-%9f  slope=%9f+-%9f\n",
            sigma2_linear_pars[source][arm][charge][0],sigma2_linear_pars_err[source][arm][charge][0],
            sigma2_linear_pars[source][arm][charge][1],sigma2_linear_pars_err[source][arm][charge][1]);
            printf("Factor: const=%9f+-%9f  slope=%9f+-%9f\n",
            factor_linear_pars[source][arm][charge][0],factor_linear_pars_err[source][arm][charge][0],
            factor_linear_pars[source][arm][charge][1],factor_linear_pars_err[source][arm][charge][1]);*/

          printf("\n--- Arm%d Charge%d ---\n",arm,charge);
          printf("Fit Chi2/ndf:   %9f\n",fit_chi2ndf[source][arm][charge]);
          printf("offset const:   %9f   +-  (%9f)\n",
              mean_linear_pars[source][arm][charge][0],mean_linear_pars_err[source][arm][charge][0]);
          printf("offset slope:   %9f   +-  (%9f)\n",
              mean_linear_pars[source][arm][charge][1],mean_linear_pars_err[source][arm][charge][1]);
          printf("Sigma1 const:   %9f   +-  (%9f)\n",
              sigma1_linear_pars[source][arm][charge][0],sigma1_linear_pars_err[source][arm][charge][0]);
          printf("Sigma1 slope:   %9f   +-  (%9f)\n",
              sigma1_linear_pars[source][arm][charge][1],sigma1_linear_pars_err[source][arm][charge][1]);
          printf("Sigma2 const:   %9f   +-  (%9f)\n",
              sigma2_linear_pars[source][arm][charge][0],sigma2_linear_pars_err[source][arm][charge][0]);
          printf("Sigma2 slope:   %9f   +-  (%9f)\n",
              sigma2_linear_pars[source][arm][charge][1],sigma2_linear_pars_err[source][arm][charge][1]);
          printf("Factor const:   %9f   +-  (%9f)\n",
              factor_linear_pars[source][arm][charge][0],factor_linear_pars_err[source][arm][charge][0]);
          printf("Factor slope:   %9f   +-  (%9f)\n",
              factor_linear_pars[source][arm][charge][1],factor_linear_pars_err[source][arm][charge][1]);
          printf("Poly Par[0]:    %9f   +-  (%9f)\n",
              poly_pars[source][arm][charge][0],poly_pars_err[source][arm][charge][0]);
          printf("Poly Par[1]:    %9f   +-  (%9f)\n",
              poly_pars[source][arm][charge][1],poly_pars_err[source][arm][charge][1]);
          printf("Poly Par[2]:    %9f   +-  (%9f)\n",
              poly_pars[source][arm][charge][2],poly_pars_err[source][arm][charge][2]);
          printf("Poly Par[3]:    %9f   +-  (%9f)\n",
              poly_pars[source][arm][charge][3],poly_pars_err[source][arm][charge][3]);
          printf("Poly Par[4]:    %9f   +-  (%9f)\n",
              poly_pars[source][arm][charge][4],poly_pars_err[source][arm][charge][4]);

          TF1 *f_2d_gaus_slice[5];
          double axis_max = 0;

          for(int region=0; region<4; region++) {
            sprintf(name,"f_2d_gaus_slice_%d",region);
            double wness_temp = h_wness_sums[source][arm][charge][region]->GetMean();
            f_2d_gaus_slice[region] = new TF1(name,double_gaus_slice,distmin[14],distmax[14],5);
            f_2d_gaus_slice[region]->SetParameter(0,mean_linear_pars[source][arm][charge][0]+
                mean_linear_pars[source][arm][charge][1]*wness_temp);
            f_2d_gaus_slice[region]->SetParameter(1,sigma1_linear_pars[source][arm][charge][0]+
                sigma1_linear_pars[source][arm][charge][1]*wness_temp);
            f_2d_gaus_slice[region]->SetParameter(2,sigma2_linear_pars[source][arm][charge][0]+
                sigma2_linear_pars[source][arm][charge][1]*wness_temp);
            f_2d_gaus_slice[region]->SetParameter(3,factor_linear_pars[source][arm][charge][0]+
                factor_linear_pars[source][arm][charge][1]*wness_temp);
            f_2d_gaus_slice[region]->SetParameter(4,
                poly_pars[source][arm][charge][0]+
                poly_pars[source][arm][charge][1]*(wness_temp)+
                poly_pars[source][arm][charge][2]*pow(wness_temp,2)+
                poly_pars[source][arm][charge][3]*pow(wness_temp,3)+
                poly_pars[source][arm][charge][4]*pow(wness_temp,4));

            c_h2_dw23_vs_wness[source][arm][charge]->cd(region+1);
            if(region==0) 
              axis_max = h_dw23_sums[source][arm][charge][region]->GetMaximum()*1.15;
            h_dw23_sums[source][arm][charge][region]->GetYaxis()->SetRangeUser(0,axis_max);
            sprintf(name,"%s %s dw_{23} (%3.1f < W_{ness} < %3.1f)",arm_label[arm],charge_label[charge],((region*2)+1)/10.0,((region*2)+3)/10.0);
            h_dw23_sums[source][arm][charge][region]->SetTitle(name);
            h_dw23_sums[source][arm][charge][region]->SetTitleSize(.09);
            h_dw23_sums[source][arm][charge][region]->GetXaxis()->SetLabelSize(.06);
            h_dw23_sums[source][arm][charge][region]->GetXaxis()->SetTitle("dw_{23}");
            h_dw23_sums[source][arm][charge][region]->GetXaxis()->SetTitleSize(.06);
            h_dw23_sums[source][arm][charge][region]->GetYaxis()->SetLabelSize(.06);
            h_dw23_sums[source][arm][charge][region]->GetYaxis()->SetTitle("Yield (2D Normalized)");
            h_dw23_sums[source][arm][charge][region]->GetYaxis()->SetTitleSize(.06);
            h_dw23_sums[source][arm][charge][region]->Draw("ep");
            //plot lower, middle, and upper slices of 2d function for this region
            f_2d_gaus_slice[region]->SetLineColor(2);
            f_2d_gaus_slice[region]->Draw("same");

            c_dw23_fits_summary->cd((arm*2+charge)+region*4+1);
            h_dw23_sums[source][arm][charge][region]->Draw("ep");
            f_2d_gaus_slice[region]->Draw("same");
          }

          c_h2_dw23_vs_wness[source][arm][charge]->cd(5);
          double wness_temp = h_wness_sums[source][arm][charge][4]->GetMean();
          TF1 *f_dw23_proj = new TF1("f_dw23_proj",double_gaus_slice,distmin[14],distmax[14],5);
          f_dw23_proj->SetParameter(0,mean_linear_pars[source][arm][charge][0]+mean_linear_pars[source][arm][charge][1]*wness_temp);
          f_dw23_proj->SetParameter(1,sigma1_linear_pars[source][arm][charge][0]+sigma1_linear_pars[source][arm][charge][1]*wness_temp);
          f_dw23_proj->SetParameter(2,sigma2_linear_pars[source][arm][charge][0]+sigma2_linear_pars[source][arm][charge][1]*wness_temp);
          f_dw23_proj->SetParameter(3,factor_linear_pars[source][arm][charge][0]+factor_linear_pars[source][arm][charge][1]*wness_temp);
          f_dw23_proj->SetParameter(4,poly_pars[source][arm][charge][0]+
              poly_pars[source][arm][charge][1]*(wness_temp)+
              poly_pars[source][arm][charge][2]*pow(wness_temp,2)+
              poly_pars[source][arm][charge][3]*pow(wness_temp,3)+
              poly_pars[source][arm][charge][4]*pow(wness_temp,4));
          sprintf(name,"%s %s Projection 0.9 < W_{ness} < 1.0",arm_label[arm],charge_label[charge]);
          h_dw23_sums[source][arm][charge][4]->SetTitle(name);
          h_dw23_sums[source][arm][charge][4]->SetTitleSize(.09);
          h_dw23_sums[source][arm][charge][4]->GetXaxis()->SetLabelSize(.06);
          h_dw23_sums[source][arm][charge][4]->GetXaxis()->SetTitle("dw_{23}");
          h_dw23_sums[source][arm][charge][4]->GetXaxis()->SetTitleSize(.06);
          h_dw23_sums[source][arm][charge][4]->GetYaxis()->SetLabelSize(.06);
          h_dw23_sums[source][arm][charge][4]->GetYaxis()->SetTitle("Yield (2D Normalized)");
          h_dw23_sums[source][arm][charge][4]->GetYaxis()->SetTitleSize(.06);
          h_dw23_sums[source][arm][charge][4]->Draw("ep");
          //f_dw23_proj->SetTitle("dw23 projection @ W_{ness}=0.95");
          f_dw23_proj->Draw("ACsame");
          f_dw23_proj->SetLineColor(2);

          c_dw23_fits_summary->cd((arm*2+charge)+17);
          h_dw23_sums[source][arm][charge][4]->Draw("ep");
          f_dw23_proj->Draw("ACsame");

          c_h2_dw23_vs_wness[source][arm][charge]->cd(6);
          TF1 *f_mean_tmp = new TF1("f_mean_temp","pol1(0)",0,1);
          f_mean_tmp->SetParameter(0,mean_linear_pars[source][arm][charge][0]);
          f_mean_tmp->SetParameter(1,mean_linear_pars[source][arm][charge][1]);
          f_mean_tmp->SetTitle("mean vs Wness");
          TGraph *g_mean_tmp = new TGraph(f_mean_tmp);
          g_mean_tmp->GetYaxis()->SetLabelSize(.06);
          g_mean_tmp->GetXaxis()->SetLabelSize(.06);
          g_mean_tmp->Draw("AC");

          c_h2_dw23_vs_wness[source][arm][charge]->cd(7);
          TF1 *f_sigma1_tmp = new TF1("f_sigma1_temp","pol1(0)",0,1);
          f_sigma1_tmp->SetParameter(0,sigma1_linear_pars[source][arm][charge][0]);
          f_sigma1_tmp->SetParameter(1,sigma1_linear_pars[source][arm][charge][1]);
          f_sigma1_tmp->SetTitle("Sigma 1 vs Wness");
          TGraph *g_sigma1_tmp = new TGraph(f_sigma1_tmp);
          g_sigma1_tmp->GetYaxis()->SetLabelSize(.06);
          g_sigma1_tmp->GetXaxis()->SetLabelSize(.06);
          g_sigma1_tmp->Draw("AC");

          c_h2_dw23_vs_wness[source][arm][charge]->cd(8);
          TF1 *f_sigma2_tmp = new TF1("f_sigma2_temp","pol1(0)",0,1);
          f_sigma2_tmp->SetParameter(0,sigma2_linear_pars[source][arm][charge][0]);
          f_sigma2_tmp->SetParameter(1,sigma2_linear_pars[source][arm][charge][1]);
          f_sigma2_tmp->SetTitle("Sigma 2 vs Wness");
          TGraph *g_sigma2_tmp = new TGraph(f_sigma2_tmp);
          g_sigma2_tmp->GetYaxis()->SetLabelSize(.06);
          g_sigma2_tmp->GetXaxis()->SetLabelSize(.06);
          g_sigma2_tmp->Draw("AC");

          c_h2_dw23_vs_wness[source][arm][charge]->cd(9);
          TF1 *f_factor_tmp = new TF1("f_factor_temp","pol1(0)",0,1);
          f_factor_tmp->SetParameter(0,factor_linear_pars[source][arm][charge][0]);
          f_factor_tmp->SetParameter(1,factor_linear_pars[source][arm][charge][1]);
          f_factor_tmp->SetTitle("Factor vs Wness");
          TGraph *g_factor_tmp = new TGraph(f_factor_tmp);
          g_factor_tmp->GetYaxis()->SetLabelSize(.06);
          g_factor_tmp->GetXaxis()->SetLabelSize(.06);
          g_factor_tmp->Draw("AC");

          sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/sbg_fit/temp/dw23_vs_wness_fit_arm%d_charge%d.png",arm,charge);
          c_h2_dw23_vs_wness[source][arm][charge]->SaveAs(name);
        }
      }
    }

    sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/sbg_fit/temp/dw23_fits_summary.png");
    c_dw23_fits_summary->SaveAs(name);
  }
  
  /*
  //construct hadronic background pdf
  h_proj_dw23[3][2][2];
  for(int source=0; source<1; source++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        sprintf(name,"h_proj_dw23_source%d_arm%d_charge%d",source,arm,charge);
        h_proj_dw23[source][arm][charge] = new TH1F(name,name,100,distmin[14],distmax[14]);

      }
    }
  }*/
  
  
  
  
/*  
  //do 1d fits on tgrapherrors
  int n_pars_1d = 4;
  double pars_1d[n_pars_1d];

  TCanvas *c_g_dw23[3][2][2];
  TGraphErrors *g_dw23[3][2][2][10];

  double dw23_entries[3][2][2][10][nhistbins];
  double dw23_axis[3][2][2][10][nhistbins];
  double dw23_err[3][2][2][10][nhistbins];
  double dw23_axis_err[3][2][2][10][nhistbins];
  
  double mean[3][2][2][10];
  double sigma1[3][2][2][10];
  double sigma2[3][2][2][10];
  double factor[3][2][2][10];

  double mean_err[3][2][2][10];
  double sigma1_err[3][2][2][10];
  double sigma2_err[3][2][2][10];
  double factor_err[3][2][2][10];

  double wness_section_arr[10] = {0.05,.15,.25,.35,.45,.55,.65,.75,.85,.95};
  double wness_err[10] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05};
  
  double linear_fit_pars[3][2][2][3][2];
  
  //DO 1D FITS ON HISTOGRAMS:
  TCanvas * c_h_dw23[3][2][2];

  for(int source=0; source < 1; source++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        sprintf(name,"c_h_dw23_source%d_arm%d_charge%d",source,arm,charge);
        c_h_dw23[source][arm][charge] = new TCanvas(name,name,1100,800);
        c_h_dw23[source][arm][charge]->Divide(4,3);
        for(int wness_section=0; wness_section<10; wness_section++) {
          TF1 *f_double_gaus_1d;
          f_double_gaus_1d = new TF1("f_double_gaus_1d",double_gaus_1d,distmin[14],distmax[14],n_pars_1d);

          if(wness_section>0) {
            c_h_dw23[source][arm][charge]->cd(wness_section);
          }

          double integral = h_dw23[source][arm][charge][wness_section]->Integral("width");
          if(integral != 0)
            h_dw23[source][arm][charge][wness_section]->Scale(1.0/integral);

          f_double_gaus_1d->SetParLimits(0,distmin[14],distmax[14]);
          f_double_gaus_1d->SetParLimits(1,.5,.05);
          f_double_gaus_1d->SetParLimits(2,.07,.01);
          f_double_gaus_1d->SetParLimits(3,-10,10);
          f_double_gaus_1d->SetParameter(0,h_dw23[source][arm][charge][wness_section]->GetMean());
          f_double_gaus_1d->SetParameter(1,h_dw23[source][arm][charge][wness_section]->GetRMS()*1.5);
          f_double_gaus_1d->SetParameter(2,h_dw23[source][arm][charge][wness_section]->GetRMS()*.9);
          f_double_gaus_1d->SetParameter(3,1);
          if(wness_section>0)
          h_dw23[source][arm][charge][wness_section]->Fit("f_double_gaus_1d","MN");
          
          f_double_gaus_1d->GetParameters(pars_1d);

          f_double_gaus_1d = new TF1("f_double_gaus_1d",double_gaus_1d,distmin[14],distmax[14],n_pars_1d);

          f_double_gaus_1d->SetParameter(0,pars_1d[0]);
          f_double_gaus_1d->SetParameter(1,pars_1d[1]);
          f_double_gaus_1d->SetParameter(2,pars_1d[2]);
          f_double_gaus_1d->SetParameter(3,pars_1d[3]);
            h_dw23[source][arm][charge][wness_section]->Fit("f_double_gaus_1d","MN");
          
          f_double_gaus_1d->GetParameters(pars_1d);
          double *par_err;
          par_err = (double*)f_double_gaus_1d->GetParErrors();
          
          mean[source][arm][charge][wness_section] = pars_1d[0];
          sigma1[source][arm][charge][wness_section] = pars_1d[1];
          sigma2[source][arm][charge][wness_section] = pars_1d[2];
          factor[source][arm][charge][wness_section] = pars_1d[3];
          
          mean_err[source][arm][charge][wness_section] = par_err[0];
          sigma1_err[source][arm][charge][wness_section] = par_err[1];
          sigma2_err[source][arm][charge][wness_section] = par_err[2];
          factor_err[source][arm][charge][wness_section] = par_err[3];
          
          if(wness_section>0) {
            sprintf(name,"dw23 fit %3.2f < Wness < %3.2f",wness_section/10.0,(wness_section+1)/10.0);
            h_dw23[source][arm][charge][wness_section]->SetTitle(name);
            h_dw23[source][arm][charge][wness_section]->Draw();
            
            TF1 *f_gaus1 = new TF1("f_gaus1","1/(sqrt(2*3.14159265)*([1]+[2]*[3]))*exp(-0.5*pow((x-[0])/[1],2))",distmin[14],distmax[14]);
            TF1 *f_gaus2 = new TF1("f_gaus2","1/(sqrt(2*3.14159265)*([1]+[2]*[3]))*[3]*exp(-0.5*pow((x-[0])/[2],2))",distmin[14],distmax[14]);
            
            for(int par=0; par<n_pars_1d; par++) {
              f_double_gaus_1d->SetParameter(par,pars_1d[par]);
              f_gaus1->SetParameter(par,pars_1d[par]);
              f_gaus2->SetParameter(par,pars_1d[par]);
            }
            
            f_double_gaus_1d->SetLineColor(2);
            f_gaus1->SetLineColor(3);
            f_gaus2->SetLineColor(4);
            f_gaus1->SetLineStyle(2);
            f_gaus2->SetLineStyle(2);
            
            f_double_gaus_1d->Draw("same");
            f_gaus1->Draw("same");
            f_gaus2->Draw("same");
            
          }
          if(wness_section == 9) {
            c_h_dw23[source][arm][charge]->cd(10);
            TGraphErrors *g_sigma1 = new TGraphErrors(10,wness_section_arr,sigma1[source][arm][charge],wness_err,sigma1_err[source][arm][charge]);
            g_sigma1->SetTitle("Sigma 1 vs Wness");
            g_sigma1->Draw("A*");
            TF1 * f_sigma1; f_sigma1 = new TF1("f_sigma1","pol1(0)",0.1,0.9);
            f_sigma1->SetParameter(0,sigma1[source][arm][charge][2]);
            f_sigma1->SetParameter(1,0);
            g_sigma1->Fit("f_sigma1","MRN");
            f_sigma1->GetParameters(linear_fit_pars[source][arm][charge][0]);
            f_sigma1->SetParameter(0,linear_fit_pars[source][arm][charge][0][0]);
            f_sigma1->SetParameter(1,linear_fit_pars[source][arm][charge][0][1]);
            g_sigma1->Fit("f_sigma1","MRN");
            f_sigma1->GetParameters(linear_fit_pars[source][arm][charge][0]);
            TF1 * f_sigma1_plot = new TF1("f_sigma1_plot","pol1(0)",0.1,0.9);
            f_sigma1_plot->SetParameter(0,linear_fit_pars[source][arm][charge][0][0]);
            f_sigma1_plot->SetParameter(1,linear_fit_pars[source][arm][charge][0][1]);
            f_sigma1_plot->Draw("same");
            
            c_h_dw23[source][arm][charge]->cd(11);
            TGraphErrors *g_sigma2 = new TGraphErrors(10,wness_section_arr,sigma2[source][arm][charge],wness_err,sigma2_err[source][arm][charge]);
            g_sigma2->SetTitle("Sigma 2 vs Wness");
            g_sigma2->Draw("A*");
            TF1 * f_sigma2; f_sigma2 = new TF1("f_sigma2","pol1(0)",0.1,0.9);
            f_sigma2->SetParameter(0,sigma2[source][arm][charge][2]);
            f_sigma2->SetParameter(1,0);
            g_sigma2->Fit("f_sigma2","MNR");
            f_sigma2->GetParameters(linear_fit_pars[source][arm][charge][1]);
            f_sigma2->SetParameter(0,linear_fit_pars[source][arm][charge][1][0]);
            f_sigma2->SetParameter(1,linear_fit_pars[source][arm][charge][1][1]);
            g_sigma2->Fit("f_sigma2","MNR");
            f_sigma2->GetParameters(linear_fit_pars[source][arm][charge][1]);
            TF1 * f_sigma2_plot = new TF1("f_sigma2_plot","pol1(0)",0.1,0.9);
            f_sigma2_plot->SetParameter(0,linear_fit_pars[source][arm][charge][1][0]);
            f_sigma2_plot->SetParameter(1,linear_fit_pars[source][arm][charge][1][1]);
            f_sigma2_plot->Draw("same");
            
            c_h_dw23[source][arm][charge]->cd(12);
            TGraphErrors *g_factor = new TGraphErrors(10,wness_section_arr,factor[source][arm][charge],wness_err,factor_err[source][arm][charge]);
            g_factor->SetTitle("Constant vs Wness");
            g_factor->Draw("A*");
            TF1 * f_factor; f_factor = new TF1("f_factor","pol1(0)",0.1,0.9);
            f_factor->SetParameter(0,factor[source][arm][charge][2]);
            f_factor->SetParameter(1,0);
            g_factor->Fit("f_factor","MNR");
            f_factor->GetParameters(linear_fit_pars[source][arm][charge][2]);
            f_factor->SetParameter(0,linear_fit_pars[source][arm][charge][2][0]);
            f_factor->SetParameter(1,linear_fit_pars[source][arm][charge][2][1]);
            g_factor->Fit("f_factor","MNR");
            f_factor->GetParameters(linear_fit_pars[source][arm][charge][2]);
            TF1 * f_factor_plot = new TF1("f_factor_plot","pol1(0)",0.1,0.9);
            f_factor_plot->SetParameter(0,linear_fit_pars[source][arm][charge][2][0]);
            f_factor_plot->SetParameter(1,linear_fit_pars[source][arm][charge][2][1]);
            f_factor_plot->Draw("same");
          }
        }
      }
    }
  }//end 1d fits on histograms
//*/

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

      double wness_mean;
      wness_mean = h_target_wness[0][arm][charge]->GetMean();
      //TF1 *f_hadbg_dw23_proj = new TF1("f_hadbg_dw23_proj",double_gaus_1d,distmin[14],distmax[14],4);
      TF1 *f_hadbg_dw23_proj = new TF1("f_hadbg_dw23_proj",double_gaus_1d,-.1,.1,4);
      f_hadbg_dw23_proj->SetParameter(0,mean_linear_pars[0][arm][charge][0]+mean_linear_pars[0][arm][charge][1]*.96);
      f_hadbg_dw23_proj->SetParameter(1,sigma1_linear_pars[0][arm][charge][0]+sigma1_linear_pars[0][arm][charge][1]*.96);
      f_hadbg_dw23_proj->SetParameter(2,sigma2_linear_pars[0][arm][charge][0]+sigma2_linear_pars[0][arm][charge][1]*.96);
      f_hadbg_dw23_proj->SetParameter(3,factor_linear_pars[0][arm][charge][0]+factor_linear_pars[0][arm][charge][1]*.96);

      //compose hadronic bkg dw23
      sprintf(name,"h_sbg_fit_dw23_source0_arm%d_charge%d",arm,charge);
      //h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]= new TH1F(name,name,1000,distmin[14],distmax[14]);
      h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]= new TH1F(name,name,nhistbins_sbg_dw23,-.1,.1);
      h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->Sumw2();
      h_sbg_fit_dw23[0][arm][charge][sbg_rpcdca_choice][sbg_fvtx_choice]->FillRandom("f_hadbg_dw23_proj",50000);


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



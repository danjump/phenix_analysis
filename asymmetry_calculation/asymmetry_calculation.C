#include <TROOT.h>
#include <TSystem.h>
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
#include "asymmetry_calculation.h"
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"
//#include "phys_hist_definitions.h"
#include "define_asymmetry_hists.h"
#include "spin_db_access/W2eGetPol.h"
#include "Math/GSLMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

double chisqr_fcn(const double *pars); //intended as a 4d function
double count_vs_asymmetry(int spin_b,int spin_y,const double * pars);

double pass_count[2][2];

void asymmetry_calculation( const std::string bkg_had_infilename) {
 
  double eta_threshold[4] = {1.1,1.4,1.8,2.6};
  //load library to access the spin database
 // gSystem->Load("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/install/lib/libW2eGetPol.la");


  char name[300]; sprintf(name,"prevent unused variable error");

  fstream linecounter;
  linecounter.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/badruns_spindb_run13.txt", fstream::in);
  char readbuffer[100];
  int count = -1;
  while(!linecounter.eof()) {
    linecounter.getline(readbuffer,100);
    count++;
  }
  linecounter.close();

  int num_badruns = count;
  int * bad_run = new int[num_badruns];

  fstream badrunlist_fstream;
  badrunlist_fstream.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/badruns_spindb_run13.txt",fstream::in);

  count=-1;
  while(!badrunlist_fstream.eof()) {
    badrunlist_fstream.getline(readbuffer,100);
    count++;
    if(!badrunlist_fstream.eof())
      bad_run[count] = atoi(readbuffer);
  }
  badrunlist_fstream.close();
  
  Int_t
  Run_Number,  Evt_Number,  triggerbit,
  clockcross,  fvtx_cone;
  
  Float_t 
  Evt_bbcZ,    Wness,
  charge,      pT,          pz,
  phi,         eta,         DG0,
  DDG0,        DG4,         chi2,
  DCA_z,       DCA_r,       dphi23,
  dw23,        Rpc1dca,     Rpc3dca,
  fvtx_dphi,   fvtx_dr,     fvtx_dtheta, fvtx_dr_dtheta;
  
  TFile *infile[11];
  std::string temp_str[11];

  temp_str[0] = bkg_had_infilename;
  
  infile[0]= new TFile(bkg_had_infilename.c_str());

  define_dw23_vs_eta_hists();

  //loop over the different input files and all events within each file to get distributions
  printf("Reading file:\n%s\n",temp_str[0].c_str());
  define_wness_tree(infile[0]);

  wness_tree->SetBranchAddress("Run_Number",&Run_Number);
  wness_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  wness_tree->SetBranchAddress("triggerbit",&triggerbit);
  wness_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  wness_tree->SetBranchAddress("clockcross",&clockcross);
  wness_tree->SetBranchAddress("Wness",&Wness);
  wness_tree->SetBranchAddress("charge",&charge);
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

  //variables to count events with each of the 4 helicity configurations 
  int asymmetry_count[2][2][4][4];
  for(int arm=0; arm<2; arm++) {
    for(int charge_index=0; charge_index<2; charge_index++) {
      for(int eta_index=0; eta_index<4; eta_index++) { // 0 is total eta range, 1-3 are segmented bins 
        for(int pat=0; pat<4; pat++) {
          asymmetry_count[arm][charge_index][eta_index][pat]=0;
        }
      }
    }
  }

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

    bool run_is_bad=false;
    for(int k=0; k<num_badruns; k++) {
      if(Run_Number==bad_run[k]) {
        run_is_bad=true;
        break;
      }
    }
    if(run_is_bad) {
      continue;
    }

    // Set arm (0 = South, 1 = North)
    int arm_index = 0;
    if		(pz < 0) arm_index = 0;
    else if (pz > 0) arm_index = 1;


    int charge_index = 0;
    if		(charge < 0) charge_index = 0;
    else if (charge > 0) charge_index = 1;
    
    bool passed_dw23_cut=false;
    if(charge <0) {
      if(-0.05 < dw23 && dw23 < 0.01)
        passed_dw23_cut=true;
    } else {
      if(-0.01 < dw23 && dw23 < 0.05)
        passed_dw23_cut=true;
    }

    int pat=0;
    if(Wness>0.92 && passed_dw23_cut ) {
      W2eGetPol spin_object;
      if(Run_Number != 391868) {
        int test = spin_object.InitPat(Run_Number);
        pat=0;
        if(test !=5) {
          pat = spin_object.GetPattern(clockcross);
          if(pat<4) {
            asymmetry_count[arm_index][charge_index][0][pat]++;

            //eta index 0 is total eta range, 1-3 are segmented bins
            if(eta >= eta_threshold[0] && eta <= eta_threshold[3]) {
              if(eta >= eta_threshold[1]) {
                if(eta >= eta_threshold[2]) {
                  asymmetry_count[arm_index][charge_index][3][pat]++;
                } else {
                  asymmetry_count[arm_index][charge_index][2][pat]++;
                }
              } else {
                asymmetry_count[arm_index][charge_index][1][pat]++;
              }
            }
            



            /*if(arm_index==1 && charge_index==1 && pat==0) {
              printf("%d\t%8d\t%d\t%d\t%f\t%f\t%f\n",
                  Run_Number,Evt_Number,clockcross,pat,dw23,pT,Wness);
            }*/
          }
        }
      }
    }
    
  } //end main event loop
  
  printf("\n\n");
  for(int eta_index=0; eta_index<4; eta_index++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        for(int pat=0; pat<4; pat++) {
          printf("Eta-bin %d Arm %d charge: %d pattern: %d count: %d\n",eta_index,arm,charge,pat,asymmetry_count[arm][charge][eta_index][pat]);
        }
      }
    }
  }

  /* Polarization configuration key:
   * blue=+ yellow=+ -> pat=0
   * blue=- yellow=+ -> pat=1
   * blue=+ yellow=- -> pat=2
   * blue=- yellow=- -> pat=3
   */

  //calculate raw asymmetries
  double raw_single_asymmetry[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ -/+ eta sign 0/1 ]
  double raw_double_asymmetry[2][2][4];    // [ arm ] [ charge ] [eta bins]

  double raw_single_asymmetry_err[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ -/+ eta sign 0/1 ]
  double raw_double_asymmetry_err[2][2][4];    // [ arm ] [ charge ] [eta bins]

  double raw_single_asymmetry_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ -/+ eta sign 0/1 ]
  double raw_double_asymmetry_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]

  double raw_single_asymmetry_err_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ -/+ eta sign 0/1 ]
  double raw_double_asymmetry_err_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]

  double numer, denom;


  for(int c=0; c<2; c++) { //charge
    for(int e=0; e<4; e++) { //eta bin (0 is full eta range. 1-3 are segmented bins)
      //Asymmetry calculations
      //south
      denom = asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][2] + asymmetry_count[0][c][e][3];

      numer = asymmetry_count[0][c][e][0] - asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][2] - asymmetry_count[0][c][e][3];
      raw_single_asymmetry[0][c][e][0] = numer/denom;

      numer = asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][1] - asymmetry_count[0][c][e][2] - asymmetry_count[0][c][e][3];
      raw_single_asymmetry[0][c][e][1] = numer/denom;

      numer = asymmetry_count[0][c][e][0] - asymmetry_count[0][c][e][1] - asymmetry_count[0][c][e][2] + asymmetry_count[0][c][e][3];
      raw_double_asymmetry[0][c][e] = numer/denom;

      //north
      denom = asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][2] + asymmetry_count[1][c][e][3];

      numer = asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][1] - asymmetry_count[1][c][e][2] - asymmetry_count[1][c][e][3];
      raw_single_asymmetry[1][c][e][0] = numer/denom;

      numer = asymmetry_count[1][c][e][0] - asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][2] - asymmetry_count[1][c][e][3];
      raw_single_asymmetry[1][c][e][1] = numer/denom;

      numer = asymmetry_count[1][c][e][0] - asymmetry_count[1][c][e][1] - asymmetry_count[1][c][e][2] + asymmetry_count[1][c][e][3];
      raw_double_asymmetry[1][c][e] = numer/denom;


      //Error calculations
      double err_denom,err_numer;
      //south
      err_denom = pow(asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][2] + asymmetry_count[0][c][e][3],4);

      err_numer = 4*((asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][2]) * (pow(asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][3],2))) -
        4*((pow(asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][2],2)) * (asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][3]));
      printf("0 %d 0 %f %f\n",c,err_numer,err_denom);
      raw_single_asymmetry_err[0][c][e][0] = sqrt(abs(err_numer/err_denom));

      err_numer = 4*((asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][1]) * (pow(asymmetry_count[0][c][e][2] + asymmetry_count[0][c][e][3],2))) -
        4*((pow(asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][1],2)) * (asymmetry_count[0][c][e][2] + asymmetry_count[0][c][e][3]));
      printf("0 %d 1 %f %f\n",c,err_numer,err_denom);
      raw_single_asymmetry_err[0][c][e][1] = sqrt(abs(err_numer/err_denom));

      err_numer = 4*((asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][3]) * (pow(asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][2],2))) -
        4*((pow(asymmetry_count[0][c][e][0] + asymmetry_count[0][c][e][3],2)) * (asymmetry_count[0][c][e][1] + asymmetry_count[0][c][e][2]));
      printf("0 %d %f %f\n",c,err_numer,err_denom);
      raw_double_asymmetry_err[0][c][e] = sqrt(abs(err_numer/err_denom));

      //north
      err_denom = pow(asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][2] + asymmetry_count[1][c][e][3],4);

      err_numer = 4*((asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][1]) * (pow(asymmetry_count[1][c][e][2] + asymmetry_count[1][c][e][3],2))) -
        4*((pow(asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][1],2)) * (asymmetry_count[1][c][e][2] + asymmetry_count[1][c][e][3]));
      printf("1 %d 0 %f %f\n",c,err_numer,err_denom);
      raw_single_asymmetry_err[1][c][e][0] = sqrt(abs(err_numer/err_denom));

      err_numer = 4*((asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][2]) * (pow(asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][3],2))) -
        4*((pow(asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][2],2)) * (asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][3]));
      printf("1 %d 1 %f %f\n",c,err_numer,err_denom);
      raw_single_asymmetry_err[1][c][e][1] = sqrt(abs(err_numer/err_denom));

      err_numer = 4*((asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][3]) * (pow(asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][2],2))) -
        4*((pow(asymmetry_count[1][c][e][0] + asymmetry_count[1][c][e][3],2)) * (asymmetry_count[1][c][e][1] + asymmetry_count[1][c][e][2]));
      printf("1 %d %f %f\n",c,err_numer,err_denom);
      raw_double_asymmetry_err[1][c][e] = sqrt(abs(err_numer/err_denom));
    }
  }

  printf("\n\n");
  char *eta_labels[4] = {(char*)"full eta range",(char*)"eta range: 1.1-1.4",(char*)"eta range: 1.4-1.8",(char*)"eta range: 1.8-2.6"};
  char *arm_name[2] = {(char*)"South",(char*)"North"};
  for(int eta_index=0; eta_index<4; eta_index++) {
    printf("------------------------------------------------------------------\n");
    printf("     Raw asymmetries for %s from numeric technique\n",eta_labels[eta_index]);
    printf("                     charge -          |          charge +\n");
    printf("            ------------------------------------------------------\n");
    for(int arm=0; arm<2; arm++) {
      printf("%s:\n",arm_name[arm]);
      printf("      A_ll:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_double_asymmetry[arm][charge][eta_index],raw_double_asymmetry_err[arm][charge][eta_index]);
      }
      printf("\n  A_l blue:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_single_asymmetry[arm][charge][eta_index][arm],raw_single_asymmetry_err[arm][charge][eta_index][arm]);
      }
      printf("\nA_l yellow:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_single_asymmetry[arm][charge][eta_index][1-arm],raw_single_asymmetry_err[arm][charge][eta_index][1-arm]);
      }
      printf("\n");
    }
    printf("------------------------------------------------------------------\n\n");

  }
  printf("\nAsymmetries from fit technique:\n");
  
  // Calculate Asymmetries with "Fit" Technique
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  //min->SetMaxFunctionCalls(1000000);
  //min->SetMaxIterations(100000);
  min->SetTolerance(0.01);

  ROOT::Math::Functor f(&chisqr_fcn,4); 
  double step[4] = {0.1,0.1,0.01,0.01}; 
  double variable[4] = {0.1,0.1,0.1,0.01};

  min->SetFunction(f);

  for(int eta_index=0; eta_index<4; eta_index++) {
    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        // Set the free variables to be minimized!
        min->SetVariable(0,"rel_lum",variable[0], step[0]);
        min->SetVariable(1,"al_b",variable[1], step[1]);
        min->SetVariable(2,"al_y",variable[2], step[2]);
        min->SetVariable(3,"a_ll",variable[3], step[3]);

        for(int spin_b=0; spin_b<=1; spin_b++) {
          for(int spin_y=0; spin_y<=1; spin_y++) {
            pass_count[spin_b][spin_y] = asymmetry_count[arm][charge][eta_index][(int)(1-spin_b)+2*(1-spin_y)];
            // asymmetry_count[arm][charge][pat] Polarization configuration key:
             // blue=+ yellow=+ -> pat=0
             // blue=- yellow=+ -> pat=1
             // blue=+ yellow=- -> pat=2
             // blue=- yellow=- -> pat=3
          }
        }

        min->Minimize(); 

        const double *par = min->X();
        const double *err = min->Errors();
        double par_err[2][4];
        for(int i=0; i<4; i++) {
          min->GetMinosError(i,par_err[0][i],par_err[1][i],0);
        }
        raw_double_asymmetry_fit[arm][charge][eta_index] = par[3];
        raw_double_asymmetry_err_fit[arm][charge][eta_index] = err[3];
        raw_single_asymmetry_fit[arm][charge][eta_index][arm] = par[1];
        raw_single_asymmetry_err_fit[arm][charge][eta_index][arm] = err[1];
        raw_single_asymmetry_fit[arm][charge][eta_index][1-arm] = par[2];
        raw_single_asymmetry_err_fit[arm][charge][eta_index][1-arm] = err[2];
        //printf("eta-bin %d arm%d charge%d \nal_b=%f-%f+%f (%f)\nal_y=%f-%f+%f (%f)  \na_ll=%f-%f+%f (%f)  \nrl_const=%f-%f+%f (%f)\n",
          //eta_index,arm,charge,
          //par[1],par_err[0][1],par_err[1][1],err[1],
          //par[2],par_err[0][2],par_err[1][2],err[2],
          //par[3],par_err[0][3],par_err[1][3],err[3],
          //par[0],par_err[0][0],par_err[1][0],err[0]);
      }
    }
    printf("\n");
  }

  for(int eta_index=0; eta_index<4; eta_index++) {
    printf("------------------------------------------------------------------\n");
    printf("         Raw asymmetries for %s from fit technique\n",eta_labels[eta_index]);
    printf("                     charge -          |          charge +\n");
    printf("            ------------------------------------------------------");
    for(int arm=0; arm<2; arm++) {
      printf("\n%s:\n",arm_name[arm]);
      printf("      A_ll:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_double_asymmetry_fit[arm][charge][eta_index],raw_double_asymmetry_err_fit[arm][charge][eta_index]);
      }
      printf("\n  A_l blue:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_single_asymmetry_fit[arm][charge][eta_index][arm],raw_single_asymmetry_err_fit[arm][charge][eta_index][arm]);
      }
      printf("\nA_l yellow:");
      for(int charge=0; charge<2; charge++) {
        printf("   %f +/- %f  ",raw_single_asymmetry_fit[arm][charge][eta_index][1-arm],raw_single_asymmetry_err_fit[arm][charge][eta_index][1-arm]);
      }
      printf("\n");
    }
    printf("------------------------------------------------------------------\n\n");

  }


  
}

double chisqr_fcn(const double *pars ) {
  //calculate chisquare
  double chisq = 0;
  double delta;
  for(int spin_b=-1; spin_b<=1; spin_b+=2) {
    for(int spin_y=-1; spin_y<=1; spin_y+=2) {
      delta  = (pass_count[(1+spin_b)/2][(1+spin_y)/2]-count_vs_asymmetry(spin_b,spin_y,pars))/sqrt(pass_count[(1+spin_b)/2][(1+spin_y)/2]);
      chisq += delta*delta;
    }
  }
  return chisq;

}

double count_vs_asymmetry(int spin_b,int spin_y,const double * pars) {
  int n_sig = 1;
  int n_bkg = 1;
  double p_b = 1;
  double p_y = 1;

  double count_eqn = pars[0]*n_bkg * (1 + n_sig/n_bkg * (spin_b*p_b*pars[1] + spin_y*p_y*pars[2] + spin_b*spin_y*p_b*p_y*pars[3]));

  return count_eqn;
}

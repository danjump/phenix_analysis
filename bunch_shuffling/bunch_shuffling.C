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
#include <TRandom3.h>
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
#include "bunch_shuffling.h"
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"
//#include "phys_hist_definitions.h"
#include "define_asymmetry_hists.h"
#include "W2eGetPol.h"
#include "spin_event.h"
#include "Math/GSLMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

double chisqr_fcn(const double *pars); //intended as a 4d function
double count_vs_asymmetry(int spin_b,int spin_y,const double * pars);
void shuffle_patterns(int unshuffled[120], int shuffled[120]);
void do_shuffling(int * unshuffled, int * shuffled, int entries);

double pass_count[2][2];
TRandom3 * randomizer = new TRandom3(9842);

void bunch_shuffling( const std::string bkg_had_infilename) {
 
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

  int good_count=0;

  printf("\nStarting Event Reading Loop...\n");
  printf("Number of events:  %d\n",entries);
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

  spin_event * spin_all = new spin_event[entries];

  //Get events from file, fill necessary histograms
  for(int i=0; i<entries; i++) {//Events loop
    //loop progress command line output
    percent_done_previous=percent_done;
    percent_done=(int)floor((float)(i+1)/(float)entries*(float)100);
    if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
      printf("%3i%% done",percent_done);
      if(percent_done==100) {
        printf("^_^");
        time( &rawtime );
        printf(" %s",ctime(&rawtime));
      } else {
        printf("...");
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
    
    int run=Run_Number;
    int evt=Evt_Number;
    
    int pat=0;
    if(Wness>0.92 && passed_dw23_cut ) {
      W2eGetPol spin_object;
      if(Run_Number != 391868) {
        int test = spin_object.InitPat(Run_Number);
        pat=0;
        if(test !=5) {
          pat = spin_object.GetPattern(clockcross);
          if(pat<4) {
            good_count++;
            //asymmetry_count[arm_index][charge_index][0][pat]++;

            //eta index 0 is total eta range, 1-3 are segmented bins
            int eta_index = 0;
            if(eta >= eta_threshold[0] && eta <= eta_threshold[3]) {
              if(eta >= eta_threshold[1]) {
                if(eta >= eta_threshold[2]) {
                  eta_index = 3;
                } else {
                  eta_index = 2;
                }
              } else {
                eta_index = 1;
              }
            }

            //fill spin_event's for this event
            //asymmetry_count[arm_index][charge_index][eta_index][pat]++;
            spin_all[i].fill(run,evt,clockcross,arm_index,charge_index,eta_index,pat);
            //printf("%d ",pat);
            



            /*if(arm_index==1 && charge_index==1 && pat==0) {
              printf("%d\t%8d\t%d\t%d\t%f\t%f\t%f\n",
                  Run_Number,Evt_Number,clockcross,pat,dw23,pT,Wness);
            }*/
          }
        }
      }
    }
    
  } //end read events loop
  
 



  //reduce spin event array to only good accepted events, and fill an array with the bunch pattern for each event
  int bunch_spins[good_count][120];
  spin_event * spin_good = new spin_event[good_count];
  int counter=0;

  percent_done=0;
  printf("\nStarting Spin Event Object Initialization Loop...\n");
  printf("Total spin object count: %d\n",good_count);
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));

  for(int i=0; i<entries; i++) {
    //loop progress command line output
    percent_done_previous=percent_done;
    percent_done=(int)floor((float)(i+1)/(float)entries*(float)100);
    if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
      printf("%3d%% done",percent_done);
      if(percent_done==100) {
        printf("^_^");
        time( &rawtime );
        printf(" %s",ctime(&rawtime));
      } else {
        printf("...");
        time( &rawtime );
        printf(" %s",ctime(&rawtime));
      }
    }

    if(spin_all[i].is_filled()) {
      spin_good[counter].copy(&spin_all[i]);
      //printf("%d ",spin_all[i].get_spin_config());

      W2eGetPol spin_object;
      spin_object.InitPat(spin_good[counter].get_run_num());
      spin_object.GetBunchPattern(bunch_spins[counter]);

      counter++;
    }
  }

  //printf("\n");

  delete [] spin_all;



  
  //Open output file for writing asymmetries
  TFile * file = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/bunch_shuffling/output/shuffled_1000.root","RECREATE");

  //Calculate actual asymmetries once before randomizing
  // NOTE: the following loop is unnecessary. It only runs it's contents once,
  // I just use it out of bad practice because I want to keep the Minimizer associated variables
  // in local scope so they don't conflict and cause a redeclaration when I run the minimizer 
  // later in the randomization loop. This is bad practice but oh well.
  
  TTree *t_real_asym = new TTree("t_real_asym","t_real_asym");
  int tmp_arm,tmp_charge,tmp_eta,which_asym;
  float asym,asym_err,scaled_asym;
  t_real_asym->Branch("arm",&tmp_arm,"arm/I");
  t_real_asym->Branch("charge",&tmp_charge,"charge/I");
  t_real_asym->Branch("eta_bin",&tmp_eta,"eta_bin/I");
  t_real_asym->Branch("which_asym",&which_asym,"which_asym/I");
  t_real_asym->Branch("asym",&asym,"asym/F");
  t_real_asym->Branch("asym_err",&asym_err,"asym_err/F");

  double raw_single_asymmetry_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ y=0 b=1 ]
  double raw_double_asymmetry_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]

  double raw_single_asymmetry_err_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ y=0 b=1 ]
  double raw_double_asymmetry_err_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]

  for(int unnecessary_loop=0; unnecessary_loop<1; unnecessary_loop++) {
    for(int i=0; i<good_count; i++) {//Asymmetry counting loop
      int arm_index = spin_good[i].get_arm();
      int charge_index = spin_good[i].get_charge_index();
      int pat = spin_good[i].get_spin_config();
      int eta_index = spin_good[i].get_eta_index();
      asymmetry_count[arm_index][charge_index][0][pat]++;

      if(eta_index>0) {
        asymmetry_count[arm_index][charge_index][eta_index][pat]++;
      }
    } //end asymmetry counting loop

    // Polarization configuration key:
    // blue=+ yellow=+ -> pat=0
    // blue=- yellow=+ -> pat=1
    // blue=+ yellow=- -> pat=2
    // blue=- yellow=- -> pat=3


    // Calculate Asymmetries with "Fit" Technique
    // this uses ROOT::Math::Minimizer to minimize a system of equations that consists
    // of 4 equations (for each condition) of asymmetry counts in terms of a_lb, a_ly, a_ll and a rel lum factor.
    // these equations are bundled in the count_vs_asymmetry() function and a chisquare of the system of equations
    // is defined in the chisqr_fcn() function which is what is actually minimized
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
          const double *err;
          err = min->Errors();
          double par_err[2][4];
          for(int i=0; i<4; i++) {
            min->GetMinosError(i,par_err[0][i],par_err[1][i],0);
          }
          
          raw_double_asymmetry_fit[arm][charge][eta_index] = par[3];
          raw_double_asymmetry_err_fit[arm][charge][eta_index] = err[3];
          raw_single_asymmetry_fit[arm][charge][eta_index][1] = par[1];    //blue
          raw_single_asymmetry_err_fit[arm][charge][eta_index][1] = err[1];
          raw_single_asymmetry_fit[arm][charge][eta_index][0] = par[2];    //yellow
          raw_single_asymmetry_err_fit[arm][charge][eta_index][0] = err[2];
          
          tmp_arm = arm;
          tmp_charge = charge;
          tmp_eta = eta_index;

          which_asym = 0;//sing asym y
          asym = par[2];
          asym_err = err[2];
          t_real_asym->Fill();

          which_asym = 1;//sing asym b
          asym = par[1];
          asym_err = err[1];
          t_real_asym->Fill();

          which_asym = 2;//double asym
          asym = par[3];
          asym_err = err[3];
          t_real_asym->Fill();

          if(eta_index==0) {
            //printf("a%dc%d y%f b%f d%f 00-%.0f 01-%.0f 10-%.0f 11-%.0f", arm, charge, par[2], par[1], par[3], pass_count[0][0],pass_count[0][1],pass_count[1][0],pass_count[1][1]);
          }
        }
      }

      if(eta_index==0) {
        //printf("\n");
      }
    }
  }






  //POLARIZATION RANDOMIZATION
  TTree * t_rand_asym = new TTree("t_rand_asym","t_rand_asym");
  int rand_index, rand_method;
  t_rand_asym->Branch("rand_method",&rand_method,"rand_method/I");
  t_rand_asym->Branch("rand_index",&rand_index,"rand_index/I");
  t_rand_asym->Branch("arm",&tmp_arm,"arm/I");
  t_rand_asym->Branch("charge",&tmp_charge,"charge/I");
  t_rand_asym->Branch("eta_bin",&tmp_eta,"eta_bin/I");
  t_rand_asym->Branch("which_asym",&which_asym,"which_asym/I");
  t_rand_asym->Branch("asym",&asym,"asym/F");
  t_rand_asym->Branch("asym_err",&asym_err,"asym_err/F");
  t_rand_asym->Branch("scaled_asym",&scaled_asym,"scaled_asym/F");
  
  // prepare histograms to be filled with asymmetry values for each randomization
  TH1F * h_single_asym_b[2][2][4]; // [arm][charge][eta bin]
  TH1F * h_single_asym_y[2][2][4]; // [arm][charge][eta bin]
  TH1F * h_double_asym[2][2][4];   // [arm][charge][eta bin]

  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      for(int eta_bin=0; eta_bin<4; eta_bin++) {
        sprintf(name,"h_single_asym_b_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_single_asym_b[arm][charge][eta_bin] = new TH1F(name,name,300,-.015,.015);

        sprintf(name,"h_single_asym_y_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_single_asym_y[arm][charge][eta_bin] = new TH1F(name,name,300,-.015,.015);

        sprintf(name,"h_double_asym_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_double_asym[arm][charge][eta_bin] = new TH1F(name,name,200,-.001,.001);
      }
    }
  }

  printf("\n");


  //how many iterations of randomization
  int num_rands = 1000;


  percent_done=0;
  percent_incriment=5;
  printf("Starting Randomization Loop...\nNumber of randomizations: %d\n",num_rands);
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));
  int over_count=0;

  for(int rand_index=0; rand_index<num_rands; rand_index++) { //randomization loop
    //loop progress command line output
    percent_done_previous=percent_done;
    percent_done=(int)floor((float)(rand_index+1)/(float)num_rands*(float)100);
    if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
      printf("%3d%% done",percent_done);
      if(percent_done==100) {
        printf("^_^");
        time( &rawtime );
        printf(" %s",ctime(&rawtime));
      } else {
        printf("...");
        time( &rawtime );
        printf(" %s",ctime(&rawtime));
      }
    }

    //First Randomization technique:
    // for each fill, shuffle the mapping of bunch->polarizations

    int shuffled_bunch_spins[120];
    int previous_bunch_spins[120];

    for(int i=0; i<120; i++) {
      previous_bunch_spins[i] = 0;
    }

    for(int arm_index = 0; arm_index<2; arm_index++) {
      for(int charge_index=0; charge_index<2; charge_index++) {
        for(int eta_index=0; eta_index<4; eta_index++) {
          for(int pat=0; pat<4; pat++) {
            asymmetry_count[arm_index][charge_index][eta_index][pat]=0;
          }
        }
      }
    }

    int good_event_pat_count=0;
    int bad_event_pat_count=0;

    for(int i=0; i<good_count; i++) {
      
      bool same_pat=true;
      for(int bunch=0; bunch<120; bunch++) {
        if(previous_bunch_spins[bunch] != bunch_spins[i][bunch]) {
          same_pat = false;
        }
      }
      
      //for each new set of runs with the same spin pattern, choose a random re-assignment of bunches
      //and corresponding spin pattern 
      if(!same_pat) {
        shuffle_patterns(bunch_spins[i],shuffled_bunch_spins);
      }
      
      //int pat = (int)randomizer->Integer(4);
      int pat = shuffled_bunch_spins[spin_good[i].get_clockcross()];
      if(pat >= 0 && pat <= 3) {

        // now add to the count of events for the proper criteria:
        int arm_index = spin_good[i].get_arm();
        int charge_index = spin_good[i].get_charge_index();
        int eta_index = spin_good[i].get_eta_index();
        asymmetry_count[arm_index][charge_index][0][pat]++;

        good_event_pat_count++;

        if(eta_index>0) {
          asymmetry_count[arm_index][charge_index][eta_index][pat]++;
        } 
      } else {
//        printf("Bad pat %d %d %d\n",i,spin_good[i].get_clockcross(),pat);
        bad_event_pat_count++;

      }

      for(int bunch=0; bunch<120; bunch++) {
        previous_bunch_spins[bunch] = bunch_spins[i][bunch];
      }
    } //end polarization randomization & asymmetry counting loop

    //printf("good pats: %d bad pats: %d\n",good_event_pat_count, bad_event_pat_count);

    // Polarization configuration key:
    // blue=+ yellow=+ -> pat=0
    // blue=- yellow=+ -> pat=1
    // blue=+ yellow=- -> pat=2
    // blue=- yellow=- -> pat=3


    // Calculate Asymmetries with "Fit" Technique
    // this uses ROOT::Math::Minimizer to minimize a system of equations that consists
    // of 4 equations (for each condition) of asymmetry counts in terms of a_lb, a_ly, a_ll and a rel lum factor.
    // these equations are bundled in the count_vs_asymmetry() function and a chisquare of the system of equations
    // is defined in the chisqr_fcn() function which is what is actually minimized
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
          const double *err;
          err = min->Errors();
          double par_err[2][4];
          for(int i=0; i<4; i++) {
            min->GetMinosError(i,par_err[0][i],par_err[1][i],0);
          }

          h_double_asym[arm][charge][eta_index]->Fill(par[3]);
          h_single_asym_b[arm][charge][eta_index]->Fill(par[1]);
          h_single_asym_y[arm][charge][eta_index]->Fill(par[2]);

          tmp_arm = arm;
          tmp_charge = charge;
          tmp_eta = eta_index;

          which_asym = 0;//sing asym y
          asym = par[2];
          asym_err = err[2];
          scaled_asym = par[2]/raw_single_asymmetry_err_fit[arm][charge][eta_index][0];
          t_rand_asym->Fill();

          which_asym = 1;//sing asym b
          asym = par[1];
          asym_err = err[1];
          scaled_asym = par[1]/raw_single_asymmetry_err_fit[arm][charge][eta_index][1];
          t_rand_asym->Fill();

          which_asym = 2;//double asym
          asym = par[3];
          asym_err = err[3];
          scaled_asym = par[3]/raw_double_asymmetry_err_fit[arm][charge][eta_index];
          t_rand_asym->Fill();
          
          if(eta_index==0) {
            //printf("a%dc%d y%f b%f l%f ", arm, charge, par[2], par[1], par[3]);
          }
        }
      }
      if(eta_index==0) {
        //printf("\n");
      }
    }
  }// end randomization loop

  printf("average over count: %f",(float)over_count/(float)num_rands);

  //save asymmetry distribution histograms
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      for(int eta_bin=0; eta_bin<4; eta_bin++) {
        sprintf(name,"c_asym_dist_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        
        h_single_asym_b[arm][charge][eta_bin]->Write();
        h_single_asym_y[arm][charge][eta_bin]->Write();
        h_double_asym[arm][charge][eta_bin]->Write();
      }
    }
  }

  t_real_asym->Write();
  t_rand_asym->Write();

  file->Close();
  
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

void shuffle_patterns(int unshuffled[120],int shuffled[120]) {
  int *good_mapping;
  int *good_spins;
  int *good_spins_shuffled;
  bool is_good_array[120];
  int good_count=0;
  for(int i=0; i<120; i++) {
    if(unshuffled[i] >= 0 && unshuffled[i] <= 3) {
      good_count++;
      is_good_array[i]=true;
    } else {
      shuffled[i] = unshuffled[i];
      is_good_array[i]=false;
    }
  }

  good_mapping = new int[good_count];
  good_spins = new int[good_count];
  good_spins_shuffled = new int[good_count];
  int good_index=0;
  for(int i=0; i<120; i++) {
    if(is_good_array[i]) {
      good_mapping[good_index] = i;
      good_spins[good_index] = unshuffled[i];
      good_index++;
    }
  }

  do_shuffling(good_spins,good_spins_shuffled,good_count);

  for(int i=0; i<good_count; i++) {
    shuffled[good_mapping[i]] = good_spins_shuffled[i];
  }

  delete [] good_mapping;
  delete [] good_spins;
  delete [] good_spins_shuffled;

  return;
}

void do_shuffling(int * unshuffled, int * shuffled, int entries) {
  if(entries==3) {
    shuffled[0] = unshuffled[1];
    shuffled[1] = unshuffled[2];
    shuffled[2] = unshuffled[0];

  } else if(entries==2) {
    shuffled[0] = unshuffled[1];
    shuffled[1] = unshuffled[0];

  } else if(entries<=1) {
    shuffled[0] = unshuffled[0];
  } else {

    int rand1 = randomizer->Integer(entries);
    int rand2 = randomizer->Integer(entries);
    while(rand1 == rand2)
      rand2 = randomizer->Integer(entries);

    shuffled[rand1] = unshuffled[rand2];
    shuffled[rand2] = unshuffled[rand1];

    int * sub_unshuffled = new int[entries-2];
    int * sub_shuffled = new int[entries-2];
    int offset=0;
    for(int i=0; i<entries-2; i++) {
      if(i+offset==rand1 || i+offset==rand2)
        offset++;
      sub_unshuffled[i] = unshuffled[i+offset];
    }
    
    do_shuffling(sub_unshuffled,sub_shuffled,entries-2);
    
    offset=0;
    for(int i=0; i<entries-2; i++) {
      if(i+offset==rand1 || i+offset==rand2)
        offset++;
      shuffled[i+offset] = sub_shuffled[i];
    }

    delete [] sub_unshuffled;
    delete [] sub_shuffled;
  }
  
  return;
}


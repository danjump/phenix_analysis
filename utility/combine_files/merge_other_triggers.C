#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TBranch.h>
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
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooFormulaVar.h>
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"

void merge_other_triggers(const std::string normal_tree_filename, const std::string other_tree_filename, const std::string treeoutfilename) {

  TFile * normal_tree_file = new TFile(normal_tree_filename.c_str());
  TFile * other_tree_file = new TFile(other_tree_filename.c_str());
  std::cout << "Normal tree input file:\n" << normal_tree_filename << std::endl << std::endl;
  std::cout << "Other tree input file:\n" << other_tree_filename << std::endl << std::endl;
  
  TFile *treeoutfile = new TFile(treeoutfilename.c_str(),"RECREATE");
  
  Int_t
    Run_Number,  Evt_Number,  triggerbit,
    clockcross,  fvtx_cone,   fvtx_tracklcone;

  Float_t
    Evt_bbcZ,    Wness,
    charge,      pT,          px,           py,         pz,
    phi,         eta,         DG0,
    DDG0,        DG4,         chi2,
    DCA_z,       DCA_r,       dphi12,       dphi23,
    dw23,        Rpc1dca,     Rpc1time,     Rpc3dca,    Rpc3time,
    Rpc1timewindow,           Rpc3timewindow,
    fvtx_dphi,   fvtx_dr,     fvtx_dtheta,
    fvtx_dr_dtheta, Rpc1x, Rpc1y, Rpc3x, Rpc3y;


  normal_tree_file->cd();
  define_basic_tree(normal_tree_file);
  TTree *normal_tree = (TTree*)basic_cuts_tree->Clone("normal_tree");
  
  other_tree_file->cd();
  define_basic_tree(other_tree_file);
  TTree *other_tree = (TTree*)basic_cuts_tree->Clone("other_tree");
  
  normal_tree->SetBranchAddress("Run_Number",&Run_Number);
  normal_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  normal_tree->SetBranchAddress("triggerbit",&triggerbit);
  normal_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  normal_tree->SetBranchAddress("clockcross",&clockcross);
  normal_tree->SetBranchAddress("Wness",&Wness);
  normal_tree->SetBranchAddress("charge",&charge);
  normal_tree->SetBranchAddress("pT",&pT);
  normal_tree->SetBranchAddress("px",&px);
  normal_tree->SetBranchAddress("pz",&py);
  normal_tree->SetBranchAddress("pz",&pz);
  normal_tree->SetBranchAddress("phi",&phi);
  normal_tree->SetBranchAddress("eta",&eta);
  normal_tree->SetBranchAddress("DG0",&DG0);
  normal_tree->SetBranchAddress("DDG0",&DDG0);
  normal_tree->SetBranchAddress("DG4",&DG4);
  normal_tree->SetBranchAddress("chi2",&chi2);
  normal_tree->SetBranchAddress("DCA_z",&DCA_z);
  normal_tree->SetBranchAddress("DCA_r",&DCA_r);
  normal_tree->SetBranchAddress("dphi12",&dphi12);
  normal_tree->SetBranchAddress("dphi23",&dphi23);
  normal_tree->SetBranchAddress("dw23",&dw23);
  normal_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  normal_tree->SetBranchAddress("Rpc1time",&Rpc1time);
  normal_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
  normal_tree->SetBranchAddress("Rpc1x",&Rpc1x);
  normal_tree->SetBranchAddress("Rpc1y",&Rpc1y);
  normal_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  normal_tree->SetBranchAddress("Rpc3time",&Rpc3time);
  normal_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
  normal_tree->SetBranchAddress("Rpc3x",&Rpc3x);
  normal_tree->SetBranchAddress("Rpc3y",&Rpc3y);
  normal_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  normal_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  normal_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  normal_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
  normal_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
  normal_tree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);
 
  other_tree->SetBranchAddress("Run_Number",&Run_Number);
  other_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  other_tree->SetBranchAddress("triggerbit",&triggerbit);
  other_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  other_tree->SetBranchAddress("clockcross",&clockcross);
  other_tree->SetBranchAddress("Wness",&Wness);
  other_tree->SetBranchAddress("charge",&charge);
  other_tree->SetBranchAddress("pT",&pT);
  other_tree->SetBranchAddress("px",&px);
  other_tree->SetBranchAddress("py",&py);
  other_tree->SetBranchAddress("pz",&pz);
  other_tree->SetBranchAddress("phi",&phi);
  other_tree->SetBranchAddress("eta",&eta);
  other_tree->SetBranchAddress("DG0",&DG0);
  other_tree->SetBranchAddress("DDG0",&DDG0);
  other_tree->SetBranchAddress("DG4",&DG4);
  other_tree->SetBranchAddress("chi2",&chi2);
  other_tree->SetBranchAddress("DCA_z",&DCA_z);
  other_tree->SetBranchAddress("DCA_r",&DCA_r);
  other_tree->SetBranchAddress("dphi12",&dphi12);
  other_tree->SetBranchAddress("dphi23",&dphi23);
  other_tree->SetBranchAddress("dw23",&dw23);
  other_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  other_tree->SetBranchAddress("Rpc1time",&Rpc1time);
  other_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
  other_tree->SetBranchAddress("Rpc1x",&Rpc1x);
  other_tree->SetBranchAddress("Rpc1y",&Rpc1y);
  other_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  other_tree->SetBranchAddress("Rpc3time",&Rpc3time);
  other_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
  other_tree->SetBranchAddress("Rpc3x",&Rpc3x);
  other_tree->SetBranchAddress("Rpc3y",&Rpc3y);
  other_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  other_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  other_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  other_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
  other_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
  other_tree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);

  treeoutfile->cd();
  define_basic_tree();
  basic_cuts_tree->SetBranchAddress("Run_Number",&Run_Number);
  basic_cuts_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  basic_cuts_tree->SetBranchAddress("triggerbit",&triggerbit);
  basic_cuts_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  basic_cuts_tree->SetBranchAddress("clockcross",&clockcross);
  basic_cuts_tree->SetBranchAddress("Wness",&Wness);
  basic_cuts_tree->SetBranchAddress("charge",&charge);
  basic_cuts_tree->SetBranchAddress("pT",&pT);
  basic_cuts_tree->SetBranchAddress("px",&px);
  basic_cuts_tree->SetBranchAddress("py",&py);
  basic_cuts_tree->SetBranchAddress("pz",&pz);
  basic_cuts_tree->SetBranchAddress("phi",&phi);
  basic_cuts_tree->SetBranchAddress("eta",&eta);
  basic_cuts_tree->SetBranchAddress("DG0",&DG0);
  basic_cuts_tree->SetBranchAddress("DDG0",&DDG0);
  basic_cuts_tree->SetBranchAddress("DG4",&DG4);
  basic_cuts_tree->SetBranchAddress("chi2",&chi2);
  basic_cuts_tree->SetBranchAddress("DCA_z",&DCA_z);
  basic_cuts_tree->SetBranchAddress("DCA_r",&DCA_r);
  basic_cuts_tree->SetBranchAddress("dphi12",&dphi12);
  basic_cuts_tree->SetBranchAddress("dphi23",&dphi23);
  basic_cuts_tree->SetBranchAddress("dw23",&dw23);
  basic_cuts_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  basic_cuts_tree->SetBranchAddress("Rpc1time",&Rpc1time);
  basic_cuts_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
  basic_cuts_tree->SetBranchAddress("Rpc1x",&Rpc1x);
  basic_cuts_tree->SetBranchAddress("Rpc1y",&Rpc1y);
  basic_cuts_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  basic_cuts_tree->SetBranchAddress("Rpc3time",&Rpc3time);
  basic_cuts_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
  basic_cuts_tree->SetBranchAddress("Rpc3x",&Rpc3x);
  basic_cuts_tree->SetBranchAddress("Rpc3y",&Rpc3y);
  basic_cuts_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  basic_cuts_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  basic_cuts_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  basic_cuts_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
  basic_cuts_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
  basic_cuts_tree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);

  int entries;
  int both_count=0;
  int other_count=0;
  int percent_done=0;
  int percent_incriment=20;
  int percent_done_previous;
  time_t rawtime;

  int unique_other=0;
  int unused_other=0;

  entries = other_tree->GetEntries();

  cout << "\nNumber of events:  " << entries << endl;
  cout << "Starting Main Event Loop...\n\n";
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));

  for(int i=0; i<entries; i++) {//Main loop
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

    other_tree->GetEntry(i);
    
    bool bit16=(triggerbit&0x00010000);
    bool bit17=(triggerbit&0x00020000); 
    bool bit18=(triggerbit&0x00040000); 
    bool bit19=(triggerbit&0x00080000); 
    bool bit20=(triggerbit&0x00100000); 
    bool bit21=(triggerbit&0x00200000); 
    bool bit22=(triggerbit&0x00400000); 
    bool bit23=(triggerbit&0x00800000); 
    bool bit24=(triggerbit&0x01000000); 
    bool bit25=(triggerbit&0x02000000); 
    bool bit26=(triggerbit&0x04000000); 
    bool bit27=(triggerbit&0x08000000); 
    bool bit28=(triggerbit&0x10000000); 

    bool ot_bits= (bit18 || bit19 || bit21);
    bool included_bits = (bit16 || bit17 || bit20 || bit22 || bit23 || bit24 || bit25 || bit26 || bit27 || bit28);

    bool proper_bits = ot_bits && !included_bits;

    if(proper_bits) {
      unique_other++;
      basic_cuts_tree->Fill();
    } else {
      unused_other++;
    }

  } //end main event loop

  entries = normal_tree->GetEntries();
  printf("\nunique_other: %d unused_other: %d normal_entries: %d combined_total: %d\n\n",unique_other,unused_other,entries,entries+unique_other);
  

  percent_done=0;
  percent_incriment=10;
  for(int i=0; i<entries; i++) {//Main loop
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
    normal_tree->GetEntry(i);

    basic_cuts_tree->Fill();
    
  } //end main event loop
    

  std::cout << "Saving new combined tree in:\n" << treeoutfilename << std::endl << std::endl;
  
  treeoutfile->cd();
  save_basic_tree();
  treeoutfile->Close();
}



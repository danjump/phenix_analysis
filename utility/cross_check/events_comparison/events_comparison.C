#include <TROOT.h>
#include <iomanip>
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
#include "events_comparison.h"
#include <vector>
#include "runevtobj.h"
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"

using namespace std;

void events_comparison() { //const std::string my_tree_filename, const std::string other_tree_filename) {
  TFile * my_tree_file = new TFile("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/current_wness_tree_data.root");//my_tree_filename.c_str());
  TFile * other_tree_file = new TFile("/phenix/u/beaumim/WanaTrees/Run13WnessTreeWithFvtxAndRpc_Data.root");//other_tree_filename.c_str());
  //std::cout << "my tree input file:\n" << my_tree_filename << std::endl << std::endl;
  //std::cout << "Other tree input file:\n" << other_tree_filename << std::endl << std::endl;
  
  char name[500];
  std::fstream in_both;
  in_both.open("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/output/in_both.txt",
      std::fstream::out | std::fstream::trunc );

  std::fstream in_mine_not_other;
  sprintf(name,"/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/output/in_daniel_not_mike.txt");
  in_mine_not_other.open(name, std::fstream::out | std::fstream::trunc );

  std::fstream in_other_not_mine;
  sprintf(name,"/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/output/in_mike_not_daniel.txt");
  in_other_not_mine.open(name, std::fstream::out | std::fstream::trunc );
  
  
  Int_t
    Run_Number,  Evt_Number,  triggerbit,
    clockcross,  fvtx_cone;

  Float_t
    Evt_bbcZ,    Wness,
    charge,      pT,          pz,
    phi,         eta,         DG0,
    DDG0,        DG4,         chi2,
    DCA_z,       DCA_r,       dphi12,       dphi23,
    dw23,        Rpc1dca,     Rpc3dca,   
    fvtx_dphi,   fvtx_dr,     fvtx_dtheta;

  my_tree_file->cd();
  TTree *temp = (TTree*) my_tree_file->Get("wness_tree");
  TTree *my_tree = (TTree*)temp->Clone("my_tree");
  
  other_tree_file->cd();
  temp = (TTree*) other_tree_file->Get("newsngmuons_basic_cut");
  TTree *other_tree = (TTree*)temp->Clone("other_tree");
 
  TTree *myout_tree = new TTree("myout_tree","myout_tree");
  TTree *otherout_tree = new TTree("otherout_tree","otherout_tree");

  my_tree->SetBranchAddress("Run_Number",&Run_Number);
  my_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  my_tree->SetBranchAddress("triggerbit",&triggerbit);
  my_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  my_tree->SetBranchAddress("clockcross",&clockcross);
  my_tree->SetBranchAddress("Wness",&Wness);
  my_tree->SetBranchAddress("charge",&charge);
  my_tree->SetBranchAddress("pT",&pT);
  my_tree->SetBranchAddress("pz",&pz);
  my_tree->SetBranchAddress("phi",&phi);
  my_tree->SetBranchAddress("eta",&eta);
  my_tree->SetBranchAddress("DG0",&DG0);
  my_tree->SetBranchAddress("DDG0",&DDG0);
  my_tree->SetBranchAddress("DG4",&DG4);
  my_tree->SetBranchAddress("chi2",&chi2);
  my_tree->SetBranchAddress("DCA_z",&DCA_z);
  my_tree->SetBranchAddress("DCA_r",&DCA_r);
  my_tree->SetBranchAddress("dphi12",&dphi12);
  my_tree->SetBranchAddress("dphi23",&dphi23);
  my_tree->SetBranchAddress("dw23",&dw23);
  my_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  my_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  my_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  my_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  my_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  my_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
 
  other_tree->SetBranchAddress("Run_Number",&Run_Number);
  other_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  other_tree->SetBranchAddress("triggerbit",&triggerbit);
  other_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  other_tree->SetBranchAddress("clockcross",&clockcross);
  other_tree->SetBranchAddress("Wness",&Wness);
  other_tree->SetBranchAddress("charge",&charge);
  other_tree->SetBranchAddress("pT",&pT);
  other_tree->SetBranchAddress("pz",&pz);
  other_tree->SetBranchAddress("phi",&phi);
  other_tree->SetBranchAddress("eta",&eta);
  other_tree->SetBranchAddress("DG0",&DG0);
  other_tree->SetBranchAddress("DDG0",&DDG0);
  other_tree->SetBranchAddress("DG4",&DG4);
  other_tree->SetBranchAddress("chi2",&chi2);
  other_tree->SetBranchAddress("DCA_z",&DCA_z);
  other_tree->SetBranchAddress("DCA_r",&DCA_r);
  other_tree->SetBranchAddress("dphi23",&dphi23);
  other_tree->SetBranchAddress("dw23",&dw23);
  other_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  other_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  other_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  other_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  other_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);

  myout_tree->Branch("Run_Number",&Run_Number,"Run_Number/I");
  myout_tree->Branch("Evt_Number",&Evt_Number,"Evt_Number/I");
  myout_tree->Branch("triggerbit",&triggerbit,"triggerbit/I");
  myout_tree->Branch("Evt_bbcZ",&Evt_bbcZ,"Evt_bbcZ/F");
  myout_tree->Branch("clockcross",&clockcross,"clockcross/I");
  myout_tree->Branch("Wness",&Wness,"Wness/F");
  myout_tree->Branch("charge",&charge,"charge/F");
  myout_tree->Branch("pT",&pT,"pT/F");
  myout_tree->Branch("pz",&pz,"pz/F");
  myout_tree->Branch("phi",&phi,"phi/F");
  myout_tree->Branch("eta",&eta,"eta/F");
  myout_tree->Branch("DG0",&DG0,"DG0/F");
  myout_tree->Branch("DDG0",&DDG0,"DDG0/F");
  myout_tree->Branch("DG4",&DG4,"DG4/F");
  myout_tree->Branch("chi2",&chi2,"chi2/F");
  myout_tree->Branch("DCA_z",&DCA_z,"DCA_z/F");
  myout_tree->Branch("DCA_r",&DCA_r,"DCA_r/F");
  myout_tree->Branch("dphi23",&dphi23,"dphi23/F");
  myout_tree->Branch("dw23",&dw23,"dw23/F");
  myout_tree->Branch("Rpc1dca",&Rpc1dca,"Rpc1dca/F");
  myout_tree->Branch("Rpc3dca",&Rpc3dca,"Rpc3dca/F");
  myout_tree->Branch("fvtx_dphi",&fvtx_dphi,"fvtx_dphi/F");
  myout_tree->Branch("fvtx_dr",&fvtx_dr,"fvtx_dr/F");
  myout_tree->Branch("fvtx_dtheta",&fvtx_dtheta,"fvtx_dtheta/F");

  otherout_tree->Branch("Run_Number",&Run_Number,"Run_Number/I");
  otherout_tree->Branch("Evt_Number",&Evt_Number,"Evt_Number/I");
  otherout_tree->Branch("triggerbit",&triggerbit,"triggerbit/I");
  otherout_tree->Branch("Evt_bbcZ",&Evt_bbcZ,"Evt_bbcZ/F");
  otherout_tree->Branch("clockcross",&clockcross,"clockcross/I");
  otherout_tree->Branch("Wness",&Wness,"Wness/F");
  otherout_tree->Branch("charge",&charge,"charge/F");
  otherout_tree->Branch("pT",&pT,"pT/F");
  otherout_tree->Branch("pz",&pz,"pz/F");
  otherout_tree->Branch("phi",&phi,"phi/F");
  otherout_tree->Branch("eta",&eta,"eta/F");
  otherout_tree->Branch("DG0",&DG0,"DG0/F");
  otherout_tree->Branch("DDG0",&DDG0,"DDG0/F");
  otherout_tree->Branch("DG4",&DG4,"DG4/F");
  otherout_tree->Branch("chi2",&chi2,"chi2/F");
  otherout_tree->Branch("DCA_z",&DCA_z,"DCA_z/F");
  otherout_tree->Branch("DCA_r",&DCA_r,"DCA_r/F");
  otherout_tree->Branch("dphi23",&dphi23,"dphi23/F");
  otherout_tree->Branch("dw23",&dw23,"dw23/F");
  otherout_tree->Branch("Rpc1dca",&Rpc1dca,"Rpc1dca/F");
  otherout_tree->Branch("Rpc3dca",&Rpc3dca,"Rpc3dca/F");
  otherout_tree->Branch("fvtx_dphi",&fvtx_dphi,"fvtx_dphi/F");
  otherout_tree->Branch("fvtx_dr",&fvtx_dr,"fvtx_dr/F");
  otherout_tree->Branch("fvtx_dtheta",&fvtx_dtheta,"fvtx_dtheta/F");

  other_tree->BuildIndex("Run_Number","Evt_Number");
  TTreeIndex *other_index_obj = (TTreeIndex*) other_tree->GetTreeIndex();

  printf("test2c\n");
  printf("test2c\n");

  my_tree->BuildIndex("Run_Number","Evt_Number");
  TTreeIndex *my_index_obj = (TTreeIndex*) my_tree->GetTreeIndex();

  printf("test3\n");
  
  time_t rawtime;

  cout << "Starting Main Event Loop...\n\n";
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));

  int my_index = 0;
  bool my_finished = false;

  int other_index = 0;
  bool other_finished = false;

  int my_count = 0;
  int other_count = 0;
  int both_count = 0;

  int percent_done=0;
  int percent_incriment=5;
  int percent_done_previous;
  int entries = my_index_obj->GetN();
  int loop_count = 0;

  while(!my_finished || !other_finished) {//Main loop
    //loop progress command line output
    percent_done_previous=percent_done;
    percent_done=(int)floor((float)(loop_count+1)/(float)entries*(float)100);
    if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
      printf("%3i%% done",percent_done);
      cout << "...";
      time( &rawtime );
      printf(" %s",ctime(&rawtime));
    }
    loop_count++;

    Long64_t other_entry = other_tree->LoadTree( other_index_obj->GetIndex()[other_index] ); 
    other_tree->GetEntry(other_entry);

    runevtobj otherobj(Run_Number,Evt_Number);

    Long64_t my_entry = my_tree->LoadTree( my_index_obj->GetIndex()[my_index] ); 
    my_tree->GetEntry(my_entry);

    runevtobj myobj(Run_Number,Evt_Number);
    
    bool in_both_trees;
    int lowest_run = 999999;
    int lowest_evt = 999999999;
    
    if(!my_finished && !other_finished) {
      if(myobj.run_num <= otherobj.run_num) {
        lowest_run = myobj.run_num;
        lowest_evt = myobj.evt_num;
        myobj.is_lowest = true;
        if(myobj.run_num == otherobj.run_num) {
          otherobj.is_lowest = true;
          in_both_trees = true;
        } else {
          in_both_trees = false;
        }
      } else {
        lowest_run = otherobj.run_num;
        lowest_evt = otherobj.evt_num;
        otherobj.is_lowest = true;
        in_both_trees = false;
      }
    } else {
      myobj.is_lowest = !my_finished;
      otherobj.is_lowest = !other_finished;
      in_both_trees = false;
    }

    if(!in_both_trees) {
      if(myobj.is_lowest) {
        my_count++;
        myout_tree->Fill();
        in_mine_not_other << std::setw(10) << Run_Number  << std::setw(10) << Evt_Number
          << std::setw(10) << triggerbit 
          << std::setw(10) << pT 
          << std::setw(10) << DCA_r 
          << std::setw(10) << DG0 
          << std::setw(10) << DDG0 
          << std::setw(10) << Rpc1dca 
          << std::setw(10) << Rpc3dca 
          << std::setw(10) << dw23 
          << std::setw(10) << eta 
          << std::setw(10) << Wness << "\n";
      } else {
        other_count++;
        other_tree->GetEntry(my_entry);
        otherout_tree->Fill();
        in_other_not_mine << std::setw(10) << Run_Number << std::setw(10) << Evt_Number 
          << std::setw(10) << triggerbit 
          << std::setw(10) << pT 
          << std::setw(10) << DCA_r 
          << std::setw(10) << DG0 
          << std::setw(10) << DDG0 
          << std::setw(10) << Rpc1dca 
          << std::setw(10) << Rpc3dca 
          << std::setw(10) << dw23 
          << std::setw(10) << eta 
          << std::setw(10) << Wness << "\n";
      }
    } else {
      both_count++;
      in_both << std::setw(10) << Run_Number << std::setw(10) << Evt_Number
          << std::setw(10) << triggerbit 
          << std::setw(10) << pT 
          << std::setw(10) << DCA_r 
          << std::setw(10) << DG0 
          << std::setw(10) << DDG0 
          << std::setw(10) << Rpc1dca 
          << std::setw(10) << Rpc3dca 
          << std::setw(10) << dw23 
          << std::setw(10) << eta 
          << std::setw(10) << Wness << "\n";
      other_tree->GetEntry(other_entry);
      in_both << std::setw(10) << Run_Number << std::setw(10) << Evt_Number
          << std::setw(10) << triggerbit 
          << std::setw(10) << pT 
          << std::setw(10) << DCA_r 
          << std::setw(10) << DG0 
          << std::setw(10) << DDG0 
          << std::setw(10) << Rpc1dca 
          << std::setw(10) << Rpc3dca 
          << std::setw(10) << dw23 
          << std::setw(10) << eta 
          << std::setw(10) << Wness << "\n";
    }

    if(my_index + myobj.is_lowest < my_index_obj->GetN()) {
      my_index += myobj.is_lowest;
    } else {
      my_finished = true;
    }

    if(other_index + otherobj.is_lowest < other_index_obj->GetN()) {
      other_index += otherobj.is_lowest;
    } else {
      other_finished = true;
    }

  } //end main event loop

  in_both.close();
  in_mine_not_other.close();
  in_other_not_mine.close();

  TFile * myout_file = new TFile(
      "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/output/tree_in_daniel_not_mike.root",
      "RECREATE");
  myout_tree->Write();
  myout_file->Close();
  TFile * otherout_file = new TFile(
      "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/output/tree_in_mike_not_daniel.root",
      "RECREATE");
  otherout_tree->Write();
  otherout_file->Close();


  printf("In both: %d In Daniel only: %d In Mike only: %d \n",both_count,my_count,other_count);

  
}



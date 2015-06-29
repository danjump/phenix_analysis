// this version of the code corresponds to output version06_01

#include "get_phys_dists.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TNtuple.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include <string>
#include <vector>
//NOTE: a number of global variables are defined in this .h:
#include "phys_hist_definitions.h"
#include "evaluate_triggers.h"


using namespace std;

Float_t get_rpc_time_window(int arm, int rpc, int run);

void get_phys_dists(const string & infilename,
const string &histoutfilename,
const string &treeoutfilename) { 
  
  //gROOT->Reset();
  gStyle->SetOptStat();
  
  //define cuts
  const int lastGap_cut  = 4;
  const int p_min_cut = 5;
  const int p_max_cut = 250;
  const int pt_min_cut = 16;
  const int pt_max_cut = 60;
  const int DG0_cut      = 20;
  const int DDG0_cut     = 9;
  const int chi2_cut = 20;
  const int rpc_dca_cut = 100;
  //cut counts
  int all_count = 0;
  int pass_lastGap_count = 0;
  int pass_pt_count = 0;
  int pass_p_count = 0;
  int pass_Evt_Nmu_count = 0;
  int pass_DG0_count = 0;
  int pass_DDG0_count = 0;
  int pass_chi2_count = 0;
  int pass_rpc_dca_count = 0;
  
  
  
  
  TFile *infile = new TFile(infilename.c_str());
  cout << "Reading file: " << endl << infilename.c_str()  << " ..." << endl;

  //Get main tree
  TTree *t_newsngmuons = (TTree*)infile->Get("newsngmuons");
  
  //Get branches and initialize related variables:	
  //   Get Event data branch
  TBranch *b_Eventdata = t_newsngmuons->GetBranch("Eventdata");
  b_Eventdata->SetAutoDelete(kTRUE);
  Float_t *Evt_bbcZ;
  Int_t *triggerbit,*triggerlive,*Run_Number,*Evt_Number,*clockcross;
  
  //   Get Muon Tracks branch
  TBranch *b__RecoTracks = t_newsngmuons->GetBranch("_RecoTracks");
  Int_t _RecoTracks;
  b__RecoTracks->SetAddress(&_RecoTracks);
  b__RecoTracks->SetAutoDelete(kTRUE);
  TBranch *b_RecoTracks = t_newsngmuons->GetBranch("RecoTracks");
  b_RecoTracks->SetAutoDelete(kTRUE);
  Float_t *charge, *px, *py, *pz, *pT, *p;
  Float_t *eta, *phi, *DG0, *DDG0;
  Float_t *DG4, *chi2, *DCA_z, *DCA_r;
  Float_t *st1y, *st1x, *st2y, *st2x, *st3y,*st3x;
  Int_t *lastGap, *Evt_Nmu;
  
  //   Get Fvtx Branch
  TBranch *b_FvtxMatch = t_newsngmuons->GetBranch("FvtxMatch");
  b_FvtxMatch->SetAutoDelete(kTRUE);
  Float_t *fvtx_dphi, *fvtx_dtheta, *fvtx_dr;
  Int_t *fvtx_tracklconebits, *fvtx_tracklcone, *fvtx_conebits,  *fvtx_cone;
  
  //   Get Various MuTr->RPC projection branches
  TBranch *b_RpcMatchSt1 = t_newsngmuons->GetBranch("RpcMatchSt1");
  b_RpcMatchSt1->SetAutoDelete(kTRUE);
  TBranch *b_RpcMatchSt3 = t_newsngmuons->GetBranch("RpcMatchSt3");
  b_RpcMatchSt3->SetAutoDelete(kTRUE);
  Float_t *Rpc3dca,*Rpc3x, *Rpc3y, *Rpc3time, *Rpc1dca, *Rpc1x, *Rpc1y, *Rpc1time;
  
  
  std::cout << "Creating histogram output file: \n" << histoutfilename << "\n";
  //TFile *histoutfile = new TFile(histoutfilename.c_str(), "RECREATE");

  // define physics hists
  //define_hists_phys();
  
  std::cout << "Creating tree output file: \n" << treeoutfilename << "\n\n";
  TFile *treeoutfile = new TFile(treeoutfilename.c_str(), "RECREATE");
  //TNtuple * basic_cuts_ntuple = new TNtuple();
  define_basic_tree();
  
  //MAIN LOOP
  int entries = t_newsngmuons->GetEntries();
  
  //   general variables
  int percent_done=0;
  int percent_incriment=20;
  int percent_done_previous;
  
  time_t rawtime;
  
  cout << "\nNumber of events:  " << entries << endl;
  cout << "Starting Main Event Loop...\n\n";
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));
  for(int i=0; i<entries; i++) { //Main loop
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
    //Set leaf addresses & get entries
    //   Eventdata
    Evt_bbcZ     = new Float_t[1];
    Evt_Number   = new Int_t[1];
    clockcross   = new Int_t[1];
    triggerbit   = new Int_t[1];
    triggerlive  = new Int_t[1];
    Run_Number   = new Int_t[1];
    ((TLeaf*)b_Eventdata->GetLeaf("Evt_bbcZ"))->SetAddress(Evt_bbcZ);
    ((TLeaf*)b_Eventdata->GetLeaf("Evt_Number"))->SetAddress(Evt_Number);
    ((TLeaf*)b_Eventdata->GetLeaf("clockcross"))->SetAddress(clockcross);
    ((TLeaf*)b_Eventdata->GetLeaf("triggerbit"))->SetAddress(triggerbit);
    ((TLeaf*)b_Eventdata->GetLeaf("triggerlive"))->SetAddress(triggerlive);
    ((TLeaf*)b_Eventdata->GetLeaf("Run_Number"))->SetAddress(Run_Number);
    b_Eventdata->GetEntry(i);
    
    //   Muon Tracks
    b__RecoTracks->GetEntry(i);
    Evt_Nmu	 = new Int_t[_RecoTracks];
    charge   = new Float_t[_RecoTracks];
    px		   = new Float_t[_RecoTracks];
    py		   = new Float_t[_RecoTracks];
    pz		   = new Float_t[_RecoTracks];
    pT		   = new Float_t[_RecoTracks];
    p		     = new Float_t[_RecoTracks];
    eta      = new Float_t[_RecoTracks];
    phi      = new Float_t[_RecoTracks];
    DG0      = new Float_t[_RecoTracks];
    DDG0     = new Float_t[_RecoTracks];
    lastGap	 = new Int_t[_RecoTracks];	
    DG4      = new Float_t[_RecoTracks];
    chi2     = new Float_t[_RecoTracks];
    DCA_r    = new Float_t[_RecoTracks];
    DCA_z    = new Float_t[_RecoTracks];
    st1y     = new Float_t[_RecoTracks];
    st1x     = new Float_t[_RecoTracks];
    st2y     = new Float_t[_RecoTracks];
    st2x     = new Float_t[_RecoTracks];
    st3y     = new Float_t[_RecoTracks];
    st3x     = new Float_t[_RecoTracks];
    
    ((TLeaf*)b_RecoTracks->GetLeaf("Evt_Nmu"))->SetAddress(Evt_Nmu);
    ((TLeaf*)b_RecoTracks->GetLeaf("charge"))->SetAddress(charge);
    ((TLeaf*)b_RecoTracks->GetLeaf("px"))->SetAddress(px);
    ((TLeaf*)b_RecoTracks->GetLeaf("py"))->SetAddress(py);
    ((TLeaf*)b_RecoTracks->GetLeaf("pz"))->SetAddress(pz);
    ((TLeaf*)b_RecoTracks->GetLeaf("pT"))->SetAddress(pT);
    ((TLeaf*)b_RecoTracks->GetLeaf("p"))->SetAddress(p);
    ((TLeaf*)b_RecoTracks->GetLeaf("eta"))->SetAddress(eta);
    ((TLeaf*)b_RecoTracks->GetLeaf("phi"))->SetAddress(phi);
    ((TLeaf*)b_RecoTracks->GetLeaf("DG0"))->SetAddress(DG0);
    ((TLeaf*)b_RecoTracks->GetLeaf("DDG0"))->SetAddress(DDG0);
    ((TLeaf*)b_RecoTracks->GetLeaf("lastGap"))->SetAddress(lastGap);
    ((TLeaf*)b_RecoTracks->GetLeaf("DG4"))->SetAddress(DG4);
    ((TLeaf*)b_RecoTracks->GetLeaf("chi2"))->SetAddress(chi2);
    ((TLeaf*)b_RecoTracks->GetLeaf("DCA_z"))->SetAddress(DCA_z);
    ((TLeaf*)b_RecoTracks->GetLeaf("DCA_r"))->SetAddress(DCA_r);
    ((TLeaf*)b_RecoTracks->GetLeaf("ySt1"))->SetAddress(st1y);
    ((TLeaf*)b_RecoTracks->GetLeaf("xSt1"))->SetAddress(st1x);
    ((TLeaf*)b_RecoTracks->GetLeaf("ySt2"))->SetAddress(st2y);
    ((TLeaf*)b_RecoTracks->GetLeaf("xSt2"))->SetAddress(st2x);
    ((TLeaf*)b_RecoTracks->GetLeaf("ySt3"))->SetAddress(st3y);
    ((TLeaf*)b_RecoTracks->GetLeaf("xSt3"))->SetAddress(st3x);

    
    b_RecoTracks->GetEntry(i);
    
    //   Different MuTr projections to rpc 3	
    //   [0]=MatchVtx	
    //   [1]=MatchSt1
    //   [2]=MatchSt3
    //   [3]=MatchMuIDdca[proj_index]  = new Float_t[_RecoTracks];
    Rpc1dca  = new Float_t[_RecoTracks];
    Rpc1x  = new Float_t[_RecoTracks];
    Rpc1y  = new Float_t[_RecoTracks];
    Rpc1time  = new Float_t[_RecoTracks];
    Rpc3dca = new Float_t[_RecoTracks];
    Rpc3x = new Float_t[_RecoTracks];
    Rpc3y = new Float_t[_RecoTracks];
    Rpc3time = new Float_t[_RecoTracks];
    ((TLeaf*)b_RpcMatchSt1->GetLeaf("Rpc1dca"))->SetAddress(Rpc1dca);
    ((TLeaf*)b_RpcMatchSt1->GetLeaf("Rpc1x"))->SetAddress(Rpc1x);
    ((TLeaf*)b_RpcMatchSt1->GetLeaf("Rpc1y"))->SetAddress(Rpc1y);
    ((TLeaf*)b_RpcMatchSt1->GetLeaf("Rpc1time"))->SetAddress(Rpc1time);
    ((TLeaf*)b_RpcMatchSt3->GetLeaf("Rpc3dca"))->SetAddress(Rpc3dca);
    ((TLeaf*)b_RpcMatchSt3->GetLeaf("Rpc3x"))->SetAddress(Rpc3x);
    ((TLeaf*)b_RpcMatchSt3->GetLeaf("Rpc3y"))->SetAddress(Rpc3y);
    ((TLeaf*)b_RpcMatchSt3->GetLeaf("Rpc3time"))->SetAddress(Rpc3time);
    b_RpcMatchSt1->GetEntry(i);
    b_RpcMatchSt3->GetEntry(i);
    
    fvtx_dphi = new Float_t[_RecoTracks];
    fvtx_dtheta = new Float_t[_RecoTracks];
    fvtx_dr = new Float_t[_RecoTracks];
    fvtx_conebits = new Int_t[_RecoTracks];
    fvtx_tracklconebits = new Int_t[_RecoTracks];
    ((TLeaf*)b_FvtxMatch->GetLeaf("fvtx_dphi"))->SetAddress(fvtx_dphi);
    ((TLeaf*)b_FvtxMatch->GetLeaf("fvtx_dtheta"))->SetAddress(fvtx_dtheta);
    ((TLeaf*)b_FvtxMatch->GetLeaf("fvtx_dr"))->SetAddress(fvtx_dr);
    ((TLeaf*)b_FvtxMatch->GetLeaf("fvtx_conebits"))->SetAddress(fvtx_conebits);
    ((TLeaf*)b_FvtxMatch->GetLeaf("fvtx_tracklconebits"))->SetAddress(fvtx_tracklconebits);
    b_FvtxMatch->GetEntry(i);
    
    fvtx_cone = new Int_t[_RecoTracks];
    fvtx_tracklcone = new Int_t[_RecoTracks];
    
    
    //loop over tracks in an event
    for(int mt_track=0; mt_track<_RecoTracks; mt_track++) { //Muon tracks loop
      
      //apply basic cuts
      all_count++;
      if (lastGap[mt_track] < lastGap_cut)   continue;  //require all muid gaps hit
      pass_lastGap_count++;
      if (p[mt_track] < p_min_cut || p[mt_track] > p_max_cut)  continue;
      pass_p_count++;
      if (pT[mt_track] < pt_min_cut || pT[mt_track] > pt_max_cut)  continue;
      pass_pt_count++;
      if (Evt_Nmu[mt_track] != 1)           continue;
      pass_Evt_Nmu_count++;
      if (DG0[mt_track] > DG0_cut)           continue;
      pass_DG0_count++;
      if (DDG0[mt_track] > DDG0_cut)         continue;
      pass_DDG0_count++;
      if (chi2[mt_track] > chi2_cut)               continue;
      pass_chi2_count++;
      
      int rpcdca_condition;
      if((Rpc1dca[mt_track] < rpc_dca_cut) && (Rpc3dca[mt_track] < rpc_dca_cut)) {
        rpcdca_condition = 2;
      } else if(Rpc1dca[mt_track] < rpc_dca_cut) {
        rpcdca_condition = 0;
      } else if(Rpc3dca[mt_track] < rpc_dca_cut) {
        rpcdca_condition = 1;
      } else {
        rpcdca_condition = -1;
        continue;
      }
      pass_rpc_dca_count++;
      
      int good_fvtx;
      //fvtx_dphi[mt_track] = fabs(fvtx_dphi[mt_track]);
      if(-1.5 < fvtx_dphi[mt_track] && fvtx_dphi[mt_track] < 1.5 && 
          0 < fvtx_dtheta[mt_track] && fvtx_dtheta[mt_track] < 1.5  &&
          0 < fvtx_dr[mt_track] && fvtx_dr[mt_track] < 100) {
        good_fvtx = 1;
      } else {
        good_fvtx = 0;
      }
      
      
      // Set arm (0 = South, 1 = North)
      int track_arm = 0;
      if		(pz[mt_track] < 0) track_arm = 0;
      else if (pz[mt_track] > 0) track_arm = 1;
      else cout << "warning: wrong track_arm info detected!" <<endl;
      
      
      int charge_index = 0;
      if(charge[mt_track] < 0) charge_index = 0;
      else if (charge[mt_track] > 0) charge_index = 1;
      
     fvtx_cone[mt_track] = 0;
      for(int test =2; test< 7;test++) {
        fvtx_cone[mt_track] +=(int)((fvtx_conebits[mt_track] >> test*4)&0xf);
        fvtx_tracklcone[mt_track] +=(int)((fvtx_tracklconebits[mt_track] >> test*4)&0xf);
      }

      Float_t Rpc1timewindow = get_rpc_time_window(track_arm,0,(int)Run_Number[mt_track]);
      Float_t Rpc3timewindow = get_rpc_time_window(track_arm,1,(int)Run_Number[mt_track]);
      
      Float_t wness = 0;
      Float_t fill_eta = fabs(eta[mt_track]);
      Float_t fill_dphi12 = (fmod(atan2(st2y[mt_track],st2x[mt_track])-atan2(st1y[mt_track],st1x[mt_track])+3*PI,2*PI)-PI);
      Float_t fill_dphi23 = (fmod(atan2(st3y[mt_track],st3x[mt_track])-atan2(st2y[mt_track],st2x[mt_track])+3*PI,2*PI)-PI);
      Float_t fill_dw23 = fill_dphi23*pT[mt_track]*sin(atan(pT[mt_track]/fabs(pz[mt_track])));
      Float_t fill_dr_dtheta = fvtx_dr[mt_track]*fvtx_dtheta[mt_track];

      // physics histograms
      basic_cuts_tree->SetBranchAddress("Run_Number", &Run_Number[mt_track]);
      basic_cuts_tree->SetBranchAddress("Evt_Number", &Evt_Number[mt_track]);
      basic_cuts_tree->SetBranchAddress("triggerbit", &triggerbit[mt_track]);
      basic_cuts_tree->SetBranchAddress("Evt_bbcZ", &Evt_bbcZ[mt_track]);
      basic_cuts_tree->SetBranchAddress("clockcross", &clockcross[mt_track]);
      basic_cuts_tree->SetBranchAddress("Wness", &wness);
      basic_cuts_tree->SetBranchAddress("charge", &charge[mt_track]);
      basic_cuts_tree->SetBranchAddress("pT", &pT[mt_track]);
      basic_cuts_tree->SetBranchAddress("px", &px[mt_track]);
      basic_cuts_tree->SetBranchAddress("py", &py[mt_track]);
      basic_cuts_tree->SetBranchAddress("pz", &pz[mt_track]);
      basic_cuts_tree->SetBranchAddress("phi", &phi[mt_track]);
      basic_cuts_tree->SetBranchAddress("eta", &fill_eta);
      basic_cuts_tree->SetBranchAddress("DG0", &DG0[mt_track]);
      basic_cuts_tree->SetBranchAddress("DDG0", &DDG0[mt_track]);
      basic_cuts_tree->SetBranchAddress("DG4", &DG4[mt_track]);
      basic_cuts_tree->SetBranchAddress("chi2", &chi2[mt_track]);
      basic_cuts_tree->SetBranchAddress("DCA_z", &DCA_z[mt_track]);
      basic_cuts_tree->SetBranchAddress("DCA_r", &DCA_r[mt_track]);
      basic_cuts_tree->SetBranchAddress("dphi12", &fill_dphi12);
      basic_cuts_tree->SetBranchAddress("dphi23", &fill_dphi23);
      basic_cuts_tree->SetBranchAddress("dw23", &fill_dw23);
      basic_cuts_tree->SetBranchAddress("Rpc1dca", &Rpc1dca[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc1x", &Rpc1x[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc1y", &Rpc1y[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc1time", &Rpc1time[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc1timewindow", &Rpc1timewindow);
      basic_cuts_tree->SetBranchAddress("Rpc3dca", &Rpc3dca[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc3x", &Rpc3x[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc3y", &Rpc3y[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc3time", &Rpc3time[mt_track]);
      basic_cuts_tree->SetBranchAddress("Rpc3timewindow", &Rpc3timewindow);
      basic_cuts_tree->SetBranchAddress("fvtx_dphi", &fvtx_dphi[mt_track]);
      basic_cuts_tree->SetBranchAddress("fvtx_dr", &fvtx_dr[mt_track]);
      basic_cuts_tree->SetBranchAddress("fvtx_dtheta", &fvtx_dtheta[mt_track]);
      basic_cuts_tree->SetBranchAddress("fvtx_dr_dtheta", &fill_dr_dtheta);
      basic_cuts_tree->SetBranchAddress("fvtx_cone", &fvtx_cone[mt_track]);
      basic_cuts_tree->SetBranchAddress("fvtx_tracklcone", &fvtx_tracklcone[mt_track]);
      
      basic_cuts_tree->Fill();
      


      //	 bool test[4] = {didrpc3bctrigger(triggerbit[0]),didrpc3atrigger(triggerbit[0]),didrpc1trigger(triggerbit[0]),didsg1trigger(triggerbit[0])};

      /*bool trigger_condition[2][5] = {
        {didrpc3bcStrigger(triggerbit[0]),
          didrpc3aStrigger(triggerbit[0]),
          didrpc1trigger(triggerbit[0]),
          didsg1Strigger(triggerbit[0]),
          didotherStrigger(triggerbit[0])},
        
        {didrpc3bcNtrigger(triggerbit[0]),
          didrpc3aNtrigger(triggerbit[0]),
          didrpc1trigger(triggerbit[0]),
          didsg1Ntrigger(triggerbit[0]),
          didotherNtrigger(triggerbit[0])} };*/


      //for (int trigs = 0; trigs < n_trig; trigs++)
      //{
      
      //if (trigger_condition[track_arm][trigs]) {
      /*h2_dg0_ddg0[track_arm][charge_index][rpcdca_condition][good_fvtx]->Fill(DDG0[mt_track],DG0[mt_track]);
      h2_chi2_dcar[track_arm][charge_index][rpcdca_condition][good_fvtx]->Fill(DCA_r[mt_track],chi2[mt_track]);
      for (int dists = 0; dists < n_dists; dists++) {
        h_phys_dists[track_arm][charge_index][rpcdca_condition][good_fvtx][dists]->Fill(phys_dist[dists]);
      }*/
      //}
      //}


    }//end muon tracks loop
    
    //clean up leaf arrays
    delete [] Evt_bbcZ;
    delete [] clockcross;
    delete [] triggerbit;
    delete [] triggerlive;
    delete [] Run_Number;
    delete [] Evt_Number;
    delete [] charge;
    delete [] px;
    delete [] py;
    delete [] pz;
    delete [] pT;
    delete [] p;  
    delete [] lastGap;    
    delete [] chi2;     
    delete [] DCA_z;     
    delete [] DCA_r;
    delete [] eta;
    delete [] phi;
    delete [] DG0;
    delete [] DDG0;
    delete [] DG4; 
    delete [] st1y;     
    delete [] st1x;     
    delete []  st2y;     
    delete []  st2x;     
    delete []  st3y;     
    delete []  st3x;   
    delete [] Rpc3dca;
    delete [] Rpc3x;
    delete [] Rpc3y;
    delete [] Rpc3time;
    delete [] Rpc1dca;
    delete [] Rpc1x;
    delete [] Rpc1y;
    delete [] Rpc1time;
    delete [] fvtx_cone;
    delete [] fvtx_conebits;
    delete [] fvtx_tracklcone;
    delete [] fvtx_tracklconebits;
    



  }//end main loop
  
  cout << endl << "CUT INFORMATION:" << endl;
  printf("Cut condition | Events Passed | %% of total | %% of previous\n");
  printf(" All tracks   | %13d | %10.2f | %10.2f\n", all_count, 100.0, 100.0);
  printf(" LastGap cut  | %13d | %10.2f | %10.2f\n", pass_lastGap_count, 100.0*(double)pass_lastGap_count/(double)all_count, 100.0*(double)pass_lastGap_count/(double)all_count);
  printf(" pT cut       | %13d | %10.2f | %10.2f\n", pass_pt_count, 100.0*(double)pass_pt_count/(double)all_count, 100.0*(double)pass_pt_count/(double)pass_lastGap_count);
  printf(" Evt_Nmu cut  | %13d | %10.2f | %10.2f\n", pass_Evt_Nmu_count, 100.0*(double)pass_Evt_Nmu_count/(double)all_count, 100.0*(double)pass_Evt_Nmu_count/(double)pass_pt_count);
  printf(" DG0 cut      | %13d | %10.2f | %10.2f\n", pass_DG0_count, 100.0*(double)pass_DG0_count/(double)all_count, 100.0*(double)pass_DG0_count/(double)pass_pt_count);
  printf(" DDG0 cut     | %13d | %10.2f | %10.2f\n", pass_DDG0_count, 100.0*(double)pass_DDG0_count/(double)all_count, 100.0*(double)pass_DDG0_count/(double)pass_DG0_count);
  printf(" chi2 cut     | %13d | %10.2f | %10.2f\n", pass_chi2_count, 100.0*(double)pass_chi2_count/(double)all_count, 100.0*(double)pass_chi2_count/(double)pass_DDG0_count);
  printf(" rpcdca cut   | %13d | %10.2f | %10.2f\n", pass_rpc_dca_count, 100.0*(double)pass_rpc_dca_count/(double)all_count, 100.0*(double)pass_rpc_dca_count/(double)pass_chi2_count);

  //write to output files
  //histoutfile->cd();
  //save_hists_phys();
  
  treeoutfile->cd();
  save_basic_tree();
  
  
  //clean-up TFile's:
  //delete histoutfile;
  delete treeoutfile;
  delete infile;
  
}//end main function


//function built from Francescas conditions for runnum-which_window correlation.
//failure modes:   -1=bad arm or rpc#   -2=bad runnum
Float_t get_rpc_time_window(int arm, int rpc, int run) {
  // definitions of the narrow time window cuts:
  //static int rpc1_narrow_low[2]={20,24};
  //static int rpc1_narrow_high[2]={27,32};
  //static int rpc3_narrow_low[2]={14,20};
  //static int rpc3_narrow_high[2]={21,26};
  
  Int_t which_window = -1;

  if(arm<0 || arm>1 || rpc<0 || rpc>1) {
    which_window=-1;
    return which_window;
  }

  if(rpc==0) {
    if(arm==0) {
      //RPC1 South
      if(run <=367239 ){
        which_window=1;
      }else if(run< 386773){
        which_window=0;
      }else if( (run >=386773 && run <=388052)|| (run>388548 &&
            run<=388745)  ||
          (run>389322 && run<=389436) || (run>390966 && run<=391100) ||
          (run>391173 && run<=391588) || (run>391817 && run<=391876) ||
          (run>392297 && run<=393167) || (run>395648 && run<=396440) ||
          (run>396768)){
        which_window=1;
      }else if( (run >388052 && run<=388548) || (run>388745 && run<=389322) ||
          (run > 389436 && run<= 390966) ||(run>391100 && run<=391173) ||
          (run>391588 && run<= 391817) || (run>391876 && run<=392297) ||
          (run>393167 && run<=395648) || (run>396440 && run<=396768)  ){
        which_window=0;
      } else {
        which_window=-2;
      }
    } else if(arm==1) {
      //RPC1 North
      if(run <=367737 ){
        // 26-32
        which_window=1;
      }else if(run< 386773){
        which_window=0;
      }else if( (run>386773 && run<=388700) ||(run>388840&& run<=391817) ||
          (run>393167 && run<=396440) || (run>397738)
          ){
        which_window=0;
      }else if( (run>388700 && run<=388840) || (run>391817 && run<=393167)||
          (run>396440 && run<= 397738)
          ){
        which_window=1;
      } else {
        which_window=-2;
      }
    }
  } else if(rpc==1) {
    if(arm==0) {
      //RPC3 South
      if(run <=365663  ||  (run >=368037  && run <=368798)){
        which_window=0;
      }else if( (run > 365663 && run <=367737) ||  (run >= 367921 && run
            <=367928)){
        which_window=1;
      }else if( (run > 386773 && run<=386844) || (run>388866 &&run<=391167) ||
          (run>391450 && run<=391863) || (run>393460 && run<=394072) ||
          (run>396280 &&run<= 396785) || (run>397317 && run<=397738)
          ){
        which_window=0;
      }else if( (run>386844 && run<=388866) || (run>391167 && run<=391450) ||
          (run>391863 && run<=393460) || (run>394072 && run<=396280)||
          (run>396785 && run<=397317) || (run>397738)
          ){
        which_window=1;
      }else if  (run ==367919){
        //rpc3_low[0] = 24;
        //rpc3_high[0] = 29;
        which_window=-2;
      } else {
        which_window=-2;
      }
    } else if(arm==1) {
      //RPC3 North
      if(run <=365152||  (run >=368037  && run <=368798)){
        which_window=0;
      }else if( (run >= 365313 && run <=367737) ||(run >= 367921 && run
            <=367928)  ){
        which_window=1;
      }else if( (run>386773 && run<= 387809) || (run>388866 && run<=388986) ||
          (run>389336 && run<=389339) || (run>394072 && run<=397738)
          ){
        which_window=0;
      }else if(  (run>387809 && run<= 388866) || (run>388986 && run<= 389336)||
          (run>389339 && run<=394072) || (run>397738)
          ){
        which_window=1;
      }else if  (run ==367919){
        //rpc3_low[1]  = 24;
        //rpc3_high[1]  = 29;
        which_window=-2;
      }else {
        which_window=-2;
      }
    }
  }

  return which_window;
}

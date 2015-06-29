#ifndef __phys_hist_definitions_h
#define __phys_hist_definitions_h

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>
#include <fstream>



//general variable initialization
//float PI = 3.141592;
float PI = TMath::Pi();
float rpc1timeshift = 6;
float rpc3timeshift = 6;
const int n_arms     = 2;		//arm
const int n_stations = 2;		//sta
const int n_octants  = 8;		//oct
const int n_halfocts = 2;		//hal
const int n_radsegments = 3;	//rad
const int n_paddles = 9;		//pad
const int n_strips   = 64;		//str
const int n_projections = 4;	//pro



const int t_gatewidth  = 4;		//timing selection window's half width
const int dca_cut_rpc1 = 15;   //maximum dca for a hit to qualify as "associated" for rpc 1
const int dca_cut_rpc3 = 15;   //maximum dca for a hit to qualify as "associated" for rpc 3
float rpc1_t_ratio_threshold = 3; //rpc timing peak to average t ratio used for excluding noisy channels
float rpc3_t_ratio_threshold = 2; //rpc timing peak to average t ratio used for excluding noisy channels

// mutr->rpc3 projection labels
//   [0]=MatchVtx
//   [1]=MatchSt1
//   [2]=MatchSt3
//   [3]=MatchMuID
const char arm_label[2] = {'S','N'};
const char * proj_label[4] = {"MatchVtx","MatchSt1","MatchSt3","MatchMuID"};
int used_proj_rpc1 = 1;
int used_proj_rpc3 = 2;

const int n_rpcdca_conditions = 3;

// physics related variables and Names
const int n_dists = 31;    //distributions
const int n_trig = 5;    //3 triggers? rpc3, rpc1 and only sg1
const char distchar[n_dists][64] = {
  "pT",//0
  "px",//1
  "py",//2
  "pz",//3
  "phi",//4
  "eta",//5
  "DG0",//6
  "DDG0",//7
  "DG4",//8
  "chi2",//9
  "DCA_z",//10
  "DCA_r",//11
  "dphi12",//12
  "dphi23",//13
  "dw23",//14
  "Rpc1dca",//15
  "Rpc1time",//16
  "Rpc3dca",//17
  "Rpc3time",//18
  "fvtx_dphi",//19
  "fvtx_dr",//20
  "fvtx_dtheta",//21
  "fvtx_dr_dtheta",//22
  "fvtx_cone",//23
  "Rpc1x",//24
  "Rpc1y",//25
  "Rpc3x",//26
  "Rpc3y",//27
  "fvtx_tracklcone",//28
  "Rpc1timewindow",//29
  "Rpc3timewindow"};//30
const char disttitles[n_dists][64] = {
  "pT [GeV]",//0
  "px",//1
  "py",//2
  "pz [GeV]",//3
  "#phi [rad]",//4
  "#eta ",//5
  "DG0 [cm]",//6
  "DDG0  [cm]",//7
  "DG4  [cm]",//8
  "#chi^{2}",//9
  "DCA_z [cm]",//10
  "DCA_r  [cm]",//11
  "dphi12",//12
  "#delta #phi_{23} [rad]",//13
  "dw{23} []",//14
  "Rpc1dca [cm]",//15
  "Rpc1time",//16
  "Rpc3dca [cm]",//17
  "Rpc3time",//18
  "fvtx_dphi [rad]",//19
  "fvtx_dr []",//20
  "fvtx_dtheta [rad]",//21
  "fvtx_dr_dtheta []",//22
  "fvtx_cone []",//23
  "Rpc1x",//24
  "Rpc1y",//25
  "Rpc3x",//26
  "Rpc3y",//27
  "fvtx_tracklcone",//28
  "Rpc1timewindow",//29
  "Rpc3timewindow"};//30
const float distmin[n_dists] = {
  0,//0
  -60,//1
  -60,//2
  -250,//3
  -3.1415,//4
  1.2,//5
  0.,//6
  0.,//7
  0.,//8
  0.,//9
  0.,//10
  0.,//11
  -0.1,//12
  -0.1,//13
  -0.3,//14
  0.,//15
  0.,//16
  0.,//17
  0.,//18
  0,//19
  0.,//20
  0.,//21
  0,//22
  0,//23
  -1000,//24
  -1000,//25
  -1000,//26
  -1000,//27
  0,//28
  0,//29
  0};//30
const float distmax[n_dists] = {
  60.,//0
  60,//1
  60,//2
  250.,//3
  3.1415,//4
  2.6,//5
  20.,//6
  9.,//7
  28.,//8
  20.,//9
  55.,//10
  30.,//11
  0.1,//12
  0.1,//13
  0.3,//14
  20.,//15
  44.,//16
  40.,//17
  44.,//18
  1.5,//19
  100.,//20
  1.5,//21
  150.,//22
  60,//23
  1000,//24
  1000,//25
  1000,//26
  1000,//27
  60,//28
  44,//29
  44,};//30

const int n_var_list = n_dists + 7;

char mu_bg_labels[9][20] = {"dy","light","onium","onlyz","openbottom","opencharm","whad","wtau","z"};

//Declare ifstream and ofstream
ifstream in;
ofstream out;

//first 2 is charge index, second 2 is fvtx bad/good
TH1F * h_phys_dists[n_arms][2][n_rpcdca_conditions][2][n_dists];
TH1F * h_fvtx_dists[n_arms][2][3];
TH1F * h_rpc_dists[n_arms][2][2];
TH2F * h2_dg0_ddg0[n_arms][2][n_rpcdca_conditions][2];
TH2F * h2_chi2_dcar[n_arms][2][n_rpcdca_conditions][2];
TH2F * h2_combined_dg0_ddg0[n_arms][2];
TH2F * h2_combined_chi2_dcar[n_arms][2];

TTree * basic_cuts_tree;
TTree * wness_tree;

void define_hists_phys(){
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      for (int rpcdca_condition = 0; rpcdca_condition < n_rpcdca_conditions; rpcdca_condition++)
      {
        char dhistoname[255];
        for (int good_fvtx=0; good_fvtx<2; good_fvtx++)
        {
          
          sprintf(dhistoname, "h2_dist_dg0_vs_ddg0_arm%d_charge%d_rpcdca%d_goodfvtx%d",arm,charge,rpcdca_condition,good_fvtx);
          h2_dg0_ddg0[arm][charge][rpcdca_condition][good_fvtx] = new TH2F(dhistoname,dhistoname,100,distmin[7],distmax[7],100,distmin[6],distmax[6]);
          
          sprintf(dhistoname, "h2_dist_chi2_dcar_arm%d_charge%d_rpcdca%d_goodfvtx%d",arm,charge,rpcdca_condition,good_fvtx);
          h2_chi2_dcar[arm][charge][rpcdca_condition][good_fvtx] = new TH2F(dhistoname,dhistoname,100,distmin[11],distmax[11],100,distmin[9],distmax[9]);
          
          for (int dists = 0; dists < n_dists; dists++)
          {
            sprintf(dhistoname, "h_dist_%s_arm%d_charge%d_rpcdca%d_goodfvtx%d",distchar[dists],arm,charge,rpcdca_condition,good_fvtx);
            h_phys_dists[arm][charge][rpcdca_condition][good_fvtx][dists] = new TH1F(dhistoname,dhistoname,100,distmin[dists],distmax[dists]);
          }
        }
      }
    }
  }

}

void define_hists_phys(TFile * infile){
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      for (int rpcdca_condition = 0; rpcdca_condition < n_rpcdca_conditions; rpcdca_condition++)
      {
        char dhistoname[255];
        for (int good_fvtx=0; good_fvtx<2; good_fvtx++)
        {
          sprintf(dhistoname, "h2_dist_dg0_vs_ddg0_arm%d_charge%d_rpcdca%d_goodfvtx%d",arm,charge,rpcdca_condition,good_fvtx);
          h2_dg0_ddg0[arm][charge][rpcdca_condition][good_fvtx] = (TH2F*)infile->Get(dhistoname);

          sprintf(dhistoname, "h2_dist_chi2_dcar_arm%d_charge%d_rpcdca%d_goodfvtx%d",arm,charge,rpcdca_condition,good_fvtx);
          h2_chi2_dcar[arm][charge][rpcdca_condition][good_fvtx] = (TH2F*)infile->Get(dhistoname);
          for (int dists = 0; dists < n_dists; dists++)
          {
            sprintf(dhistoname, "h_dist_%s_arm%d_charge%d_rpcdca%d_goodfvtx%d",distchar[dists],arm,charge,rpcdca_condition,good_fvtx);
            h_phys_dists[arm][charge][rpcdca_condition][good_fvtx][dists] = (TH1F*)infile->Get(dhistoname);
          }
        }
      }
    }
  }

}

void save_hists_phys() {
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      for (int rpcdca_condition = 0; rpcdca_condition < n_rpcdca_conditions; rpcdca_condition++)
      {
        for(int good_fvtx=0; good_fvtx<2; good_fvtx++)
        {
          h2_dg0_ddg0[arm][charge][rpcdca_condition][good_fvtx]->Write();

          h2_chi2_dcar[arm][charge][rpcdca_condition][good_fvtx]->Write();
          for (int dists = 0; dists < n_dists; dists++)
          {
            h_phys_dists[arm][charge][rpcdca_condition][good_fvtx][dists]->Write();
          }
        }
      }
    }
  }
}


void define_hists_wness(){
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      char dhistoname[255];

      sprintf(dhistoname, "h2_combined_dg0_vs_ddg0_arm%d_charge%d",arm,charge);
      h2_combined_dg0_ddg0[arm][charge] = new TH2F(dhistoname,dhistoname,100,distmin[7],distmax[7],100,distmin[6],distmax[6]);

      sprintf(dhistoname, "h2_combined_chi2_dcar_arm%d_charge%di",arm,charge);
      h2_combined_chi2_dcar[arm][charge] = new TH2F(dhistoname,dhistoname,100,distmin[11],distmax[11],100,distmin[9],distmax[9]);

      sprintf(dhistoname, "h_rpc1dca_arm%d_charge%d",arm,charge);
      h_rpc_dists[arm][charge][0] = new TH1F(dhistoname,dhistoname,100,distmin[15],distmax[15]);

      sprintf(dhistoname, "h_rpc3dca_arm%d_charge%d",arm,charge);
      h_rpc_dists[arm][charge][1] = new TH1F(dhistoname,dhistoname,100,distmin[17],distmax[17]);

      sprintf(dhistoname, "h_fvtx_dphi_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][0] = new TH1F(dhistoname,dhistoname,100,distmin[19],distmax[19]);

      sprintf(dhistoname, "h_fvtx_drxdtheta_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][1] = new TH1F(dhistoname,dhistoname,100,distmin[22],distmax[22]);

      sprintf(dhistoname, "h_fvtx_cone_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][2] = new TH1F(dhistoname,dhistoname,60,distmin[23],distmax[23]);

    }
  }

}

void define_hists_wness(TFile * infile){
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      char dhistoname[255];

      sprintf(dhistoname, "h2_combined_dg0_vs_ddg0_arm%d_charge%d",arm,charge);
      h2_combined_dg0_ddg0[arm][charge] = (TH2F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h2_combined_chi2_dcar_arm%d_charge%di",arm,charge);
      h2_combined_chi2_dcar[arm][charge] = (TH2F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_rpc1dca_arm%d_charge%d",arm,charge);
      h_rpc_dists[arm][charge][0] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_rpc3dca_arm%d_charge%d",arm,charge);
      h_rpc_dists[arm][charge][1] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_fvtx_dphi_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][0] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_fvtx_drxdtheta_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][1] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_fvtx_cone_arm%d_charge%d",arm,charge);
      h_fvtx_dists[arm][charge][2] = (TH1F*)infile->Get(dhistoname);


    }
  }

}

void save_hists_wness() {
  for (int arm = 0; arm < n_arms; arm ++ )
  {
    for (int charge = 0; charge < 2; charge++)
    {
      h2_combined_dg0_ddg0[arm][charge] ->Write();

      h2_combined_chi2_dcar[arm][charge]->Write();

      h_rpc_dists[arm][charge][0]->Write();

      h_rpc_dists[arm][charge][1]->Write();

      h_fvtx_dists[arm][charge][0]->Write();

      h_fvtx_dists[arm][charge][1]->Write();

      h_fvtx_dists[arm][charge][2]->Write();
    }

  }
}

void define_basic_tree() {
  
  basic_cuts_tree = new TTree("basic_cuts_tree","basic_cuts_tree");
  Int_t temp_int;
  Float_t temp_float;
  basic_cuts_tree->Branch("Run_Number", &temp_int, "Run_Number/I");
  basic_cuts_tree->Branch("Evt_Number", &temp_int, "Evt_Number/I");
  basic_cuts_tree->Branch("triggerbit", &temp_int, "triggerbit/I");
  basic_cuts_tree->Branch("Evt_bbcZ", &temp_float, "Evt_bbcZ/F");
  basic_cuts_tree->Branch("clockcross", &temp_int, "clockcross/I");
  basic_cuts_tree->Branch("Wness", &temp_float, "Wness/F");
  basic_cuts_tree->Branch("charge", &temp_float, "charge/F");
  basic_cuts_tree->Branch("pT", &temp_float, "pT/F");
  basic_cuts_tree->Branch("px", &temp_float, "px/F");
  basic_cuts_tree->Branch("py", &temp_float, "py/F");
  basic_cuts_tree->Branch("pz", &temp_float, "pz/F");
  basic_cuts_tree->Branch("phi", &temp_float, "phi/F");
  basic_cuts_tree->Branch("eta", &temp_float, "eta/F");
  basic_cuts_tree->Branch("DG0", &temp_float, "DG0/F");
  basic_cuts_tree->Branch("DDG0", &temp_float, "DDG0/F");
  basic_cuts_tree->Branch("DG4", &temp_float, "DG4/F");
  basic_cuts_tree->Branch("chi2", &temp_float, "chi2/F");
  basic_cuts_tree->Branch("DCA_z", &temp_float, "DCA_z/F");
  basic_cuts_tree->Branch("DCA_r", &temp_float, "DCA_r/F");
  basic_cuts_tree->Branch("dphi12", &temp_float, "dphi12/F");
  basic_cuts_tree->Branch("dphi23", &temp_float, "dphi23/F");
  basic_cuts_tree->Branch("dw23", &temp_float, "dw23/F");
  basic_cuts_tree->Branch("Rpc1dca", &temp_float, "Rpc1dca/F");
  basic_cuts_tree->Branch("Rpc1x", &temp_float, "Rpc1x/F");
  basic_cuts_tree->Branch("Rpc1y", &temp_float, "Rpc1y/F");
  basic_cuts_tree->Branch("Rpc1time", &temp_float, "Rpc1time/F");
  basic_cuts_tree->Branch("Rpc1timewindow", &temp_float, "Rpc1timewindow/F");
  basic_cuts_tree->Branch("Rpc3dca", &temp_float, "Rpc3dca/F");
  basic_cuts_tree->Branch("Rpc3x", &temp_float, "Rpc3x/F");
  basic_cuts_tree->Branch("Rpc3y", &temp_float, "Rpc3y/F");
  basic_cuts_tree->Branch("Rpc3time", &temp_float, "Rpc3time/F");
  basic_cuts_tree->Branch("Rpc3timewindow", &temp_float, "Rpc3timewindow/F");
  basic_cuts_tree->Branch("fvtx_dphi", &temp_float, "fvtx_dphi/F");
  basic_cuts_tree->Branch("fvtx_dr", &temp_float, "fvtx_dr/F");
  basic_cuts_tree->Branch("fvtx_dtheta", &temp_float, "fvtx_dtheta/F");
  basic_cuts_tree->Branch("fvtx_dr_dtheta", &temp_float, "fvtx_dr_dtheta/F");
  basic_cuts_tree->Branch("fvtx_cone", &temp_int, "fvtx_cone/I");
  basic_cuts_tree->Branch("fvtx_tracklcone", &temp_int, "fvtx_tracklcone/I");

  basic_cuts_tree->ResetBranchAddresses();
  
}

void define_basic_tree(TFile * infile) {
  basic_cuts_tree = (TTree*) infile->Get("basic_cuts_tree");
}

void save_basic_tree() {
  basic_cuts_tree->Write();
}

void define_wness_tree() {
  
  wness_tree = new TTree("wness_tree","wness_tree");
  Int_t temp_int;
  Float_t temp_float;
  wness_tree->Branch("Run_Number", &temp_int, "Run_Number/I");
  wness_tree->Branch("Evt_Number", &temp_int, "Evt_Number/I");
  wness_tree->Branch("triggerbit", &temp_int, "triggerbit/I");
  wness_tree->Branch("Evt_bbcZ", &temp_float, "Evt_bbcZ/F");
  wness_tree->Branch("clockcross", &temp_int, "clockcross/I");
  wness_tree->Branch("Wness", &temp_float, "Wness/F");
  wness_tree->Branch("charge", &temp_float, "charge/F");
  wness_tree->Branch("pT", &temp_float, "pT/F");
  wness_tree->Branch("px", &temp_float, "px/F");
  wness_tree->Branch("py", &temp_float, "py/F");
  wness_tree->Branch("pz", &temp_float, "pz/F");
  wness_tree->Branch("phi", &temp_float, "phi/F");
  wness_tree->Branch("eta", &temp_float, "eta/F");
  wness_tree->Branch("DG0", &temp_float, "DG0/F");
  wness_tree->Branch("DDG0", &temp_float, "DDG0/F");
  wness_tree->Branch("DG4", &temp_float, "DG4/F");
  wness_tree->Branch("chi2", &temp_float, "chi2/F");
  wness_tree->Branch("DCA_z", &temp_float, "DCA_z/F");
  wness_tree->Branch("DCA_r", &temp_float, "DCA_r/F");
  wness_tree->Branch("dphi12", &temp_float, "dphi12/F");
  wness_tree->Branch("dphi23", &temp_float, "dphi23/F");
  wness_tree->Branch("dw23", &temp_float, "dw23/F");
  wness_tree->Branch("Rpc1dca", &temp_float, "Rpc1dca/F");
  wness_tree->Branch("Rpc1x", &temp_float, "Rpc1x/F");
  wness_tree->Branch("Rpc1y", &temp_float, "Rpc1y/F");
  wness_tree->Branch("Rpc1time", &temp_float, "Rpc1time/F");
  wness_tree->Branch("Rpc1timewindow", &temp_float, "Rpc1timewindow/F");
  wness_tree->Branch("Rpc3dca", &temp_float, "Rpc3dca/F");
  wness_tree->Branch("Rpc3x", &temp_float, "Rpc3x/F");
  wness_tree->Branch("Rpc3y", &temp_float, "Rpc3y/F");
  wness_tree->Branch("Rpc3time", &temp_float, "Rpc3time/F");
  wness_tree->Branch("Rpc3timewindow", &temp_float, "Rpc3timewindow/F");
  wness_tree->Branch("fvtx_dphi", &temp_float, "fvtx_dphi/F");
  wness_tree->Branch("fvtx_dr", &temp_float, "fvtx_dr/F");
  wness_tree->Branch("fvtx_dtheta", &temp_float, "fvtx_dtheta/F");
  wness_tree->Branch("fvtx_dr_dtheta", &temp_float, "fvtx_dr_dtheta/F");
  wness_tree->Branch("fvtx_cone", &temp_int, "fvtx_cone/I");
  wness_tree->Branch("fvtx_tracklcone", &temp_int, "fvtx_tracklcone/I");

  wness_tree->ResetBranchAddresses();
  
}

void define_wness_tree(TFile * infile) {
  wness_tree = (TTree*) infile->Get("wness_tree");
}

void save_wness_tree() {
  wness_tree->Write();
}

#endif

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
const int n_dists = 19;    //distributions
const int n_trig = 5;    //3 triggers? rpc3, rpc1 and only sg1
const char distchar[n_dists][64] = {"pT","pz","phi","eta",
  "DG0","DDG0","DG4","chi2",
  "DCA_z","DCA_r","dphi23","dw23",
  "Rpc1dca","Rpc3dca",
  "fvtx_dphi","fvtx_dr","fvtx_dtheta","fvtx_dr_dtheta",
  "fvtx_cone"};
const char disttitles[n_dists][64] = {"pT [GeV]","pz [GeV]","#phi [rad]","#eta ",
  "DG0 [cm]","DDG0  [cm]","DG4  [cm]","#chi^{2}",
  "DCA_z [cm]","DCA_r  [cm]","#delta #phi_{23} [rad]","dw{23} []",
  "Rpc1dca [cm]","Rpc3dca [cm]",
  "fvtx_dphi [rad]","fvtx_dr []","fvtx_dtheta [rad]","fvtx_dr_dtheta []",
  "fvtx_cone []"};
const float distmin[n_dists] = {0,-250,-3.1415,1.2,
  0.,0.,0.,0.,
  0.,0.,-0.1,-0.3,
  0.,0.,
  0,0.,0.,0,
  0};
const float distmax[n_dists] = {60.,250.,3.1415,2.6,
  20.,9.,28.,20.,
  55.,30.,0.1,0.3,
  20.,40.,
  1.5,100.,1.5,150.,
  60};

const int n_var_list = n_dists + 7;

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

TNtuple * basic_cuts_ntuple;
TNtuple * wness_ntuple;

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
          h2_dg0_ddg0[arm][charge][rpcdca_condition][good_fvtx] = new TH2F(dhistoname,dhistoname,100,distmin[5],distmax[5],100,distmin[4],distmax[4]);
          
          sprintf(dhistoname, "h2_dist_chi2_dcar_arm%d_charge%d_rpcdca%d_goodfvtx%d",arm,charge,rpcdca_condition,good_fvtx);
          h2_chi2_dcar[arm][charge][rpcdca_condition][good_fvtx] = new TH2F(dhistoname,dhistoname,100,distmin[9],distmax[9],100,distmin[7],distmax[7]);
          
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
      h2_combined_dg0_ddg0[arm][charge] = new TH2F(dhistoname,dhistoname,100,distmin[5],distmax[5],100,distmin[4],distmax[4]);

      sprintf(dhistoname, "h2_combined_chi2_dcar_arm%d_charge%di",arm,charge);
      h2_combined_chi2_dcar[arm][charge] = new TH2F(dhistoname,dhistoname,100,distmin[9],distmax[9],100,distmin[7],distmax[7]);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[12],arm,charge);
      h_rpc_dists[arm][charge][0] = new TH1F(dhistoname,dhistoname,100,distmin[12],distmax[12]);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[13],arm,charge);
      h_rpc_dists[arm][charge][1] = new TH1F(dhistoname,dhistoname,100,distmin[13],distmax[13]);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[14],arm,charge);
      h_fvtx_dists[arm][charge][0] = new TH1F(dhistoname,dhistoname,100,distmin[14],distmax[14]);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[17],arm,charge);
      h_fvtx_dists[arm][charge][1] = new TH1F(dhistoname,dhistoname,100,distmin[17],distmax[17]);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[18],arm,charge);
      h_fvtx_dists[arm][charge][2] = new TH1F(dhistoname,dhistoname,100,distmin[18],distmax[18]);

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

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[12],arm,charge);
      h_rpc_dists[arm][charge][0] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[13],arm,charge);
      h_rpc_dists[arm][charge][1] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[14],arm,charge);
      h_fvtx_dists[arm][charge][0] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[17],arm,charge);
      h_fvtx_dists[arm][charge][1] = (TH1F*)infile->Get(dhistoname);

      sprintf(dhistoname, "h_%s_arm%d_charge%d",distchar[18],arm,charge);
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

void define_basic_ntuple() {
  char var_name[n_var_list][50];
  sprintf(var_name[0],"Run_Number");
  sprintf(var_name[1],"Evt_Number");
  sprintf(var_name[2],"triggerbit");
  sprintf(var_name[3],"Evt_bbcZ");
  sprintf(var_name[4],"SpinX_ID");
  sprintf(var_name[5],"Wness");
  sprintf(var_name[6],"charge");
  for(int dist=0; dist<n_dists; dist++) {
    sprintf(var_name[7+dist],"%s",distchar[dist]);
  }
  char var_list[500];
  for(int var=0; var<n_var_list; var++) {
    if(var==0) {
      sprintf(var_list,"%s",var_name[var]);
    } else {
      sprintf(var_list,"%s:%s",var_list,var_name[var]);
    }
  }
  
  basic_cuts_ntuple = new TNtuple("basic_cuts_ntuple","basic_cuts_ntuple",var_list);
}

void define_basic_ntuple(TFile * infile) {
  basic_cuts_ntuple = (TNtuple*) infile->Get("basic_cuts_ntuple");
}

void save_basic_ntuple() {
  basic_cuts_ntuple->Write();
}

void define_wness_ntuple() {
  char var_name[n_var_list][50];
  sprintf(var_name[0],"Run_Number");
  sprintf(var_name[1],"Evt_Number");
  sprintf(var_name[2],"triggerbit");
  sprintf(var_name[3],"Evt_bbcZ");
  sprintf(var_name[4],"SpinX_ID");
  sprintf(var_name[5],"Wness");
  sprintf(var_name[6],"charge");
  for(int dist=0; dist<n_dists; dist++) {
    sprintf(var_name[7+dist],"%s",distchar[dist]);
  }
  char var_list[500];
  for(int var=0; var<n_var_list; var++) {
    if(var==0) {
      sprintf(var_list,"%s",var_name[var]);
    } else {
      sprintf(var_list,"%s:%s",var_list,var_name[var]);
    }
  }
  
  wness_ntuple = new TNtuple("wness_ntuple","wness_ntuple",var_list);
}

void define_wness_ntuple(TFile * infile) {
  wness_ntuple = (TNtuple*) infile->Get("wness_ntuple");
//  wness_ntuple = (TTree*) infile->Get("newsngmuons_basic_cut");
}

void save_wness_ntuple() {
  wness_ntuple->Write();
}

#endif

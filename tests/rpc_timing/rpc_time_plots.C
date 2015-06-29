#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include <string>
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/phys_hist_definitions.h"

using namespace std;

void rpc_time_plots() {
  TFile *data = new TFile("/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_total_wness_refrun367593.root");
  TTree *temp_tree = (TTree*)data->Get("wness_tree");
  TTree *data_tree = (TTree*)temp_tree->Clone("data_tree");

  TFile *sim = new TFile("/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/current_wness_tree_sig_sim.root");
  temp_tree = (TTree*)sim->Get("wness_tree");
  TTree *sim_tree = (TTree*)temp_tree->Clone("sim_tree");
  
  TH1F *h_sim[2];
  TH1F *h_data_low[2];
  TH1F *h_data_hi[2];
  
  h_sim[0] = new TH1F("h_sim_rpc1","h_sim_rpc1",44,0,44);
  sim_tree->Project("h_sim_rpc1","Rpc1time","Rpc1time<100");
  h_data_low[0] = new TH1F("h_data_low_rpc1","h_data_low_rpc1",44,0,44);
  data_tree->Project("h_data_low_rpc1","Rpc1time","Rpc1time<100 && Wness<.92");
  h_data_hi[0] = new TH1F("h_data_hi_rpc1","h_data_hi_rpc1",44,0,44);
  data_tree->Project("h_data_hi_rpc1","Rpc1time","Rpc1time<100 && Wness>.92");

  h_sim[1] = new TH1F("h_sim_rpc3","h_sim_rpc3",44,0,44);
  sim_tree->Project("h_sim_rpc3","Rpc3time","Rpc3time<100");
  h_data_low[1] = new TH1F("h_data_low_rpc3","h_data_low_rpc3",44,0,44);
  data_tree->Project("h_data_low_rpc3","Rpc3time","Rpc3time<100 && Wness<.92");
  h_data_hi[1] = new TH1F("h_data_hi_rpc3","h_data_hi_rpc3",44,0,44);
  data_tree->Project("h_data_hi_rpc3","Rpc3time","Rpc3time<100 && Wness>.92");

  TCanvas *c1;
  c1= new TCanvas("Wness","Wness",800,600);
  
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  
  double data_low_integral[2];
  double data_hi_integral[2];
  double sim_integral[2];

  //data_low_integral[0] = h_data_low[0]->Integral("width");
  //data_hi_integral[0] = h_data_hi[0]->Integral("width");
  //sim_integral[0] = h_sim[0]->Integral("width");
  data_low_integral[0] = h_data_low[0]->GetMaximum();
  data_hi_integral[0] = h_data_hi[0]->GetMaximum();
  sim_integral[0] = h_sim[0]->GetMaximum();

  h_data_low[0]->Scale(1.0/data_low_integral[0]);
  h_data_hi[0]->Scale(1.0/data_hi_integral[0]);
  h_sim[0]->Scale(1.0/sim_integral[0]);

  h_data_low[0]->SetLineColor(kRed);
  h_data_hi[0]->SetLineColor(kBlue);
  h_sim[0]->SetLineColor(kBlack);
  h_sim[0]->SetTitle("");

  h_sim[0]->GetXaxis()->SetTitle("RpcTimeBin");
  h_sim[0]->GetXaxis()->SetTitleSize(.06);
  h_sim[0]->GetXaxis()->SetTitleOffset(0.85);
  h_sim[0]->GetXaxis()->SetLabelSize(.05);
  h_sim[0]->GetYaxis()->SetTitle("Normalized Yields");
  h_sim[0]->GetYaxis()->SetTitleSize(.06);
  h_sim[0]->GetYaxis()->SetTitleOffset(0.85);
  h_sim[0]->GetYaxis()->SetLabelSize(.05);
  

  h_sim[0]->SetTitle("RPC1");
  h_sim[0]->Draw();
  h_data_low[0]->Draw("SAME");
  h_data_hi[0]->Draw("SAME");

  TLegend * leg = new TLegend(0.5,0.75,0.9,0.9);
  leg->AddEntry(h_data_low[0],"Data < .92 Wness");
  leg->AddEntry(h_data_hi[0],"Data > .92 Wness");
  leg->AddEntry(h_sim[0],"W->#mu Sim");
  leg->Draw();

  TCanvas *c2;
  c2= new TCanvas("rpc3","rpc3",800,600);
  //data_low_integral[1] = h_data_low[1]->Integral("width");
  //data_hi_integral[1] = h_data_hi[1]->Integral("width");
  //sim_integral[1] = h_sim[1]->Integral("width");
  data_low_integral[1] = h_data_low[1]->GetMaximum();
  data_hi_integral[1] = h_data_hi[1]->GetMaximum();
  sim_integral[1] = h_sim[1]->GetMaximum();

  h_data_low[1]->Scale(1.0/data_low_integral[1]);
  h_data_hi[1]->Scale(1.0/data_hi_integral[1]);
  h_sim[1]->Scale(1.0/sim_integral[1]);

  h_data_low[1]->SetLineColor(kRed);
  h_data_hi[1]->SetLineColor(kBlue);
  h_sim[1]->SetLineColor(kBlack);
  h_sim[1]->SetTitle("");

  h_sim[1]->GetXaxis()->SetTitle("RpcTimeBin");
  h_sim[1]->GetXaxis()->SetTitleSize(.06);
  h_sim[1]->GetXaxis()->SetTitleOffset(0.85);
  h_sim[1]->GetXaxis()->SetLabelSize(.05);
  h_sim[1]->GetYaxis()->SetTitle("Normalized Yields");
  h_sim[1]->GetYaxis()->SetTitleSize(.06);
  h_sim[1]->GetYaxis()->SetTitleOffset(0.85);
  h_sim[1]->GetYaxis()->SetLabelSize(.05);
  

  h_sim[1]->SetTitle("RPC3");
  h_sim[1]->Draw();
  h_data_low[1]->Draw("SAME");
  h_data_hi[1]->Draw("SAME");

  TLegend * leg2 = new TLegend(0.5,0.75,0.9,0.9);
  leg2->AddEntry(h_data_low[1],"Data < .92 Wness");
  leg2->AddEntry(h_data_hi[1],"Data > .92 Wness");
  leg2->AddEntry(h_sim[1],"W->#mu Sim");
  leg2->Draw();

  Int_t
    Run_Number,  Evt_Number,  triggerbit,
    clockcross,  fvtx_cone, fvtx_tracklcone;

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


  data_tree->SetBranchAddress("Wness",&Wness);
  data_tree->SetBranchAddress("px",&px);
  data_tree->SetBranchAddress("py",&py);
  data_tree->SetBranchAddress("dphi12",&dphi12);
  data_tree->SetBranchAddress("Rpc1time",&Rpc1time);
  data_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
  data_tree->SetBranchAddress("Rpc1x",&Rpc1x);
  data_tree->SetBranchAddress("Rpc1y",&Rpc1y);
  data_tree->SetBranchAddress("Rpc3time",&Rpc3time);
  data_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
  data_tree->SetBranchAddress("Rpc3x",&Rpc3x);
  data_tree->SetBranchAddress("Rpc3y",&Rpc3y);
  data_tree->SetBranchAddress("Run_Number", &Run_Number);
  data_tree->SetBranchAddress("Evt_Number", &Evt_Number);
  data_tree->SetBranchAddress("triggerbit", &triggerbit);
  data_tree->SetBranchAddress("Evt_bbcZ", &Evt_bbcZ);
  data_tree->SetBranchAddress("clockcross", &clockcross);
  data_tree->SetBranchAddress("data", &data);
  data_tree->SetBranchAddress("charge", &charge);
  data_tree->SetBranchAddress("pT", &pT);
  data_tree->SetBranchAddress("pz", &pz);
  data_tree->SetBranchAddress("phi", &phi);
  data_tree->SetBranchAddress("eta", &eta);
  data_tree->SetBranchAddress("DG0", &DG0);
  data_tree->SetBranchAddress("DDG0", &DDG0);
  data_tree->SetBranchAddress("DG4", &DG4);
  data_tree->SetBranchAddress("chi2", &chi2);
  data_tree->SetBranchAddress("DCA_z", &DCA_z);
  data_tree->SetBranchAddress("DCA_r", &DCA_r);
  data_tree->SetBranchAddress("dphi23", &dphi23);
  data_tree->SetBranchAddress("dw23", &dw23);
  data_tree->SetBranchAddress("Rpc1dca", &Rpc1dca);
  data_tree->SetBranchAddress("Rpc3dca", &Rpc3dca);
  data_tree->SetBranchAddress("fvtx_dphi", &fvtx_dphi);
  data_tree->SetBranchAddress("fvtx_dr", &fvtx_dr);
  data_tree->SetBranchAddress("fvtx_dtheta", &fvtx_dtheta);
  data_tree->SetBranchAddress("fvtx_dr_dtheta", &fvtx_dr_dtheta);
  data_tree->SetBranchAddress("fvtx_cone", &fvtx_cone);
  data_tree->SetBranchAddress("fvtx_tracklcone", &fvtx_tracklcone);


  TFile *outfile = new TFile("rpc_t_diff.root","recreate");

  TH1F *h_rpc_t_diff[2][4][3];
  TH2F *h2_rpc_t_diff[2][3][4];
  TH1F *h_rpc1_t[2][4][3];
  TH1F *h_rpc3_t[2][4][3];

  //TH2F *h_delta_vs_rpcdca[2];
  //h_delta_vs_rpcdca[0] = new TH2F("h_deltat_vs_rpcdca1","h_deltat_vs_rpcdca1",80,0,20,44,-22,22);
  //h_delta_vs_rpcdca[0] = new TH2F("h_deltat_vs_rpcdca3","h_deltat_vs_rpcdca3",80,0,40,44,-22,22);
  char *rpc_cond[4] = {"r1only","r3only","r1and3","rall"};
  char *vs_var[4] = {"DG0","DDG0","chi2","DCA_r"};
  int var_num[4] = {6,7,9,11};
  char *wness_cond[3] = {"wlt92","wgt92","wall"};
  char name[50];
  for(int f=0; f<2; f++) {
    for(int wc=0; wc<3; wc++) {
      for(int rc=0; rc<4; rc++) {
        sprintf(name,"h_%s_rpc_t_diff_%s_%s",(f==0)?"data":"wsim",rpc_cond[rc],wness_cond[wc]);
        h_rpc_t_diff[f][rc][wc] = (TH1F*) new TH1F(name,name,44,-22,22);
        sprintf(name,"h_%s_rpc1_t_%s_%s",(f==0)?"data":"wsim",rpc_cond[rc],wness_cond[wc]);
        h_rpc1_t[f][rc][wc] = (TH1F*) new TH1F(name,name,44,0,44);
        sprintf(name,"h_%s_rpc3_t_%s_%s",(f==0)?"data":"wsim",rpc_cond[rc],wness_cond[wc]);
        h_rpc3_t[f][rc][wc] = (TH1F*) new TH1F(name,name,44,0,44);
        
        sprintf(name,"h2_%s_rpc_t_diff_vs_%s_%s",(f==0)?"data":"wsim",vs_var[rc],wness_cond[wc]);
        h2_rpc_t_diff[f][wc][rc] = (TH2F*) new TH2F(name,name,44,-22,22,distmax[var_num[rc]]*2,distmin[var_num[rc]],distmax[var_num[rc]]);
      }
    }

    int entries = (f==0)?data_tree->GetEntries():sim_tree->GetEntries();
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

      float corrected_rpc1_time;
      float corrected_rpc3_time;

      if(f==0) {
        data_tree->GetEntry(i);

        float rpc1correction,rpc3correction;
        if(Rpc1timewindow<0) {
          rpc1correction = 0;
        } else {
          rpc1correction = rpc1timeshift*Rpc1timewindow;
        }

        if(Rpc3timewindow<0) {
          rpc3correction = 0;
        } else {
          rpc3correction = rpc3timeshift*Rpc3timewindow;
        }

        corrected_rpc1_time = Rpc1time - rpc1correction;
        corrected_rpc3_time = Rpc3time - rpc3correction;
      } else {
        sim_tree->GetEntry(i);

        corrected_rpc1_time = Rpc1time;
        corrected_rpc3_time = Rpc3time-27;
      }

      int rpcdca_condition = -1;
      if((Rpc1dca < 100) && (Rpc3dca < 100)) {
        rpcdca_condition = 2;
      } else if(Rpc1dca < 100) {
        rpcdca_condition = 0;
      } else if(Rpc3dca < 100) {
        rpcdca_condition = 1;
      } else {
        //cout << "continuing" << endl;
        continue;
      }

      int wness_condition = (Wness<.92)?0:1;

      float rpc_t_diff = corrected_rpc3_time - corrected_rpc1_time;
      //if(f==1) cout<< rpc_t_diff << endl;
      h_rpc_t_diff[f][rpcdca_condition][wness_condition]->Fill(rpc_t_diff);
      h_rpc_t_diff[f][3][wness_condition]->Fill(rpc_t_diff);
      h_rpc_t_diff[f][rpcdca_condition][2]->Fill(rpc_t_diff);
      h_rpc_t_diff[f][3][2]->Fill(rpc_t_diff);

      h2_rpc_t_diff[f][wness_condition][0]->Fill(rpc_t_diff,DG0);
      h2_rpc_t_diff[f][2][0]->Fill(rpc_t_diff,DG0);
      h2_rpc_t_diff[f][wness_condition][1]->Fill(rpc_t_diff,DDG0);
      h2_rpc_t_diff[f][2][1]->Fill(rpc_t_diff,DDG0);
      h2_rpc_t_diff[f][wness_condition][2]->Fill(rpc_t_diff,chi2);
      h2_rpc_t_diff[f][2][2]->Fill(rpc_t_diff,chi2);
      h2_rpc_t_diff[f][wness_condition][3]->Fill(rpc_t_diff,DCA_r);
      h2_rpc_t_diff[f][2][3]->Fill(rpc_t_diff,DCA_r);

      h_rpc1_t[f][rpcdca_condition][wness_condition]->Fill(corrected_rpc1_time);
      h_rpc1_t[f][3][wness_condition]->Fill(corrected_rpc1_time);
      h_rpc1_t[f][rpcdca_condition][2]->Fill(corrected_rpc1_time);
      h_rpc1_t[f][3][2]->Fill(corrected_rpc1_time);

      h_rpc3_t[f][rpcdca_condition][wness_condition]->Fill(corrected_rpc3_time);
      h_rpc3_t[f][3][wness_condition]->Fill(corrected_rpc3_time);
      h_rpc3_t[f][rpcdca_condition][2]->Fill(corrected_rpc3_time);
      h_rpc3_t[f][3][2]->Fill(corrected_rpc3_time);

    }
  }

  TCanvas *c3;
  c3 = (TCanvas*) new TCanvas("c3","c3",800,600);
  h_rpc_t_diff[0][2][0]->SetLineColor(kRed);
  h_rpc_t_diff[0][2][0]->DrawCopy();
  h_rpc_t_diff[0][2][1]->SetLineColor(kBlue);
  h_rpc_t_diff[0][2][1]->DrawCopy("same");

  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  h_rpc_t_diff[0][2][0]->SetLineColor(kRed);
  float scale = h_rpc_t_diff[0][2][0]->GetMaximum();
  h_rpc_t_diff[0][2][0]->Scale(1.0/scale);
  h_rpc_t_diff[0][2][0]->DrawCopy();
  h_rpc_t_diff[0][2][1]->SetLineColor(kBlue);
  scale = h_rpc_t_diff[0][2][1]->GetMaximum();
  h_rpc_t_diff[0][2][1]->Scale(1.0/scale);
  h_rpc_t_diff[0][2][1]->DrawCopy("same");

  TCanvas *c5;
  c5 = (TCanvas*) new TCanvas("c5","c5",800,600);
  h_rpc1_t[0][3][0]->SetLineColor(kRed);
  h_rpc1_t[0][3][0]->DrawCopy();
  h_rpc1_t[0][3][1]->SetLineColor(kBlue);
  h_rpc1_t[0][3][1]->DrawCopy("same");
  h_rpc1_t[1][3][2]->SetLineColor(kBlack);
  h_rpc1_t[1][3][2]->DrawCopy("same");
  TLegend * leg3 = new TLegend(0.5,0.75,0.9,0.9);
  leg3->AddEntry(h_rpc1_t[0][3][0],"Data < .92 Wness");
  leg3->AddEntry(h_rpc1_t[0][3][1],"Data > .92 Wness");
  leg3->AddEntry(h_rpc1_t[1][3][2],"W->#mu Sim");
  leg3->Draw();

  TCanvas *c6;
  c6 = (TCanvas*) new TCanvas("c6","c6",800,800);
  c6->Divide(2,2);
  c6->cd(1);
  h2_rpc_t_diff[0][0][0]->DrawCopy();
  h2_rpc_t_diff[0][1][0]->SetMarkerColor(kRed);
  h2_rpc_t_diff[0][1][0]->DrawCopy("same");
  c6->cd(2);
  h2_rpc_t_diff[0][0][1]->DrawCopy();
  h2_rpc_t_diff[0][1][1]->SetMarkerColor(kRed);
  h2_rpc_t_diff[0][1][1]->DrawCopy("same");
  c6->cd(3);
  h2_rpc_t_diff[0][0][2]->DrawCopy();
  h2_rpc_t_diff[0][1][2]->SetMarkerColor(kRed);
  h2_rpc_t_diff[0][1][2]->DrawCopy("same");
  c6->cd(4);
  h2_rpc_t_diff[0][0][3]->DrawCopy();
  h2_rpc_t_diff[0][1][3]->SetMarkerColor(kRed);
  h2_rpc_t_diff[0][1][3]->DrawCopy("same");



  outfile->Write();
  outfile->Close();
  
}


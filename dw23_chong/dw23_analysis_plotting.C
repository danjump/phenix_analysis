#include <TF2.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>

void dw23_analysis_plotting() {

  TFile * f_sbg_had_dw23 = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_chong/output.root","READ");
  TCanvas * c_target_dw23 = new TCanvas("c_target_dw23","c_target_dw23",800,600);
  c_target_dw23->Divide(2,2);
  TCanvas * c_dw23_gpr_slices[2][2];


  TF1 * f_double_gaus_coax_target[2][2];
  TH1F * h_double_gaus_coax_target[2][2];
  TF1 * f_double_gaus_offset_target[2][2];
  TH1F * h_double_gaus_offset_target[2][2];
  TH1D * h_gpr_target[2][2]; 
  TH1D * h_gpr_slices[2][2][8]; 
  TH2F * h2_gpr_dw23_vs_wness[2][2];
  TH2F * h2_dw23_vs_wness[2][2];
  TH1D * h_dw23_slices[2][2][8];

  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      c_dw23_gpr_slices[a][c] = new TCanvas(Form("c_dw23_gpr_slices_a%d_c%d",a,c),Form("c_dw23_gpr_slices_a%d_c%d",a,c),1100,600);
      c_dw23_gpr_slices[a][c]->Divide(4,2);
      //sprintf(name,"h_sbg_fit_dw23_source0_a%d_c%d",arm,charge);
      //h_sbg_fit_dw23[a][c] = new TH1F(
      
      c_target_dw23->cd(1+a*2+c);
      f_double_gaus_coax_target[a][c] = (TF1*)f_sbg_had_dw23->Get(Form("f_double_gaus_coax_target_a%d_c%d",a,c)); 
      h_double_gaus_coax_target[a][c] = new TH1F(Form("h_double_gaus_coax_target_a%d_c%d",a,c),Form("h_double_gaus_coax_target_a%d_c%d",a,c),120,-.3,.3);
      h_double_gaus_coax_target[a][c]->FillRandom(f_double_gaus_coax_target[a][c]->GetName(),100000); 
      h_double_gaus_coax_target[a][c]->Scale(1.0/h_double_gaus_coax_target[a][c]->Integral("width"));
      h_double_gaus_coax_target[a][c]->SetLineColor(1); 
      h_double_gaus_coax_target[a][c]->SetTitle(".92 < Wness < 1");
      h_double_gaus_coax_target[a][c]->Draw();
      
      f_double_gaus_offset_target[a][c] = (TF1*)f_sbg_had_dw23->Get(Form("f_double_gaus_offset_target_a%d_c%d",a,c)); 
      h_double_gaus_offset_target[a][c] = new TH1F(Form("h_double_gaus_offset_target_a%d_c%d",a,c),Form("h_double_gaus_offset_target_a%d_c%d",a,c),120,-.3,.3);
      h_double_gaus_offset_target[a][c]->FillRandom(f_double_gaus_offset_target[a][c]->GetName(),100000); 
      h_double_gaus_offset_target[a][c]->Scale(1.0/h_double_gaus_offset_target[a][c]->Integral("width"));
      h_double_gaus_offset_target[a][c]->SetLineColor(2); 
      h_double_gaus_offset_target[a][c]->Draw("same");
     
      //c_target_dw23_gpr->cd(1+a*2+c);
      h2_gpr_dw23_vs_wness[a][c] = (TH2F*)f_sbg_had_dw23->Get(Form("h2_gpr_dw23_vs_wness_a%d_c%d",a,c)); 
      int nwness = h2_gpr_dw23_vs_wness[a][c]->GetYaxis()->GetNbins();
      h_gpr_target[a][c] = h2_gpr_dw23_vs_wness[a][c]->ProjectionX(Form("h_gpr_target_a%d_c%d",a,c),nwness*.92,nwness);
      h_gpr_target[a][c]->Scale(1.0/h_gpr_target[a][c]->Integral("width"));
      h_gpr_target[a][c]->SetLineColor(4); 
      h_gpr_target[a][c]->DrawCopy("same hist"); 
      h_gpr_target[a][c]->SetFillColor(4); 
      h_gpr_target[a][c]->SetFillStyle(3002); 
      h_gpr_target[a][c]->Draw("same e3"); 

      h2_dw23_vs_wness[a][c] = (TH2F*)f_sbg_had_dw23->Get(Form("h2_dw23_vs_wness_data_a%d_c%d",a,c)); 
      nwness = h2_dw23_vs_wness[a][c]->GetYaxis()->GetNbins();

      for(int i=1; i<9; i++) {
        c_dw23_gpr_slices[a][c]->cd(i);
        h_dw23_slices[a][c][i-1] = h2_dw23_vs_wness[a][c]->ProjectionY(Form("h_dw23_slices_a%d_c%d_s%d",a,c,i),nwness*(i/10.0),nwness*((i+1.0)/10.0));
        h_dw23_slices[a][c][i-1]->Scale(120/30*200/50*4*4);
        h_dw23_slices[a][c][i-1]->SetLineColor(2);
        h_dw23_slices[a][c][i-1]->SetTitle(Form("%.1f < Wness < %.1f",i/10.,i/10.+.1));
        h_dw23_slices[a][c][i-1]->Draw("");
        h_gpr_slices[a][c][i-1] = h2_gpr_dw23_vs_wness[a][c]->ProjectionX(Form("h_gpr_slices_a%d_c%d_s%d",a,c,i),nwness*(i/10.0),nwness*((i+1.0)/10.0));
        h_gpr_slices[a][c][i-1]->SetLineColor(4);
        h_gpr_slices[a][c][i-1]->DrawCopy("hist same");
        h_gpr_slices[a][c][i-1]->SetFillColor(4);
        h_gpr_slices[a][c][i-1]->SetFillStyle(3002);
        h_gpr_slices[a][c][i-1]->Draw("e3 same");
      }
    }
  }
  //f_sbg_had_dw23->Close();

}

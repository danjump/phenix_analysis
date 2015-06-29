#include <TCanvas.h>
#include <TFile.h>
#include <TNtuple.h>

void plot_wness(char * sig_infilename = (char*)"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_wness_refrun367593.root", 
    char * bkg_infilename = (char*)"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/backup/run13_data_wness_refrun367593.root") {
  
  TFile * sig_infile = new TFile(sig_infilename);
  TFile * bkg_infile = new TFile(bkg_infilename);

  TCanvas * c_wness = new TCanvas("c_wness","c_wness",1000,1000);
  c_wness->Divide(2,2);

  TNtuple * wness_ntuple;
  
  sig_infile->cd();
  wness_ntuple = (TNtuple*)sig_infile->Get("wness_ntuple");

  c_wness->cd(1);
  gPad->SetLogy();
  test->Draw("Wness","pz<0 && charge<0");

  c_wness->cd(2);
  gPad->SetLogy();
  test->Draw("Wness","pz<0 && charge>0");
  
  c_wness->cd(3);
  gPad->SetLogy();
  test->Draw("Wness","pz>0 && charge<0");
  
  c_wness->cd(4);
  gPad->SetLogy();
  test->Draw("Wness","pz>0 && charge>0");

  bkg_infile->cd();
  wness_ntuple = (TNtuple*)bkg_infile->Get("wness_ntuple");

  c_wness->cd(1);
  gPad->SetLogy();
  wness_ntuple->Draw("Wness","pz<0 && charge<0","same");

  c_wness->cd(2);
  gPad->SetLogy();
  wness_ntuple->Draw("Wness","pz<0 && charge>0","same");
  
  c_wness->cd(3);
  gPad->SetLogy();
  wness_ntuple->Draw("Wness","pz>0 && charge<0","same");
  
  c_wness->cd(4);
  gPad->SetLogy();
  wness_ntuple->Draw("Wness","pz>0 && charge>0","same");
}


#include "../../get_phys_dists/phys_hist_definitions.h"
#include "../../get_phys_dists/version_number.h"
#include <TCanvas.h>

void plot_wness_pdf_hists(char * sig_histinfilename, char * bkg_histinfilename) {
  
  TFile * histinfile[2];
  histinfile[0] = new TFile(sig_histinfilename);
  histinfile[1] = new TFile(bkg_histinfilename);
  
  //[arm][charge][goodfvtx]
  TCanvas * distributions[8][2];

  char name[500];
  for(int file=0; file<2; file++) {
    for(int i=0; i<8; i++) {
      sprintf(name,"c_file%d_dist%d",file,i);
      if(i>=2 && i<=4) {
        distributions[i][file] = new TCanvas(name,name,1000,1600);
        distributions[i][file]->Divide(2,4);
      } else {
        distributions[i][file] = new TCanvas(name,name,1000,800);
        distributions[i][file]->Divide(2,2);
      }
    }
  }

  for(int file = 0; file<2; file++) {
    define_hists_wness(histinfile[file]);

    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        distributions[0][file]->cd(1+arm+2*charge);
        h2_combined_dg0_ddg0[arm][charge]->Draw("COLZ");

        distributions[1][file]->cd(1+arm+2*charge);
        h2_combined_chi2_dcar[arm][charge]->Draw("COLZ");

        distributions[2][file]->cd(1+arm+2*charge);
        h_rpc_dists[arm][charge][0][0]->Draw();

        distributions[3][file]->cd(1+arm+2*charge);
        h_rpc_dists[arm][charge][1][1]->Draw();

        distributions[4][file]->cd(1+arm+2*charge);
        h_rpc_dists[arm][charge][2][0]->Draw();
        distributions[4][file]->cd(1+4+arm+2*charge);
        h_rpc_dists[arm][charge][2][1]->Draw();

        distributions[5][file]->cd(1+arm+2*charge);
        h_fvtx_dists[arm][charge][0]->Draw();
        distributions[6][file]->cd(1+arm+2*charge);
        h_fvtx_dists[arm][charge][1]->Draw();
        distributions[7][file]->cd(1+arm+2*charge);
        h_fvtx_dists[arm][charge][2]->Draw();
      }
    }

  }

  char plot_path[500];
  sprintf(plot_path,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/phys_dists/version%02d",version_number);
  for(int file=0; file<2; file++) {
    sprintf(name,"%s/wness_hists/dg0_vs_ddg0.png",plot_path);
    distributions[0][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/chi2_v_dcar.png",plot_path);
    distributions[1][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/rpcdca_1only.png",plot_path);
    distributions[2][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/rpcdca_3only.png",plot_path);
    distributions[3][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/rpcdca_1&3.png",plot_path);
    distributions[4][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/fvtx_dphi.png",plot_path);
    distributions[5][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/fvtx_dtheta_dr.png",plot_path);
    distributions[6][file]->SaveAs(name);
    sprintf(name,"%s/wness_hists/fvtx_cone.png",plot_path);
    distributions[7][file]->SaveAs(name);
  }

}

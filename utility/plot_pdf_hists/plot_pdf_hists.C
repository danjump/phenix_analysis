#include "../../get_phys_dists/phys_hist_definitions.h"
#include "../../get_phys_dists/version_number.h"
#include <TCanvas.h>

void plot_pdf_hists(char * sig_histinfilename, char * bkg_histinfilename) {
  
  TFile * histinfile[2];
  histinfile[0] = new TFile(sig_histinfilename);
  histinfile[1] = new TFile(bkg_histinfilename);
  
  //[arm][charge][goodfvtx]
  TCanvas * distributions[n_dists+2];

  char name[500];

  for(int i=0; i<n_dists+2; i++) {
    sprintf(name,"c_dist%d",i);
    distributions[i] = new TCanvas(name,name,1000,800);
    distributions[i]->Divide(3,2);
  }

  for(int file = 0; file<2; file++) {
    define_hists_phys(histinfile[file]);

    for(int arm=0; arm<2; arm++) {
      for(int charge=0; charge<2; charge++) {
        for(int goodfvtx=0; goodfvtx<2; goodfvtx++) {
          for(int rpcdca_condition=0; rpcdca_condition<3; rpcdca_condition++) {
            h2_dg0_ddg0[arm][charge][rpcdca_condition][goodfvtx]->Draw("COLZ");

            distributions[arm][charge][goodfvtx][1]->cd(rpcdca_condition+1+3*file);
            h2_chi2_dcar[arm][charge][rpcdca_condition][goodfvtx]->Draw("COLZ");
            for (int dists = 0; dists < n_dists; dists++) {
              distributions[2+dists]->cd(1+4*file+2*arm+charge);
              h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][dists]->SetLineColor(1+3*goodfvtx+rpcdca_condition);

              if(goodfvtx==0 && rpcdca_condition==0) {
              h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][dists]->Draw();
            }
          }
        }
      }
    }

  }

  char plot_path[500];
  sprintf(plot_path,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/phys_dists/version%02d",version_number);
  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      for(int goodfvtx=0; goodfvtx<2; goodfvtx++) {
        sprintf(name,"%s/arm%d_charge%d/dg0_vs_ddg0_fvtxstate%d.png",plot_path,arm,charge,goodfvtx);
        distributions[arm][charge][goodfvtx][0]->SaveAs(name);

        sprintf(name,"%s/arm%d_charge%d/chi2_vs_dcar_fvtxstate%d.png",plot_path,arm,charge,goodfvtx);
        distributions[arm][charge][goodfvtx][1]->SaveAs(name);

        for(int i=0; i<n_dists; i++) {
          sprintf(name,"%s/arm%d_charge%d/%s_fvtxstate%d.png",plot_path,arm,charge,distchar[i],goodfvtx);
          distributions[arm][charge][goodfvtx][i+2]->SaveAs(name);
        }
      }
    }
  }

}

#include <TH2F.h>
#include <../get_phys_dists/phys_hist_definitions.h>


TH2F * h2_dw23_vs_eta[11][2][2][10];
TH2F * h2_dw23_vs_wness[11][2][2];
TH1F * h_wness[11][2][2];
TH1F * h_dw23[11][2][2][10];
TH1F * h_wness_sections[11][2][2][10];
TH1F * h_eta[11][2][2][10];
TH1F * h_sbg_fit_eta[11][2][2][4][2];
TH1F * h_sbg_fit_dw23[11][2][2][4][2];
TH1F * h_target_wness[11][2][2];

int nhistbins_1dwness = 16;
int nhistbins_2dwness = 100;
int nhistbins_2ddw23 = 40;
int nhistbins_dw23eta_fit = 40;

int nhistbins_sbg_eta = 30;
int nhistbins_sbg_dw23 = 30;

void define_dw23_vs_eta_hists() {
  
  char histname[500];
  
  for(int source=0; source < 11; source++) {
    for(int arm=0; arm<n_arms; arm++) {
      for(int charge=0; charge<2; charge++) {
        
        sprintf(histname,"h_target_wness_source%d_arm%d_charge%d",source,arm,charge);
        h_target_wness[source][arm][charge] = new TH1F(histname,histname,300,0,1);

        for(int rpcdca_condition=0; rpcdca_condition<4; rpcdca_condition++) {
          for(int goodfvtx=0; goodfvtx<2; goodfvtx++) {

            if(source==0) {
              sprintf(histname,"h_sbg_fit_eta_source%d_arm%d_charge%d_rpcdca%d_fvtx%d",source,arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_eta,1.1,2.6);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();

              sprintf(histname,"h_sbg_fit_dw23_source%d_arm%d_charge%d_rpcdca%d_fvtx%d",source,arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_dw23,-.1,.1);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();
            }else if(source==1) {
              sprintf(histname,"h_sbg_fit_eta_source%d_arm%d_charge%d_rpcdca%d_fvtx%d",source,arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_eta,1.1,2.6);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();

              sprintf(histname,"h_sbg_fit_dw23_source%d_arm%d_charge%d_rpcdca%d_fvtx%d",source,arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_dw23,-.1,.1);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();
            } else {
              sprintf(histname,"h_sbg_fit_eta_%s_arm%d_charge%d_rpcdca%d_fvtx%d",mu_bg_labels[source-2],arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_eta,1.1,2.6);
              h_sbg_fit_eta[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();

              sprintf(histname,"h_sbg_fit_dw23_%s_arm%d_charge%d_rpcdca%d_fvtx%d",mu_bg_labels[source-2],arm,charge,rpcdca_condition,goodfvtx);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx] = new TH1F(histname,histname,nhistbins_sbg_dw23,-.1,.1);
              h_sbg_fit_dw23[source][arm][charge][rpcdca_condition][goodfvtx]->Sumw2();
            }
          }
        }

        sprintf(histname,"h2_dw23_vs_wness_source%d_arm%d_charge%d",source,arm,charge);
        h2_dw23_vs_wness[source][arm][charge] = new TH2F(histname,histname,nhistbins_2dwness,.1,.9,nhistbins_2ddw23,distmin[14],distmax[14]);
        
        sprintf(histname,"h_wness_source%d_arm%d_charge%d",source,arm,charge);
        h_wness[source][arm][charge] = new TH1F(histname,histname,nhistbins_1dwness,.1,.9);
        
        for(int wness_section=0; wness_section<10; wness_section++) {
          
          sprintf(histname,"h2_dw23_vs_eta_source%d_arm%d_charge%d_wnesscut%d",source,arm,charge,wness_section);
          h2_dw23_vs_eta[source][arm][charge][wness_section] = new TH2F(histname,histname,nhistbins_2ddw23,distmin[5],distmax[5],nhistbins_2ddw23,distmin[14],distmax[14]);
          
          sprintf(histname,"h_dw23_source%d_arm%d_charge%d_wnesscut%d",source,arm,charge,wness_section);
          h_dw23[source][arm][charge][wness_section] = new TH1F(histname,histname,nhistbins_2ddw23,distmin[14],distmax[14]);
          
          sprintf(histname,"h_wness_sections_source%d_arm%d_charge%d_wnesscut%d",source,arm,charge,wness_section);
          h_wness_sections[source][arm][charge][wness_section] = new TH1F(histname,histname,nhistbins_2dwness,0,1);
          
          sprintf(histname,"h_eta_source%d_arm%d_charge%d_wnesscut%d",source,arm,charge,wness_section);
          h_eta[source][arm][charge][wness_section] = new TH1F(histname,histname,nhistbins_2ddw23,distmin[5],distmax[5]);
        }
      }
    }
  }

}

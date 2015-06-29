#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
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
#include <vector>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooFormulaVar.h>
#include "wness_preselection.h"
#include "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/rpc_cluster_get_phys_dists/phys_hist_definitions.h"

void wness_preselection(const std::string inputtreefilename, const std::string bkg_treeinfilename,
    const std::string sig_treeinfilename, const std::string treeoutfilename,
    const std::string bkg_histoutfilename, const std::string sig_histoutfilename, 
    const std::string bkg_wnesshistoutfilename, const std::string sig_wnesshistoutfilename) { 

  bool do_rpc_cluster_cut = true;
  int rpc_cluster_cut = 3;

  TFile *treeinfile[2];
  treeinfile[0] = new TFile(bkg_treeinfilename.c_str());
  treeinfile[1] = new TFile(sig_treeinfilename.c_str());
  std::cout << "Using data (bkg) input file:\n" << bkg_treeinfilename << std::endl << std::endl;
  std::cout << "Using simulation (sig) input file:\n\n" << sig_treeinfilename << std::endl << std::endl;
  
  TFile *histoutfile[2];
  histoutfile[0] = new TFile(bkg_histoutfilename.c_str(),"RECREATE");
  histoutfile[1] = new TFile(sig_histoutfilename.c_str(),"RECREATE");
  std::cout << "Data (bkg) hists write file:\n" << bkg_histoutfilename << std::endl << std::endl;
  std::cout << "Simulation (sig) hists write file:\n\n" << sig_histoutfilename << std::endl << std::endl;
  
  TFile *wnesshistoutfile[2];
  wnesshistoutfile[0] = new TFile(bkg_wnesshistoutfilename.c_str(),"RECREATE");
  wnesshistoutfile[1] = new TFile(sig_wnesshistoutfilename.c_str(),"RECREATE");
  std::cout << "Data (bkg) wness hists write file:\n" << bkg_wnesshistoutfilename << std::endl << std::endl;
  std::cout << "Simulation (sig) wness hists write file:\n\n" << sig_wnesshistoutfilename << std::endl << std::endl;
 
  TFile *input_file = new TFile(inputtreefilename.c_str());
  
  Int_t
  Run_Number,  Evt_Number,  triggerbit,
  clockcross,  fvtx_cone, fvtx_tracklcone,
  rpc_awayclusters3;
  
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

  int entries,percent_done,percent_incriment,percent_done_previous;
  time_t rawtime;

  
  std::cout << "Setting likelihood distributions...";
  
  RooRealVar v_dg0("v_dg0","v_dg0",distmin[6],distmax[6]);
  RooRealVar v_ddg0("v_ddg0","v_ddg0",distmin[7],distmax[7]);
  RooRealVar v_chi2("v_chi2","v_chi2",distmin[9],distmax[9]);
  RooRealVar v_dcar("v_dcar","v_dcar",distmin[11],distmax[11]);
  RooRealVar v_rpc1dca("v_rpc1dca","v_rpc1dca",distmin[15],distmax[15]);
  RooRealVar v_rpc3dca("v_rpc3dca","v_rpc3dca",distmin[17],distmax[17]);
  RooRealVar v_fvtx_dphi("v_fvtx_dphi","v_fvtx_dphi",distmin[19],distmax[19]);
  RooRealVar v_fvtx_dr_dtheta("v_fvtx_dr_dtheta","v_fvtx_dr_dtheta",distmin[22],distmax[22]);
  RooRealVar v_fvtx_cone("v_fvtx_cone","v_fvtx_cone",distmin[23],distmax[23]);
  
  //[sig/bkg][arm][charge][rpcdca_condition][fvtx_Status]
  RooDataHist
  *h_dg0_vs_ddg0[2][n_arms][2][n_rpcdca_conditions][2],
  *h_chi2_vs_dcar[2][n_arms][2][n_rpcdca_conditions][2],
  *h_rpc1dca[2][n_arms][2][n_rpcdca_conditions][2],
  *h_rpc3dca[2][n_arms][2][n_rpcdca_conditions][2],
  *h_fvtx_dphi[2][n_arms][2][n_rpcdca_conditions],
  *h_fvtx_dr_dtheta[2][n_arms][2][n_rpcdca_conditions],
  *h_fvtx_cone[2][n_arms][2][n_rpcdca_conditions];
  
  RooHistPdf
  *pdf_dg0_vs_ddg0[2][n_arms][2][n_rpcdca_conditions][2],
  *pdf_chi2_vs_dcar[2][n_arms][2][n_rpcdca_conditions][2],
  *pdf_rpc1dca[2][n_arms][2][n_rpcdca_conditions][2],
  *pdf_rpc3dca[2][n_arms][2][n_rpcdca_conditions][2],
  *pdf_fvtx_dphi[2][n_arms][2][n_rpcdca_conditions],
  *pdf_fvtx_dr_dtheta[2][n_arms][2][n_rpcdca_conditions],
  *pdf_fvtx_cone[2][n_arms][2][n_rpcdca_conditions];
  
  RooProdPdf *likelihood[2][n_arms][2][n_rpcdca_conditions][2];


  for(int bgsig=0; bgsig<2; bgsig++) {
    printf("test1\n");
    define_basic_tree(treeinfile[bgsig]);
    define_hists_phys();
    define_hists_wness();

    printf("test2\n");
    basic_cuts_tree->SetBranchAddress("Run_Number",&Run_Number);
    basic_cuts_tree->SetBranchAddress("Evt_Number",&Evt_Number);
    basic_cuts_tree->SetBranchAddress("triggerbit",&triggerbit);
    basic_cuts_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
    basic_cuts_tree->SetBranchAddress("clockcross",&clockcross);
    basic_cuts_tree->SetBranchAddress("Wness",&Wness);
    basic_cuts_tree->SetBranchAddress("charge",&charge);
    basic_cuts_tree->SetBranchAddress("pT",&pT);
    basic_cuts_tree->SetBranchAddress("px",&px);
    basic_cuts_tree->SetBranchAddress("py",&py);
    basic_cuts_tree->SetBranchAddress("pz",&pz);
    basic_cuts_tree->SetBranchAddress("phi",&phi);
    basic_cuts_tree->SetBranchAddress("eta",&eta);
    basic_cuts_tree->SetBranchAddress("DG0",&DG0);
    basic_cuts_tree->SetBranchAddress("DDG0",&DDG0);
    basic_cuts_tree->SetBranchAddress("DG4",&DG4);
    basic_cuts_tree->SetBranchAddress("chi2",&chi2);
    basic_cuts_tree->SetBranchAddress("DCA_z",&DCA_z);
    basic_cuts_tree->SetBranchAddress("DCA_r",&DCA_r);
    basic_cuts_tree->SetBranchAddress("dphi12",&dphi12);
    basic_cuts_tree->SetBranchAddress("dphi23",&dphi23);
    basic_cuts_tree->SetBranchAddress("dw23",&dw23);
    basic_cuts_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
    basic_cuts_tree->SetBranchAddress("Rpc1time",&Rpc1time);
    basic_cuts_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
    basic_cuts_tree->SetBranchAddress("Rpc1x",&Rpc1x);
    basic_cuts_tree->SetBranchAddress("Rpc1y",&Rpc1y);
    basic_cuts_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
    basic_cuts_tree->SetBranchAddress("Rpc3time",&Rpc3time);
    basic_cuts_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
    basic_cuts_tree->SetBranchAddress("Rpc3x",&Rpc3x);
    basic_cuts_tree->SetBranchAddress("Rpc3y",&Rpc3y);
    basic_cuts_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
    basic_cuts_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
    basic_cuts_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
    basic_cuts_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
    basic_cuts_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
    basic_cuts_tree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);
    basic_cuts_tree->SetBranchAddress("rpc_awayclusters3",&rpc_awayclusters3);

    printf("test3\n");

    entries = basic_cuts_tree->GetEntries();

    percent_done=0;
    percent_incriment=20;

    cout << "\nNumber of events:  " << entries << endl;
    cout << "Starting Main Event Loop...\n\n";
    time(&rawtime);
    printf("Start time:  %s",ctime(&rawtime));

    for(int i=0; i<entries; i++) {//start event loop for getting histograms
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
      basic_cuts_tree->GetEntry(i);
      
      if(do_rpc_cluster_cut) {
        if( rpc_awayclusters3 > rpc_cluster_cut) continue;
      }
      
      // Set arm (0 = South, 1 = North)
      int arm = 0;
      if		(pz < 0) arm = 0;
      else if (pz > 0) arm = 1;


      int charge_index = 0;
      if		(charge < 0) charge_index = 0;
      else if (charge > 0) charge_index = 1;

      int rpcdca_condition = -1;
      if((Rpc1dca < 100) && (Rpc3dca < 100)) {
        rpcdca_condition = 2;
      } else if(Rpc1dca < 100) {
        rpcdca_condition = 0;
      } else if(Rpc3dca < 100) {
        rpcdca_condition = 1;
      } else {
        continue;
      }

      int goodfvtx;
      fvtx_dphi = fabs(fvtx_dphi);
      if(0 < fvtx_dphi && fvtx_dphi < 1.5 &&
          0 < fvtx_dtheta && fvtx_dtheta < 1.5  &&
          0 < fvtx_dr && fvtx_dr < 100) {
        bool fvtx_outofbounds = false;
        if(distmin[19] > fvtx_dphi || fvtx_dphi > distmax[19]) {
          printf("fvtx_dphi: %f  ",fvtx_dphi);
          fvtx_outofbounds = true;
        }
        if(distmin[20] > fvtx_dr || fvtx_dr > distmax[20]) {
          printf("fvtx_dr: %f  ",fvtx_dr);
          fvtx_outofbounds = true;
        }
        if(distmin[21] > fvtx_dtheta || fvtx_dtheta > distmax[21]) {
          printf("fvtx_dtheta: %f  ",fvtx_dtheta);
          fvtx_outofbounds = true;
        }
        if(distmin[22] > fvtx_dr_dtheta || fvtx_dr_dtheta > distmax[22]) {
          printf("fvtx_dr_dtheta: %f",fvtx_dr_dtheta);
          fvtx_outofbounds = true;
        }

        if(fvtx_outofbounds) {
          printf("\n");
          // continue;
        }

        goodfvtx = 1;

      } else {
        goodfvtx = 0;
      }

      float phys_dists[n_dists];
      phys_dists[0] = pT;
      phys_dists[1] = px;
      phys_dists[2] = py;
      phys_dists[3] = pz;
      phys_dists[4] = phi;
      phys_dists[3] = eta;
      phys_dists[6] = DG0;
      phys_dists[7] = DDG0;
      phys_dists[8] = DG4;
      phys_dists[9] = chi2;
      phys_dists[10] = DCA_z;
      phys_dists[11] = DCA_r;
      phys_dists[12] = dphi12;
      phys_dists[13] = dphi23;
      phys_dists[14] = dw23;
      phys_dists[15] = Rpc1dca;
      phys_dists[16] = Rpc1time;
      phys_dists[17] = Rpc3dca;
      phys_dists[18] = Rpc3time;
      phys_dists[19] = fabs(fvtx_dphi);
      phys_dists[20] = fvtx_dr;
      phys_dists[21] = fvtx_dtheta;
      phys_dists[22] = fvtx_dr_dtheta;
      phys_dists[23] = fvtx_cone;
      phys_dists[24] = Rpc1x;
      phys_dists[25] = Rpc1y;
      phys_dists[26] = Rpc3x;
      phys_dists[27] = Rpc3y;
      phys_dists[28] = fvtx_tracklcone;
      phys_dists[29] = Rpc1timewindow;
      phys_dists[30] = Rpc3timewindow;
      phys_dists[31] = rpc_awayclusters3;

      h2_dg0_ddg0[arm][charge_index][rpcdca_condition][goodfvtx]->Fill(DDG0,DG0);
      h2_chi2_dcar[arm][charge_index][rpcdca_condition][goodfvtx]->Fill(DCA_r,chi2);
      for (int dists = 0; dists < n_dists; dists++) {
        h_phys_dists[arm][charge_index][rpcdca_condition][goodfvtx][dists]->Fill(phys_dists[dists]);
      }

      h2_combined_dg0_ddg0[arm][charge_index]->Fill(DDG0,DG0);
      h2_combined_chi2_dcar[arm][charge_index]->Fill(DCA_r,chi2);

      if(rpcdca_condition != 1) {
        h_rpc_dists[arm][charge_index][0]->Fill(Rpc1dca);
      } if(rpcdca_condition != 0) {
        h_rpc_dists[arm][charge_index][1]->Fill(Rpc3dca);
      }
      if(goodfvtx) {
        h_fvtx_dists[arm][charge_index][0]->Fill(fvtx_dphi);
        h_fvtx_dists[arm][charge_index][1]->Fill(fvtx_dr_dtheta);
        h_fvtx_dists[arm][charge_index][2]->Fill(fvtx_cone);
      }

    } //end event loop for getting histograms
    
    histoutfile[bgsig]->cd();
    save_hists_phys();

    wnesshistoutfile[bgsig]->cd();
    save_hists_wness();
    
    for(int arm=0; arm<n_arms; arm++) {
      for(int charge=0; charge<2; charge++) {
        for(int rpcdca_condition=0; rpcdca_condition<n_rpcdca_conditions; rpcdca_condition++) {
          for(int goodfvtx=0; goodfvtx<2; goodfvtx++) {
            char name[200];

            sprintf(name,"h_dg0_vs_ddg0_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            h_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,RooArgList(v_ddg0,v_dg0),h2_combined_dg0_ddg0[arm][charge]);
//            h_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,RooArgList(v_ddg0,v_dg0),h2_dg0_ddg0[arm][charge][rpcdca_condition][goodfvtx]);
            sprintf(name,"pdf_dg0_vs_ddg0_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooHistPdf(name,name,RooArgList(v_ddg0,v_dg0),*h_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx]);

            sprintf(name,"h_chi2_vs_dcar_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            h_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,RooArgList(v_dcar,v_chi2),h2_combined_chi2_dcar[arm][charge]);
//            h_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,RooArgList(v_dcar,v_chi2),h2_chi2_dcar[arm][charge][rpcdca_condition][goodfvtx]);
            sprintf(name,"pdf_chi2_vs_dcarbgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooHistPdf(name,name,RooArgList(v_dcar,v_chi2),*h_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx]);

            sprintf(name,"h_rpc1dca_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            h_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,v_rpc1dca,h_rpc_dists[arm][charge][0]);
//            h_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,v_rpc1dca,h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][12]);
            sprintf(name,"pdf_rpc1dca_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            pdf_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooHistPdf(name,name,v_rpc1dca,*h_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx]);

            sprintf(name,"h_rpc3dca_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            h_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,v_rpc3dca,h_rpc_dists[arm][charge][1]);
//            h_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooDataHist(name,name,v_rpc3dca,h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][13]);
            sprintf(name,"pdf_rpc3dca_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            pdf_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooHistPdf(name,name,v_rpc3dca,*h_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx]);

            sprintf(name,"likelihood_bgsig%d_arm%d_charge%d_rpcdca%d_goodfvtx%d",bgsig,arm,charge,rpcdca_condition,goodfvtx);
            if(goodfvtx==0) {
              if(rpcdca_condition == 0) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                    RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx]),0);

              } else if(rpcdca_condition == 1) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                    RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx]),0);

              } else if(rpcdca_condition == 2) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                   RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                              *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                              *pdf_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                              *pdf_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx]),0);
              }
            } else {

              sprintf(name,"h_fvtx_dphi_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
              h_fvtx_dphi[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_dphi,h_fvtx_dists[arm][charge][0]);
//              h_fvtx_dphi[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_dphi,h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][14]);
              sprintf(name,"pdf_fvtx_dphi_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
              pdf_fvtx_dphi[bgsig][arm][charge][rpcdca_condition] = new RooHistPdf(name,name,v_fvtx_dphi,*h_fvtx_dphi[bgsig][arm][charge][rpcdca_condition]);

              sprintf(name,"h_fvtx_dr_dtheta_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
              h_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_dr_dtheta,h_fvtx_dists[arm][charge][1]);
//              h_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_dr_dtheta,h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][17]);
              sprintf(name,"pdf_fvtx_dr_dtheta_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
              pdf_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition] = new RooHistPdf(name,name,v_fvtx_dr_dtheta,*h_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition]);
              
             sprintf(name,"h_fvtx_cone_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
             h_fvtx_cone[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_cone,h_fvtx_dists[arm][charge][2]);
//               h_fvtx_cone[bgsig][arm][charge][rpcdca_condition] = new RooDataHist(name,name,v_fvtx_cone,h_phys_dists[arm][charge][rpcdca_condition][goodfvtx][18]);
             sprintf(name,"pdf_fvtx_cone_bgsig%d_arm%d_charge%d_rpcdca%d",bgsig,arm,charge,rpcdca_condition);
             pdf_fvtx_cone[bgsig][arm][charge][rpcdca_condition] = new RooHistPdf(name,name,v_fvtx_cone,*h_fvtx_cone[bgsig][arm][charge][rpcdca_condition]);

              if(rpcdca_condition == 0) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                    RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_fvtx_dphi[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_cone[bgsig][arm][charge][rpcdca_condition]),0);

              } else if(rpcdca_condition == 1) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                    RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_fvtx_dphi[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_cone[bgsig][arm][charge][rpcdca_condition]),0);

              } else if(rpcdca_condition == 2) {
                likelihood[bgsig][arm][charge][rpcdca_condition][goodfvtx] = new RooProdPdf(name,name,
                    RooArgList(*pdf_dg0_vs_ddg0[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_chi2_vs_dcar[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc1dca[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_rpc3dca[bgsig][arm][charge][rpcdca_condition][goodfvtx],
                               *pdf_fvtx_dphi[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_dr_dtheta[bgsig][arm][charge][rpcdca_condition],
                               *pdf_fvtx_cone[bgsig][arm][charge][rpcdca_condition]),0);

              }
            }
          }
        }
      }
    }
    printf("end first loop\n");
  }
  
  std::cout << "done!\n\n";
  
  std::cout << "Preparing for loop through events to calculate wness...\n";



  TFile *treeoutfile = new TFile(treeoutfilename.c_str(),"RECREATE");
  define_wness_tree();

  percent_done=0;
  percent_incriment=20;

  define_basic_tree(input_file);

  basic_cuts_tree->SetBranchAddress("Run_Number",&Run_Number);
  basic_cuts_tree->SetBranchAddress("Evt_Number",&Evt_Number);
  basic_cuts_tree->SetBranchAddress("triggerbit",&triggerbit);
  basic_cuts_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  basic_cuts_tree->SetBranchAddress("clockcross",&clockcross);
  basic_cuts_tree->SetBranchAddress("Wness",&Wness);
  basic_cuts_tree->SetBranchAddress("charge",&charge);
  basic_cuts_tree->SetBranchAddress("pT",&pT);
  basic_cuts_tree->SetBranchAddress("px",&px);
  basic_cuts_tree->SetBranchAddress("py",&py);
  basic_cuts_tree->SetBranchAddress("pz",&pz);
  basic_cuts_tree->SetBranchAddress("phi",&phi);
  basic_cuts_tree->SetBranchAddress("eta",&eta);
  basic_cuts_tree->SetBranchAddress("DG0",&DG0);
  basic_cuts_tree->SetBranchAddress("DDG0",&DDG0);
  basic_cuts_tree->SetBranchAddress("DG4",&DG4);
  basic_cuts_tree->SetBranchAddress("chi2",&chi2);
  basic_cuts_tree->SetBranchAddress("DCA_z",&DCA_z);
  basic_cuts_tree->SetBranchAddress("DCA_r",&DCA_r);
  basic_cuts_tree->SetBranchAddress("dphi12",&dphi12);
  basic_cuts_tree->SetBranchAddress("dphi23",&dphi23);
  basic_cuts_tree->SetBranchAddress("dw23",&dw23);
  basic_cuts_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  basic_cuts_tree->SetBranchAddress("Rpc1time",&Rpc1time);
  basic_cuts_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
  basic_cuts_tree->SetBranchAddress("Rpc1x",&Rpc1x);
  basic_cuts_tree->SetBranchAddress("Rpc1y",&Rpc1y);
  basic_cuts_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  basic_cuts_tree->SetBranchAddress("Rpc3time",&Rpc3time);
  basic_cuts_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
  basic_cuts_tree->SetBranchAddress("Rpc3x",&Rpc3x);
  basic_cuts_tree->SetBranchAddress("Rpc3y",&Rpc3y);
  basic_cuts_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  basic_cuts_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  basic_cuts_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  basic_cuts_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
  basic_cuts_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);
  basic_cuts_tree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);
  basic_cuts_tree->SetBranchAddress("rpc_awayclusters3",&rpc_awayclusters3);

  entries = basic_cuts_tree->GetEntries();
  
  cout << "\nNumber of events:  " << entries << endl;
  cout << "Starting Main Event Loop...\n\n";
  time(&rawtime);
  printf("Start time:  %s",ctime(&rawtime));

  
  for(int i=0; i<entries; i++) {//Main loop
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
    basic_cuts_tree->GetEntry(i);
    // Set arm (0 = South, 1 = North)
    int arm = 0;
    if		(pz < 0) arm = 0;
    else if (pz > 0) arm = 1;
    
    
    int charge_index = 0;
    if		(charge < 0) charge_index = 0;
    else if (charge > 0) charge_index = 1;

    if(do_rpc_cluster_cut) {
      if(rpc_awayclusters3 > rpc_cluster_cut) continue;
    }
    
    int rpcdca_condition = -1;
    if((Rpc1dca < 100) && (Rpc3dca < 100)) {
      rpcdca_condition = 2;
    } else if(Rpc1dca < 100) {
      rpcdca_condition = 0;
    } else if(Rpc3dca < 100) {
      rpcdca_condition = 1;
    } else {
      continue;
    }

    int goodfvtx;
    fvtx_dphi = fabs(fvtx_dphi);
    if(0 < fvtx_dphi && fvtx_dphi < 1.5 &&
         0 < fvtx_dtheta && fvtx_dtheta < 1.5  &&
         0 < fvtx_dr && fvtx_dr < 100) {
      bool fvtx_outofbounds = false;
      if(distmin[19] > fvtx_dphi || fvtx_dphi > distmax[19]) {
        printf("fvtx_dphi: %f  ",fvtx_dphi);
        fvtx_outofbounds = true;
      }
      if(distmin[20] > fvtx_dr || fvtx_dr > distmax[20]) {
        printf("fvtx_dr: %f  ",fvtx_dr);
        fvtx_outofbounds = true;
      }
      if(distmin[21] > fvtx_dtheta || fvtx_dtheta > distmax[21]) {
        printf("fvtx_dtheta: %f  ",fvtx_dtheta);
        fvtx_outofbounds = true;
      }
      if(distmin[22] > fvtx_dr_dtheta || fvtx_dr_dtheta > distmax[22]) {
        printf("fvtx_dr_dtheta: %f",fvtx_dr_dtheta);
        fvtx_outofbounds = true;
      }

      if(fvtx_outofbounds) {
        printf("\n");
       // continue;
      }

      goodfvtx = 1;
      
    } else {
      goodfvtx = 0;
    }
    
    v_dg0     = DG0;
    v_ddg0    = DDG0;
    v_chi2    = chi2;
    v_dcar    = DCA_r;
    v_rpc1dca = Rpc1dca;
    v_rpc3dca = Rpc3dca;

    //printf("dg0:%f ddg0:%f chi2:%f dcar%f rpc1dca%f rpc3dca%f ",DG0,DDG0,chi2,DCA_r,Rpc1dca,Rpc3dca);
    //printf("dphi%f dr%f dtheta%f drdtheta%f ",fvtx_dphi,fvtx_dr,fvtx_dtheta,fvtx_dr_dtheta);
    //printf("goodfvtx%d rpc_dca%d\n",goodfvtx,rpcdca_condition);
    RooArgSet *var_set = new RooArgSet();
    if(goodfvtx) {
      v_fvtx_dphi = fvtx_dphi;
      v_fvtx_dr_dtheta = fvtx_dr_dtheta;
      v_fvtx_cone = fvtx_cone;
      if(rpcdca_condition==0) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc1dca,v_fvtx_dphi,v_fvtx_dr_dtheta,v_fvtx_cone);
      } else if(rpcdca_condition==1) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc3dca,v_fvtx_dphi,v_fvtx_dr_dtheta,v_fvtx_cone);
      } else if(rpcdca_condition==2) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc1dca,v_rpc3dca,v_fvtx_dphi,v_fvtx_dr_dtheta,v_fvtx_cone);
      }

    } else {
      if(rpcdca_condition==0) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc1dca);
      } else if(rpcdca_condition==1) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc3dca);
      } else if(rpcdca_condition==2) {
        var_set = new RooArgSet(v_dg0,v_ddg0,v_chi2,v_dcar,v_rpc1dca,v_rpc3dca);
      }

    }
    

    double bkg_like = likelihood[0][arm][charge_index][rpcdca_condition][goodfvtx]->getVal(var_set);
    double sig_like = likelihood[1][arm][charge_index][rpcdca_condition][goodfvtx]->getVal(var_set);

    Wness  = sig_like / ( bkg_like + sig_like);

    wness_tree->SetBranchAddress("px",&px);
    wness_tree->SetBranchAddress("py",&py);
    wness_tree->SetBranchAddress("dphi12",&dphi12);
    wness_tree->SetBranchAddress("Rpc1time",&Rpc1time);
    wness_tree->SetBranchAddress("Rpc1timewindow",&Rpc1timewindow);
    wness_tree->SetBranchAddress("Rpc1x",&Rpc1x);
    wness_tree->SetBranchAddress("Rpc1y",&Rpc1y);
    wness_tree->SetBranchAddress("Rpc3time",&Rpc3time);
    wness_tree->SetBranchAddress("Rpc3timewindow",&Rpc3timewindow);
    wness_tree->SetBranchAddress("Rpc3x",&Rpc3x);
    wness_tree->SetBranchAddress("Rpc3y",&Rpc3y);
    wness_tree->SetBranchAddress("Run_Number", &Run_Number);
    wness_tree->SetBranchAddress("Evt_Number", &Evt_Number);
    wness_tree->SetBranchAddress("triggerbit", &triggerbit);
    wness_tree->SetBranchAddress("Evt_bbcZ", &Evt_bbcZ);
    wness_tree->SetBranchAddress("clockcross", &clockcross);
    wness_tree->SetBranchAddress("Wness", &Wness);
    wness_tree->SetBranchAddress("charge", &charge);
    wness_tree->SetBranchAddress("pT", &pT);
    wness_tree->SetBranchAddress("pz", &pz);
    wness_tree->SetBranchAddress("phi", &phi);
    wness_tree->SetBranchAddress("eta", &eta);
    wness_tree->SetBranchAddress("DG0", &DG0);
    wness_tree->SetBranchAddress("DDG0", &DDG0);
    wness_tree->SetBranchAddress("DG4", &DG4);
    wness_tree->SetBranchAddress("chi2", &chi2);
    wness_tree->SetBranchAddress("DCA_z", &DCA_z);
    wness_tree->SetBranchAddress("DCA_r", &DCA_r);
    wness_tree->SetBranchAddress("dphi23", &dphi23);
    wness_tree->SetBranchAddress("dw23", &dw23);
    wness_tree->SetBranchAddress("Rpc1dca", &Rpc1dca);
    wness_tree->SetBranchAddress("Rpc3dca", &Rpc3dca);
    wness_tree->SetBranchAddress("fvtx_dphi", &fvtx_dphi);
    wness_tree->SetBranchAddress("fvtx_dr", &fvtx_dr);
    wness_tree->SetBranchAddress("fvtx_dtheta", &fvtx_dtheta);
    wness_tree->SetBranchAddress("fvtx_dr_dtheta", &fvtx_dr_dtheta);
    wness_tree->SetBranchAddress("fvtx_cone", &fvtx_cone);
    wness_tree->SetBranchAddress("fvtx_tracklcone", &fvtx_tracklcone);
    wness_tree->SetBranchAddress("rpc_awayclusters3", &rpc_awayclusters3);
    
    wness_tree->Fill();
    
  } //end main event loop
    
  std::cout << "Saving new wness tree in:\n" << treeoutfilename << std::endl << std::endl;
  
  treeoutfile->cd();
  save_wness_tree();
}



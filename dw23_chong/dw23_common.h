#ifndef DW23_COMMON_H
#define DW23_COMMON_H

#include <TROOT.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>
#include <TH2F.h>
#include "universal_constants.h"

int nbins_wness=50;
int nbins_dw23=30;

//Utility function for opening a file and filling a 2d dw23_vs_wness histogram
void get_dw23_from_file(int fileindex, bool do_rpc_cluster, int cluster_cut, TH2F *h2_dw23_vs_wness[][2],TH1F *h_wness_slices[][2][10],TH1F *h_dw23_slices[][2][10]) {

  char * wnessfile[11] = {
    (char*)"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_total_wness_refrun367593.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_w_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_dy_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_light_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_onium_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_onlyz_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_openbottom_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_opencharm_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_whad_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_wtau_combined_wness.root",
    (char*)"/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_z_combined_wness.root"};

  char * rpccluster_wnessfile[11] = {
    (char*)Form("/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_data_total_wness_refrun367593.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_w_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_dy_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_light_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_onium_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_onlyz_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_openbottom_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_opencharm_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_whad_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_wtau_combined_wness.root",cluster_cut),
    (char*)Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_z_combined_wness.root",cluster_cut)};
  
  TFile *infile;
  infile = TFile::Open("/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/nocut/run13_data_total_wness_refrun367593.root" );
  /*if(do_rpc_cluster) {
    infile = TFile::Open( rpccluster_wnessfile[fileindex] );
    cout << "Opening file:\n" << rpccluster_wnessfile[fileindex] << endl;
  } else{
    infile = TFile::Open( wnessfile[fileindex] );
    cout << "Opening file:\n" << wnessfile[fileindex] << endl;
  }*/

  TTree *wness_tree = (TTree*)infile->Get("wness_tree");
  
  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      h2_dw23_vs_wness[a][c] = new TH2F(
          Form("h2_dw23_vs_wness_%s_a%d_c%d",processname[fileindex],a,c),
          Form("h2_dw23_vs_wness_%s_a%d_c%d",processname[fileindex],a,c),
          nbins_wness,0,1,nbins_dw23,-.3,.3);
      h2_dw23_vs_wness[a][c]->Sumw2();
      
      TCut armcut,chargecut,clustercut;
      if(a) armcut = "pz > 0"; else armcut= "pz < 0";
      if(c) chargecut = "charge > 0"; else chargecut = "charge < 0";
      clustercut = "";//if(do_rpc_cluster) clustercut = Form("rpc_awayclusters3 <=%d",cluster_cut); else clustercut = "";

      wness_tree->Project(h2_dw23_vs_wness[a][c]->GetName(),"dw23:Wness",armcut+chargecut+clustercut);

      for(int w=0; w<10; w++) {
        h_dw23_slices[a][c][w] = new TH1F(
            Form("h_dw23_slice%d_a%d_c%d",w,a,c),
            Form("h_dw23_slice%d_a%d_c%d",w,a,c),
            nbins_dw23,-.3,.3);
        h_dw23_slices[a][c][w]->Sumw2();
        
        h_wness_slices[a][c][w] = new TH1F(
            Form("h_wness_slice%d_a%d_c%d",w,a,c),
            Form("h_wness_slice%d_a%d_c%d",w,a,c),
            nbins_wness,0,1);
        h_wness_slices[a][c][w]->Sumw2();

        TCut wnesscut;
        if(w!=9)
          wnesscut = Form("Wness > %f && Wness < %f", .1*w,.1*(w+1));
        else
          wnesscut = Form("Wness > %f && Wness < %f", .92,1.0);

        wness_tree->Project(h_dw23_slices[a][c][w]->GetName(),"dw23",armcut+chargecut+wnesscut+clustercut);
        wness_tree->Project(h_wness_slices[a][c][w]->GetName(),"Wness",armcut+chargecut+wnesscut+clustercut);
      }
    }
  }

}


#endif


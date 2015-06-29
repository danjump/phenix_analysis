#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

void run_totaltrig_sbg_fit_rpccluster(int cluster_cut=0) {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/sbg_fit";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libsbg_fit.so",code_dir);
  gSystem->Load(lib_file);
  sbg_fit(
      Form("/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_data_total_wness_refrun367593.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_w_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_dy_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_light_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_onium_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_onlyz_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_openbottom_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_opencharm_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_whad_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_wtau_combined_wness.root",cluster_cut),
      Form("/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/cut_%d/run13_simulation_z_combined_wness.root",cluster_cut));
  cout << "done\n";
}

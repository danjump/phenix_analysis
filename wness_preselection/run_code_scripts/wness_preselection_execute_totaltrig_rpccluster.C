
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooFormulaVar.h>

//simulation reference run numbers: 367466 367593 368630
void wness_preselection_execute_totaltrig_rpccluster(int which_ref_run_index, int run_over_which) {
  //run_over_which is 0-2 for simulation ref run 0-2 and 12 or 13 for data
  
  char inputtreefilename[500], bkg_treeinfilename[500], sig_treeinfilename[500],
       treeoutfilename[500], bkg_histoutfilename[500], sig_histoutfilename[500], 
       bkg_wnesshistoutfilename[500], sig_wnesshistoutfilename[500];
  
  int ref_run[3];
  ref_run[0] = 367466;
  ref_run[1] = 367593;
  ref_run[2] = 368630;

  int which_year = 13; 
  
  bool changed=false;
  
  if(run_over_which>10) {
    sprintf(inputtreefilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_data_total_trees_combined.root");
    sprintf(treeoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/run13_data_total_wness_refrun%d.root",ref_run[which_ref_run_index]);
  } else {
    sprintf(inputtreefilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_simulation_w_combined_tree_combined.root");
    sprintf(treeoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/run13_simulation_w_combined_wness.root");
  }
  
  sprintf(bkg_treeinfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_data_total_trees_combined.root");
  sprintf(bkg_histoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_data_total_hists_combined.root");
  sprintf(bkg_wnesshistoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/run13_data_total_hists_combined.root");

  sprintf(sig_treeinfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_simulation_w_combined_tree_combined.root");
  sprintf(sig_histoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_simulation_w_hists_combined.root");
  sprintf(sig_wnesshistoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_wness_tree/run13_simulation_w_hists_combined.root");

  
  
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/wness_preselection";
  char lib_file[300];
	sprintf(lib_file,"%s/install/lib/libwness_preselection.so",code_dir);
	gSystem->Load(lib_file);
	wness_preselection(inputtreefilename,bkg_treeinfilename,sig_treeinfilename,treeoutfilename,bkg_histoutfilename,sig_histoutfilename,bkg_wnesshistoutfilename,sig_wnesshistoutfilename);
}

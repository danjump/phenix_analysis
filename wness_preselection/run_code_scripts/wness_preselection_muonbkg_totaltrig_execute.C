
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooFormulaVar.h>

//simulation reference run numbers: 367466 367593 368630
void wness_preselection_muonbkg_totaltrig_execute_rpccluster(int which_ref_run_index, int run_over_which) {
  
  char inputtreefilename[500], bkg_treeinfilename[500], sig_treeinfilename[500],
       treeoutfilename[500], bkg_histoutfilename[500], sig_histoutfilename[500], 
       bkg_wnesshistoutfilename[500], sig_wnesshistoutfilename[500];
  
  int ref_run[3];
  ref_run[0] = 367466;
  ref_run[1] = 367593;
  ref_run[2] = 368630;

  int which_year = 13;

  char * label[9] = {"dy",
    "light","onium","onlyz",
    "openbottom","opencharm",
    "whad","wtau","z"};
  
  bool changed=false;
  
  sprintf(inputtreefilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_phys_dists/run%d_simulation_%s_combined_tree_combined.root",
      which_year,which_year,label[run_over_which]);
  sprintf(treeoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_wness_tree/run%d_simulation_%s_combined_wness.root",
      which_year,which_year,label[run_over_which]);


  sprintf(bkg_treeinfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_phys_dists/run%d_data_total_trees_combined.root",
      which_year,which_year);
  sprintf(bkg_histoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_clusterphys_dists/run%d_data_total_hists_combined.root",
      which_year,which_year);
  sprintf(bkg_wnesshistoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_clusterwness_tree/run%d_data_total_hists_combined.root",
      which_year,which_year);

  
  sprintf(sig_treeinfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_phys_dists/run%d_simulation_w_combined_tree_combined.root",
      which_year,which_year);
  sprintf(sig_histoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_phys_dists/run%d_simulation_w_hists_combined.root",
      which_year,which_year);
  sprintf(sig_wnesshistoutfilename,"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run%d/rpc_cluster_wness_tree/run%d_simulation_w_hists_combined.root",
      which_year,which_year);
  
  
  
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/wness_preselection";
  char lib_file[300];
	sprintf(lib_file,"%s/install/lib/libwness_preselection.so",code_dir);
	gSystem->Load(lib_file);
	wness_preselection(inputtreefilename,bkg_treeinfilename,sig_treeinfilename,treeoutfilename,bkg_histoutfilename,sig_histoutfilename,bkg_wnesshistoutfilename,sig_wnesshistoutfilename);
}

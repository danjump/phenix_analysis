#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

void run_totaltrig_eta_dw23_fit() {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/eta_dw23_fit";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libeta_dw23_fit.so",code_dir);
  gSystem->Load(lib_file);
  eta_dw23_fit(
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_dy_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_light_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onium_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onlyz_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_openbottom_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_opencharm_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_whad_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_wtau_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_z_total_wness_refrun367593.root");
  cout << "done\n";
}

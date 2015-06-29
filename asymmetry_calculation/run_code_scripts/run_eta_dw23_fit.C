#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

void run_eta_dw23_fit() {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/eta_dw23_fit";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libeta_dw23_fit.so",code_dir);
  gSystem->Load(lib_file);
  eta_dw23_fit(
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_dy_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_light_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onium_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onlyz_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_openbottom_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_opencharm_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_whad_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_wtau_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_z_wness_refrun367593.root");
  cout << "done\n";
}

#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

void run_totaltrig_dw23_iteration() {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_iteration";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libdw23_iteration.so",code_dir);
  gSystem->Load(lib_file);
  dw23_iteration(
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_total_wness_refrun367593.root",
      /*"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_dy_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_light_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onium_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_onlyz_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_openbottom_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_opencharm_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_whad_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_wtau_total_wness_refrun367593.root",
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_z_total_wness_refrun367593.root"*/
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_w_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_dy_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_light_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_onium_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_onlyz_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_openbottom_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_opencharm_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_whad_combined_wness.root",
            //"/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_sim_wtau_total_wness_refrun367593.root",
            "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_wtau_combined_wness.root",
      "/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_simulation_z_combined_wness.root");
  cout << "done\n";
}

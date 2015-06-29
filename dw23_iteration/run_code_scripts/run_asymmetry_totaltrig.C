#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>

void run_asymmetry_totaltrig() {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libasymmetry_calculation.so",code_dir);
  gSystem->Load(lib_file);
  asymmetry_calculation(
      "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_total_wness_refrun367593.root");
  cout << "done\n";
}

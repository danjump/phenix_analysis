
/*void rpc_cluster_get_phys_dists_execute(char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_367466_combined.root",
        char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367466_combined_hists.root",
        char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367466_combined_tree.root") {//*/
/*void rpc_cluster_get_phys_dists_execute(char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_367593_combined.root",
        char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367593_combined_hists.root",
        char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367593_combined_tree.root") {//*/
/*void rpc_cluster_get_phys_dists_execute(char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_368630_combined.root",
        char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_368630_combined_hists.root",
        char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_368630_combined_tree.root") {//*/
void run_rpc_cluster_get_phys_dists() {
int run_over_which=13;
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/rpc_cluster_get_phys_dists";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/librpc_cluster_get_phys_dists.so",code_dir);
  gSystem->Load(lib_file);

  if(run_over_which==12) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_data_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_data_hists_combined.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_data_trees_combined.root";
  } else if(run_over_which==13) {
    //char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files/run13_data_combined.root";
    char * infilename = "/direct/phenix+spin2/rseidl/taxi/Run13pp510Muon/3437/data/398149_0.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/temprun13_data_hists_combined.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/temprun13_data_trees_combined.root";
  } else if(run_over_which==3) {  
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_368630_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_368630_combined_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_368630_combined_tree.root";
  } else if(run_over_which==2) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_367593_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367593_combined_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367593_combined_tree.root";
  } else if(run_over_which==1) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulation_367466_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367466_combined_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulation_367466_combined_tree.root";
  } else if(run_over_which==0) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simulations_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulations_combined_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/rpc_cluster_phys_dists/run12_simulations_combined_tree.root";
  } else {
    char * infilename = "";
    char * histoutfilename = "";
    char * treeoutfilename = "";
  }
  rpc_cluster_get_phys_dists(infilename,histoutfilename,treeoutfilename);

}

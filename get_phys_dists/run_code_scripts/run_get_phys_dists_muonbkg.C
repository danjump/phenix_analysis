void run_get_phys_dists_muonbkg(int run_over_which) {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libget_phys_dists.so",code_dir);
  gSystem->Load(lib_file);

  if(run_over_which==0) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_dy_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_dy_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_dy_tree.root";
  } else if(run_over_which==1) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_light_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_light_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_light_tree.root";
  } else if(run_over_which==2) {  
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_onium_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_onium_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_onium_tree.root";
  } else if(run_over_which==3) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_onlyz_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_onlyz_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_onlyz_tree.root";
  } else if(run_over_which==4) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_openbottom_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_openbottom_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_openbottom_tree.root";
  } else if(run_over_which==5) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_opencharm_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_opencharm_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_opencharm_tree.root";
  } else if(run_over_which==6) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_whad_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_whad_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_whad_tree.root";
  } else if(run_over_which==7) {  
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_wtau_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_wtau_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_wtau_tree.root";
  } else if(run_over_which==8) {
    char * infilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files/run12_simu367593_z_combined.root";
    char * histoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_z_hists.root";
    char * treeoutfilename = "/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/phys_dists/run12_sim367593_z_tree.root";
  } else {
    char * infilename = "";
    char * histoutfilename = "";
    char * treeoutfilename = "";
  }
  get_phys_dists(infilename,histoutfilename,treeoutfilename);

}


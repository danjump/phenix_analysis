

void get_phys_dists_execute(char * infilename, char * histoutfilename, char * treeoutfilename) {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/libget_phys_dists.so",code_dir);
  gSystem->Load(lib_file);

  get_phys_dists(infilename,histoutfilename,treeoutfilename);

}

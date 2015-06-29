
//example manual execution:
//  root 'rpc_cluster_get_phys_dists_execute.C("/direct/phenix+spin/phnxsp01/rseidl/Run13pp510Muon/ana103/398146.root","/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13datarp/hists/hist_398146.root","/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13datarp/trees/tree_398146.root")'

void rpc_cluster_get_phys_dists_execute(char * infilename, char * histoutfilename, char * treeoutfilename) {
  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/rpc_cluster_get_phys_dists";
  char lib_file[300];
  sprintf(lib_file,"%s/install/lib/librpc_cluster_get_phys_dists.so",code_dir);
  gSystem->Load(lib_file);

  rpc_cluster_get_phys_dists(infilename,histoutfilename,treeoutfilename);

}

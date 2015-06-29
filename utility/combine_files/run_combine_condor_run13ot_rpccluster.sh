
#outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files"
outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists"
label="run13_dataot_trees"
ls_dir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13dataot/trees"
search_prefix=""

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix

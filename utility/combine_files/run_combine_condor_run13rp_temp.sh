
#outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files"
outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/temp_phys_dists"
mkdir -p $outfile_path
label="run13_datarp_trees"
ls_dir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13datarp/trees"
search_prefix=""

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix


#outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files"
outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists"
label="run13_data_trees"
ls_dir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13data/trees"
search_prefix=""

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix

which_sim="wtau"
outfile_path="/direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists"
label="run13_simulation_${which_sim}_combined_tree"
ls_dir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13_sim_${which_sim}/trees"
search_prefix="tree"

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix


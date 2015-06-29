
which_sim="w"
outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/summed_event_files"
label="run12_sim367593_${which_sim}"
ls_dir="/direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/${which_sim}_sum"
search_prefix="pdst"

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix

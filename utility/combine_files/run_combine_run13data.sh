
#outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files"
outfile_path="/direct/phenix+scratch/danielj/temp"
label="run13_data"
ls_dir="/direct/phenix+spin2/rseidl/Run13pp510Muon/2103/data"
search_prefix="3"

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix

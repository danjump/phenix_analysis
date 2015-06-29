#!/bin/bash

#outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/summed_event_files"
outfile_path="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists"
label="run13_data_ot_trees"
ls_dir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13dataot/trees"
search_prefix="tree"

hadd_combine_files.sh $outfile_path $label $ls_dir $search_prefix

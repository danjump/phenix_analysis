#!/bin/bash

macro=/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/run_code_scripts/process_pdst_execute.C
outfilepath=/direct/phenix+scratch/danielj/run13analysis_output/physdisthists

sim_path=/direct/phenix+spin2/rseidl/Wsims/run12sim/muonsims/simdata/pytune100_367593/old2/w_sum
sim_run=`echo $sim_path | awk -F/ '{print $9}' | cut -c11-`
for file in `ls $sim_path/pdst_w_*.root`
do

  segment=`basename $file .root| sed 's/pdst_w_//g'`
  echo "Running over sim_segment $segment..."
  
    root -l -b -q $macro\(\"$file\",\"$outfilepath/sim_physhists-$sim_run-$segment.root\"\)
  echo "Done with segment $segment!"
  echo ""


done

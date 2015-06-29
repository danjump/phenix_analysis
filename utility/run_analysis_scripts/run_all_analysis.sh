#!/bin/bash

which_year=13

code_root_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis"

#accept optional starting point argument. set to 0 if no argument given
starting_point=${1:-0}

if [ $starting_point -lt 1 ]; then   #sp=0
  echo "Starting on data condor"
elif [ $starting_point -lt 2 ]; then #sp=1
  echo "Starting on OT condor"
elif [ $starting_point -lt 3 ]; then #sp=2
  echo "Starting on merging trigger files"
elif [ $starting_point -lt 4 ]; then #sp=3
  echo "Starting on wness calculations"
fi

echo ""
echo pausing 5 sec...
sleep 5
echo ""

#first run get_phys_dists through condor on all data files
#you will want to make sure the filelist is up to date in the following directory:
# /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/run_code_scripts/condor/filelists
if [ $starting_point -lt 1 ]; then
  echo "Running normal data condor...`date`"
  ${code_root_dir}/get_phys_dists/run_code_scripts/condor/data_run13_condor.sh
  echo "Finished submitting normal data condor...`date`"
fi
if [ $starting_point -lt 2 ]; then
  echo "Running ot data condor...`date`"
  ${code_root_dir}/get_phys_dists/run_code_scripts/condor/data_run13_ot_condor.sh
  echo "Finished submitting ot data condor...`date`"
fi

echo ""

if [ $starting_point -lt 2 ]; then
  condor_done=0
  while [ $condor_done -lt 1 ]; do
    condor_done=`condor_status -submitters | grep danielj | tail -n 1 | awk '{if ($2 == 0 && $3 == 0) {print 1} else {print 0}}'`
    running=`condor_status -submitters | grep danielj | tail -n 1 | awk '{print $2}'`
    waiting=`condor_status -submitters | grep danielj | tail -n 1 | awk '{print $3}'`
    echo "Waiting for condor to finish($running jobs running, $waiting jobs waiting)...`date`"
    if [ $condor_done -lt 1 ]; then
      sleep 30
    fi
  done

  echo ""
  echo pausing 5 sec...
  sleep 5
  echo ""
fi

#this puts lots of output files on the scratch directory that need to be combined. the default output location is:
# /direct/phenix+scratch/danielj/condor/get_phys_dists/run13data/trees

#combine the condor get_phys_dists output files into one file with the hadd_combine_files.sh script found here:
# /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/combine_files
if [ $starting_point -lt 1 ]; then
  echo "hadd_combining data phys_dist output...`date`"
  ${code_root_dir}/utility/combine_files/hadd_combine_files.sh \
    /direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists \
    run13_data_trees \
    /direct/phenix+scratch/danielj/condor/get_phys_dists/run13data/trees \
    tree
  echo"Finished combining data phys_dists...`date`"
fi

if [ $starting_point -lt 2 ]; then
  echo "hadd_combining OT phys_dist output...`date`"
  ${code_root_dir}/utility/combine_files/hadd_combine_files.sh \
    /direct/phenix+spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists \
    run13_data_ot_trees \
    /direct/phenix+scratch/danielj/condor/get_phys_dists/run13dataot/trees \
    tree
  echo"Finished combining OT phys_dists...`date`"
  
  echo ""
  echo pausing 5 sec...
  sleep 5
  echo ""
fi


#now combine the normal data set with the other trigger data set with the merge_other_triggers.C macro, 
# run by the run_merge_other_triggs.sh script both found in:
# /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/combine_files
if [ $starting_point -lt 3 ]; then
  echo "Merging different trigger files...`date`"
  ${code_root_dir}/utility/combine_files/run_merge_other_triggs.sh
  echo "Finished merging triggers...`date`"

  echo ""
  echo pausing 5 sec...
  sleep 5
  echo ""
fi

#this puts the final phys_dists data file here:
# /phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/phys_dists/run13_data_total_trees_combined_test.root

#now do wness preselection by running the following macros:

if [ $starting_point -lt 4 ]; then
  echo "Starting wness calculations...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute_totaltrig.C\(${which_year},1,${which_year}\)
  echo "wness calculation 1 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute_totaltrig.C\(${which_year},1,1\)
  echo "wness calculation 2 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,0\)
  echo "wness calculation 3 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,1\)
  echo "wness calculation 4 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,2\)
  echo "wness calculation 5 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,3\)
  echo "wness calculation 6 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,4\)
  echo "wness calculation 7 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,5\)
  echo "wness calculation 8 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,6\)
  echo "wness calculation 9 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,7\)
  echo "wness calculation 10 of 11 done...`date`"
  root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_muonbkg_totaltrig_execute.C\(${which_year},1,8\)
  echo "All wness calculations done!...`date`"
fi
#do eta vs dw23 fit




#!/bin/bash

which_year=$1

code_root_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis"

#do wness preselection

root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},0,${which_year}\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},0,0\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},1,${which_year}\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},1,1\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},2,${which_year}\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute.C\(${which_year},2,2\)

#do eta vs dw23 fit




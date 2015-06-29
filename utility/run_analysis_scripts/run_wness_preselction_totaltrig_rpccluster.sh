#!/bin/bash


code_root_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis"

#do wness preselection

root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute_totaltrig_rpccluster.C\(1,13\)
root -b -q ${code_root_dir}/wness_preselection/run_code_scripts/wness_preselection_execute_totaltrig_rpccluster.C\(1,1\)

#do eta vs dw23 fit




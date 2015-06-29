#!/bin/bash

code_root_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis"

normal_tree="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_data_trees_combined.root"
other_tree="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_dataot_trees_combined.root"
tree_outfile="/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/rpc_cluster_phys_dists/run13_data_total_trees_combined.root"

root -l -b -q ${code_root_dir}/utility/combine_files/merge_other_triggers_rpccluster.C+\(\"$normal_tree\",\"$other_tree\",\"$tree_outfile\"\)

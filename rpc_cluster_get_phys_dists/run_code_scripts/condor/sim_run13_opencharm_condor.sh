#!/bin/bash
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""

base_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/rpc_cluster_get_phys_dists/run_code_scripts/condor"

filelist="${base_dir}/filelists/filelist_run13_py100_sim393888_opencharm.txt"
logdir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13_sim_opencharm/logs"
mkdir -p $logdir
histdir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13_sim_opencharm/hists"
mkdir -p $histdir
treedir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run13_sim_opencharm/trees"
mkdir -p $treedir

count=0
count_limit=35
#/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/rpc_cluster_get_phys_dists/run_code_scripts/condor/generate_condorcmd.sh $count_limit
args=""
label=""
tmp=""
for file in `cat $filelist`
do
  echo "file: "$file
  tmp=`basename $file .root`
  echo "run: "$tmp
  histoutfile="$histdir/hist_$tmp.root"
  echo "output files:\n$histoutfile"
  treeoutfile="$treedir/tree_$tmp.root"
  echo "$treeoutfile"

  count=$((count+1))
  args="$args$file $histoutfile $treeoutfile "

  if [ "$count" -eq "1" ]; then
    label="$tmp-"
  fi

  if [ "$count" -eq "$count_limit" ]; then
    label="$label$tmp"
    condor_submit \
      -a "output = $logdir/$label.out" \
      -a "error  = $logdir/$label.err" \
      -a "Log    = $logdir/$label.log" \
      -a "Arguments = $args" \
      ${base_dir}/condor.job
    
    count=0
    args=""
  fi

done

if [ "$count" -gt "0" ]; then

  label="$label$tmp"
  condor_submit \
    -a "output = $logdir/$label.out" \
    -a "error  = $logdir/$label.err" \
    -a "Log    = $logdir/$label.log" \
    -a "Arguments = $args" \
    ${base_dir}/condor.job

fi

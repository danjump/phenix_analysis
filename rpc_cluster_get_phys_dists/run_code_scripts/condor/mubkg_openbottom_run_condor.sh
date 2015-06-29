#!/bin/bash

run_year="12"
run_over_which="sim367593_openbottom"
filelist="filelists/filelist_run${run_year}_${run_over_which}.txt"
logdir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run${run_year}_${run_over_which}/logs"
mkdir -p $logdir
histdir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run${run_year}_${run_over_which}/hists"
mkdir -p $histdir
treedir="/direct/phenix+scratch/danielj/condor/rpc_cluster_get_phys_dists/run${run_year}_${run_over_which}/trees"
mkdir -p $treedir

count=0
count_limit=25
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
      condor.job
    
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
    condor.job

fi

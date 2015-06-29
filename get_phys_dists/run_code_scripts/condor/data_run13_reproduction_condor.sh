#!/bin/bash
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""

base_dir="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/run_code_scripts/condor"

filelist="${base_dir}/filelists/filelist_run13pp510_data_reproduction.txt"
logdir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13datarp/logs"
mkdir -p $logdir
histdir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13datarp/hists"
mkdir -p $histdir
treedir="/direct/phenix+scratch/danielj/condor/get_phys_dists/run13datarp/trees"
mkdir -p $treedir

count=0
count_limit=35
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

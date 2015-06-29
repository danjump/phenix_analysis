#!/bin/bash

count_limit=${1}
#Manually create a condor.cmd file corresponding to the desired count_limit
end=$[3*$count_limit]

cmd_file="/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/run_code_scripts/condor/condor.cmd"
rm $cmd_file

echo "#!/bin/bash" >> $cmd_file
echo "" >> $cmd_file
echo "macro=\"/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/get_phys_dists/run_code_scripts/get_phys_dists_execute.C\"" >> $cmd_file

for ((i=3; i<=$end; i=i+3))
do
  echo "" >> $cmd_file
  in=$[$i-2]
  hout=$[i-1]
  echo "infile=\${${in}}" >> $cmd_file
  echo "hist_outfile=\${${hout}}" >> $cmd_file
  echo "tree_outfile=\${${i}}" >> $cmd_file

  echo "" >> $cmd_file
  echo "root -l -b -q \$macro\\(\\\"\$infile\\\",\\\"\$hist_outfile\\\",\\\"\$tree_outfile\\\"\\)" >> $cmd_file
done

chmod a+x $cmd_file

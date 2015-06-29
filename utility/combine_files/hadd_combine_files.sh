#!/bin/bash

final_outfile_path=$1
outfile_base=$2
ls_dir=$3
search_prefix=$4

final_outfilebase="${final_outfile_path}/${outfile_base}"
temp_outfilebase="/direct/phenix+scratch/danielj/hadd_temp/${outfile_base}"

files_per_hadd=200

filecount=0
hadd_group_count=0
hadd_groups=0
hadd_file_list=""
temp_part_file_list=""
#ls_command="ls $ls_dir/$search_prefix*.root"
ls_command="find $ls_dir -name $search_prefix*.root"

for file in `$ls_command`
do
  filecount=$(($filecount + 1))
  hadd_group_count=$(($hadd_group_count + 1))
  hadd_file_list="$hadd_file_list $file"
  
  if [ "$hadd_group_count" -eq "$files_per_hadd" ]; then
    echo "hadding $hadd_group_count files in the '$hadd_groups'th group"
    echo ""
    temp_part_filename="${temp_outfilebase}_temp_part${hadd_groups}.root"
    hadd -f $temp_part_filename $hadd_file_list
    temp_part_file_list="$temp_part_file_list $temp_part_filename"
    
    hadd_groups=$(($hadd_groups+1))
    hadd_group_count=0
    hadd_file_list=""
    echo ""
  fi
  
done

if [ "$hadd_group_count" -gt "0" ]; then
  echo "hadding $hadd_group_count files in the '$hadd_groups'th (last) group"
  echo ""
  temp_part_filename="${temp_outfilebase}_temp_part${hadd_groups}.root"
  hadd -f $temp_part_filename $hadd_file_list
  temp_part_file_list="$temp_part_file_list $temp_part_filename"
  echo ""
fi

echo ""
echo "done with initial hadding, now combining those files..."

temp_part_filename="${final_outfilebase}_combined.root"
hadd -f $temp_part_filename $temp_part_file_list

echo "done! removing temporary files..."

rm ${temp_outfilebase}_temp_part*

echo "all finished!! ^_^"


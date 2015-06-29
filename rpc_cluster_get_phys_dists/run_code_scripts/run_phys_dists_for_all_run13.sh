#!/bin/bash

run13_data_dir=""
temp_ouput_dir=""
combined_output_dir=""

ls_criteria="${run13_data_dir}/*.root | grep -v se-"
outfilebase=$1
ls_selection=$2

files_per_hadd=60

filecount=0
hadd_group_count=0
hadd_groups=0
hadd_file_list=""
temp_part_file_list=""
ls_command="ls ${ls_selection}"

for file in `$ls_command`
do
  filecount=$(($filecount + 1))
  hadd_group_count=$(($hadd_group_count + 1))
  hadd_file_list="$hadd_file_list $file"

  if [ "$hadd_group_count" -eq "$files_per_hadd" ]; then
  echo "hadding $hadd_group_count files in the '$hadd_groups'th group"
  echo ""
  temp_part_filename="${outfilebase}_temp_part${hadd_groups}.root"
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
  temp_part_filename="${outfilebase}_temp_part${hadd_groups}.root"
  hadd -f $temp_part_filename $hadd_file_list
  temp_part_file_list="$temp_part_file_list $temp_part_filename"
  echo ""
fi

echo ""
echo "done with initial hadding, now combining those files..."

temp_part_filename="${outfilebase}_combined.root"
hadd -f $temp_part_filename $temp_part_file_list

echo "done! removing temporary files..."

rm ${outfilebase}_temp_part*

echo "all finished!! ^_^"

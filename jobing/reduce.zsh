#!/bin/zsh


job=$1
dtags_file=$2

job_dir=test/tests/outdir_"$job"/


mkdir $job_dir/merged-sets

echo hadd-ing jobs of $dtags_file dtags in $job_dir to merged-sets/


for d in `cat $dtags_file`                                                          
do
  echo $d
  hadd $job_dir/merged-sets/$d.root $job_dir/"$d"*root
done


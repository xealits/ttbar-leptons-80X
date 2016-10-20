#!/bin/zsh

job=$1
dtags_file=$2

merged_dir=test/tests/outdir_test_v11_ttbar_"$job"/merged-sets/



find_existing_template_dtag() {                                                        
  job_dir=$2

  for d in `cat $1`
  do
    if [ -e $job_dir/$d.root ]
    then
        echo $d
        break
    fi
  done
}


distrs=${*:3}
for distr in ${*:3}
do
  echo extracting $distr from $merged_dir
  templ_dtag=`find_existing_template_dtag $dtags_file $merged_dir`
  echo template dtag is $templ_dtag
  root -l -q 'job-reduce/csv_distr_TH1D.cc("'$distr'", "'$merged_dir'", "'$templ_dtag'", true)' > "$job"_$distr

  for d in `cat $dtags_file`
  do
    if [ -e "$merged_dir/$d.root" ]
    then
        echo $d
        root -l -q 'job-reduce/csv_distr_TH1D.cc("'$distr'", "'$merged_dir'", "'$d'", false)' >> "$job"_$distr
    fi
  done
done


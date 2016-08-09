#!/bin/bash


raws_dir="$1"



#distr_names=`cat distrs`
#counters_names=`cat counters counters.multi*`

echo making job reduce directories


# DIRECTORIES
mkdir $raws_dir/merged-sets
mkdir $raws_dir/merged-sets/jobsums

# SCRIPTS
#cp job-reduce/processing-counters_merging-sets.R      $raws_dir/merged-sets/jobsums
#cp job-reduce/processing-distributions_merging-sets.R $raws_dir/merged-sets/jobsums
#cp job-reduce/processing_stitch_weightflow.R          $raws_dir/merged-sets/jobsums

# No copy, access the scripts via path in env var, obtained by readlink -e

SCRIPTS_SUM_COUNTERS=`readlink -e job-reduce/processing-counters_merging-sets.R`
SCRIPTS_SUM_DISTR=`readlink -e job-reduce/processing-distributions_merging-sets.R`
SCRIPTS_STITCH_WEIGHTFLOW=`readlink -e job-reduce/processing_stitch_weightflow.R`


# THE OUTPUT DATA
#cp job-reduce/distrs_* job-reduce/counters counters.multi* $raws_dir/merged-sets/
#cp job-reduce/distrs_* job-reduce/counters counters.multi* $raws_dir/merged-sets/jobsums

# No copy, access output names via env var:
# store full path to the job-reduce/ dir and access files in there

JOBREDUCE_DIR=`readlink -e job-reduce/`
#DISTRS=`cat job-reduce/distrs_*`
#COUNTERS=`cat job-reduce/counters`
#COUNTERS_MULTISELECT=`cat counters.multi*`


# +++++++++ merge procedures
function Merge_distr {
   # Description: 
   echo $1
   distrs_headername=$2
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $1
   grep --no-filename "^$1:content," ../*csv | sed "s/^$1:content,//" >> $1 
   echo "[Merged distr]"
}  

function Merge_counter {
   # Description: 
   echo $1
   echo "type,dtag,job_num,$1" > $1
   grep --no-filename "^$1," ../*csv | sed "s/^$1,//" >> $1 
   echo "[Merged counter]"
}



# +++++++++ sum procedures
function Sum_distr {
   # Description:
   echo $1
   input=../"$1"
   output="$1"
   # ./processing-distributions_merging-sets.R $input $output
   $SCRIPTS_SUM_DISTR $input $output
   echo "[Summed distr]"
}

function Sum_counter {
   # Description:
   echo $1
   input=../"$1"
   output="$1"
   #./processing-counters_merging-sets.R $input $output
   $SCRIPTS_SUM_COUNTERS $input $output
   echo "[Summed counter]"
}



# ------------------------- MERGE

cd $raws_dir/merged-sets
echo merging jobs in $raws_dir

echo Pile-Up Distributions


distrs=`cat $JOBREDUCE_DIR/distrs_pileup`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_pileup_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done








echo




# ||||||||||||||||||||||||||||||||| SUM


cd jobsums
echo summing jobs in $raws_dir/merged-sets

echo Pile-Up Distributions


for t in `cat  $JOBREDUCE_DIR/distrs_pileup `
do
   Sum_distr $t
done

echo







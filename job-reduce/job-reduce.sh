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
   # 2 input parameters
   #   $1 distr name
   #   $2 distr header name
   distrs_headername=$2

   { # adding the red collor to error stream
   echo Merge_distr $1 $2
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $1
   grep --no-filename "^$1:content," ../*csv | sed "s/^$1:content,//" >> $1 
   echo "[Merged distr]"
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}  

function Merge_counter {
   # Description: 
   # 1 input parameter
   #   $1 counter name

   { # adding the red collor to error stream
   echo Merge_counter $1
   echo "type,dtag,job_num,$1" > $1
   grep --no-filename "^$1," ../*csv | sed "s/^$1,//" >> $1 
   echo "[Merged counter]"
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}



# +++++++++ sum procedures
function Sum_distr {
   # Description:
   echo Sum_distr $1
   input=../"$1"
   output="$1"

   { # adding the red collor to error stream
   # ./processing-distributions_merging-sets.R $input $output
   $SCRIPTS_SUM_DISTR $input $output
   echo "[Summed distr]"
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}

function Sum_counter {
   # Description:
   echo Sum_counter $1
   input=../"$1"
   output="$1"

   { # adding the red collor to error stream
   #./processing-counters_merging-sets.R $input $output
   $SCRIPTS_SUM_COUNTERS $input $output
   echo "[Summed counter]"
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}


# ============ both procedures

function Reduce_distr {
   # Description: 
   # 3 input parameters
   #   $1 counter name
   #   $2 distr header name
   #   $3 directory with jobs outputs

   { # adding the red collor to error stream
   echo in $3 reducing $1
   cd $3/merged-sets
   Merge_distr $1 $2
   cd jobsums
   Sum_distr $1
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}

function Reduce_counter {
   # Description: 
   # 2 input parameters
   #   $1 counter name
   #   $2 directory with jobs outputs

   { # adding the red collor to error stream
   echo in $2 reducing $1
   cd $2/merged-sets
   Merge_counter $1
   cd jobsums
   Sum_counter $1
   } 2> >(while read line; do echo -e "\e[01;31m$line\e[0m" >&2; done)
}




# ------------------------- MERGE

cd $raws_dir/merged-sets
echo merging jobs in $raws_dir

echo Distributions

# messy merging now
# distrs_eta_headername  
#distrs_n
#distrs_n_headername
#distrs_pileup
#distrs_pileup_headername
#distrs_pt_headername
#distrs_pt

# header format in output files:
#header,pu, ...

# sed "s/^$t:header,[^,]*,[^,]*,[^,]*,/type,dtag,job_num,/"

#for t in `cat distrs`

distrs=`cat $JOBREDUCE_DIR/distrs_pt`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_pt_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done


distrs=`cat $JOBREDUCE_DIR/distrs_energy`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_energy_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done

distrs=`cat $JOBREDUCE_DIR/distrs_eta`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_eta_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done





distrs=`cat $JOBREDUCE_DIR/distrs_n`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_n_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done


distrs=`cat $JOBREDUCE_DIR/distrs_pileup`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_pileup_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done



distrs=`cat $JOBREDUCE_DIR/distrs_particleids`
distrs_headername=`cat $JOBREDUCE_DIR/distrs_particleids_headername`
echo $distrs_headername

for t in $distrs
do
  Merge_distr $t $distrs_headername
done






echo


echo Counters

for t in `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multi*`
do
  Merge_counter $t
done

echo





# ||||||||||||||||||||||||||||||||| SUM


cd jobsums
echo summing jobs in $raws_dir/merged-sets

echo Distributions


for t in `cat $JOBREDUCE_DIR/distrs_pt $JOBREDUCE_DIR/distrs_n $JOBREDUCE_DIR/distrs_pileup $JOBREDUCE_DIR/distrs_energy $JOBREDUCE_DIR/distrs_eta`
do
   Sum_distr $t
done


for t in `cat  $JOBREDUCE_DIR/distrs_eta`
do
   Sum_distr $t
done


echo




echo Counters

for t in `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multi*`
do
   Sum_counter $t
done

echo





# ................................... STITCH WEIGHT_FLOW

echo stitching weight_flows

echo SingleElectron
$SCRIPTS_STITCH_WEIGHTFLOW weightflow_el `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multiselect.el`

echo SingleMuon
$SCRIPTS_STITCH_WEIGHTFLOW weightflow_mu `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multiselect.mu`

echo EMu
$SCRIPTS_STITCH_WEIGHTFLOW weightflow_emu `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multiselect.emu`

echo MuMu
$SCRIPTS_STITCH_WEIGHTFLOW weightflow_mumu `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multiselect.mumu`

echo EE
$SCRIPTS_STITCH_WEIGHTFLOW weightflow_ee `cat $JOBREDUCE_DIR/counters $JOBREDUCE_DIR/counters.multiselect.ee`




#!/bin/bash


raws_dir="$1"



#distr_names=`cat distrs`
#counters_names=`cat counters counters.multi*`

echo making job reduce directories


# DIRECTORIES
mkdir $raws_dir/merged-sets
mkdir $raws_dir/merged-sets/jobsums



# source definitions
. job-reduce/reduce_defs.sh



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




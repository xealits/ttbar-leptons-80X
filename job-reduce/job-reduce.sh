#!/bin/bash

raws_dir="$1"



#distr_names=`cat distrs`
#counters_names=`cat counters counters.multi*`

echo making job reduce directories


# DIRECTORIES
mkdir $raws_dir/merged-sets
mkdir $raws_dir/merged-sets/jobsums

# SCRIPTS
#cp merge_jobs.sh $raws_dir/merged-sets/
#cp sum_jobs.sh $raws_dir/merged-sets/jobsums
cp processing-counters_merging-sets.R      $raws_dir/merged-sets/jobsums
cp processing-distributions_merging-sets.R $raws_dir/merged-sets/jobsums
cp processing_stitch_weightflow.R          $raws_dir/merged-sets/jobsums
#cp processing-counters_merging-multiselect.R $raws_dir/merged-sets/jobsums
#cp stitch_weightflows.sh $raws_dir/merged-sets/jobsums

# THE OUTPUT DATA
cp distrs_* counters counters.multi* $raws_dir/merged-sets/
cp distrs_* counters counters.multi* $raws_dir/merged-sets/jobsums
#cp dtags $raws_dir/merged-sets/jobsums






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

distrs=`cat distrs_pt`
distrs_headername=`cat distrs_pt_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done


distrs=`cat distrs_energy`
distrs_headername=`cat distrs_energy_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done

distrs=`cat distrs_eta`
distrs_headername=`cat distrs_eta_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done





distrs=`cat distrs_n`
distrs_headername=`cat distrs_n_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done


distrs=`cat distrs_pileup`
distrs_headername=`cat distrs_pileup_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done



distrs=`cat distrs_particleids`
distrs_headername=`cat distrs_particleids_headername`
echo $distrs_headername

for t in $distrs
do
   echo $t
   echo "type,dtag,job_num," `grep -m 1 --no-filename "^header,$distrs_headername," ../*.csv | head -n 1 | sed "s/^header,$distrs_headername,//"` > $t
   grep --no-filename "^$t:content," ../*csv | sed "s/^$t:content,//" >> $t 
   echo "[Merged]"
done






echo


echo Counters

for t in `cat counters counters.multi*`
do
   echo $t
   echo "type,dtag,job_num,$t" > $t
   grep --no-filename "^$t," ../*csv | sed "s/^$t,//" >> $t 
   echo "[Merged]"
done

echo





# ||||||||||||||||||||||||||||||||| SUM

cd jobsums
echo summing jobs in $raws_dir/merged-sets

echo Distributions


for t in `cat distrs_pt distrs_n distrs_pileup distrs_energy distrs_eta`
do
   echo $t
   input=../"$t"
   output="$t"
   ./processing-distributions_merging-sets.R $input $output
   echo "[Merged jobs]"
done

echo




echo Counters

for t in `cat counters counters.multi*`
do
   echo $t
   input=../"$t"
   output="$t"
   ./processing-counters_merging-sets.R $input $output
   echo "[Merged jobs]"
done

echo





# ................................... STITCH WEIGHT_FLOW

echo stitching weight_flows
echo SingleElectron
./processing_stitch_weightflow.R weightflow_el `cat counters counters.multiselect.el`

echo SingleMuon
./processing_stitch_weightflow.R weightflow_mu `cat counters counters.multiselect.mu`

echo EMu
./processing_stitch_weightflow.R weightflow_emu `cat counters counters.multiselect.emu`

echo MuMu
./processing_stitch_weightflow.R weightflow_mumu `cat counters counters.multiselect.mumu`

echo EE
./processing_stitch_weightflow.R weightflow_ee `cat counters counters.multiselect.ee`


#counters  counters.multisel.el  counters.multisel.mu  counters.multiselect.ee  counters.multiselect.emu  counters.multiselect.mumu


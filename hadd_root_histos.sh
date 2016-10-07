mc_dtags="MC2016_noHLT_WJets_amcatnlo MC2016_noHLT_DYJetsToLL_10to50_amcatnlo MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo MC2016_noHLT_QCD_HT-100-200 MC2016_noHLT_QCD_HT-200-300 MC2016_noHLT_QCD_HT-300-500 MC2016_noHLT_QCD_HT-500-700 MC2016_noHLT_QCD_HT-700-1000 MC2016_noHLT_QCD_HT-1000-1500 MC2016_noHLT_QCD_HT-1500-2000 MC2016_noHLT_QCD_HT-2000-Inf"
data_dtags="Data13TeV_JetHT2016D_PromptReco_v2 Data13TeV_SingleMuon2016D_PromptRecoV2"
dtags="$data_dtags $mc_dtags"

task=j5.6
task_dir=/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/tests/outdir_test_jettaufakes_"$task"/

#hadd j5.6_MC2016_noHLT_WJets.root /afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/tests/outdir_test_jettaufakes_j5.6/MC2016_noHLT_WJets_amcatnlo_*root
#hadd j5.6_MC2016_noHLT_DTJets_10to50.root /afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/tests/outdir_test_jettaufakes_j5.6/MC2016_noHLT_DYJetsToLL_10to50_amcatnlo_*root
#hadd j5.6_MC2016_noHLT_DTJets_50toInf.root /afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/tests/outdir_test_jettaufakes_j5.6/MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo_*root

echo "from $task_dir"

for dtag in $data_dtags
do
  echo "[ $dtag ]"
  hadd "$task"_"$dtag".root $task_dir/$dtag*root
done


# imputs: job_output_directory dtag
# output: job_output_directory/summary_job -- merged dtag jobs in the output directory

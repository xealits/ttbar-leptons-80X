# maybe I just should make a bunch of scripts with these
# just..lets keep it together here for now

#just in case
SHELL:=/bin/bash

reduce_ntuples:
	# job output filenames:
	# outdir/v10.8/Data13TeV_SingleElectron2016B_03Feb2017_ver2_0._.root
	# no need for mc_channel in ntuples
	#mkdir outdir/${job}/merged-sets
	echo hadd-ing jobs of jobing/dtags/${dtags} dtags in outdir/${job} to merged-sets/
	for pair in `cat jobing/dtags/${dtags}`; \
	do                                           \
	  echo $$pair - `echo $$pair | cut -d, -f1` _ `echo $$pair | cut -d, -f2` ;               \
	  hadd outdir/${job}/merged-sets/"`echo $$pair | cut -d, -f1`".root outdir/${job}/"`echo $$pair | cut -d, -f1`"*[0-9].root ; \
	done

reduce:
	# job output filenames:
	# outdir/v10.8/Data13TeV_SingleElectron2016B_03Feb2017_ver2_0._.root
	#mkdir outdir/${job}/merged-sets
	echo hadd-ing jobs of jobing/dtags/${dtags} dtags in outdir/${job} to merged-sets/
	for pair in `cat jobing/dtags/${dtags}`; \
	do                                           \
	  echo $$pair - `echo $$pair | cut -d, -f1` _ `echo $$pair | cut -d, -f2` ;                                     \
	  if ls -U outdir/${job}/"`echo $$pair | cut -d, -f1`"*.`echo $$pair | cut -d, -f2`.root 1> /dev/null ;      \
	  then   \
	     hadd outdir/${job}/merged-sets/"`echo $$pair | cut -d, -f1`"_`echo $$pair | cut -d, -f2`.root outdir/${job}/"`echo $$pair | cut -d, -f1`"*.`echo $$pair | cut -d, -f2`.root ; \
	  fi ; \
	done

#for new job format:
#$ for d in `cat jobing/dtags/dtag_xsecs_ttbar`                                                                                                                                                                                    
#do
#  dtag=`echo $d | sed "s/,.*$//"`
#  echo $dtag
#  chan=`echo $d | sed "s/.*,//"`
#  echo $chan
#  hadd outdir/v10.8/merged-sets/"$dtag"_"$chan".root outdir/v10.8/$dtag*$chan*root
#done


extract_all: extract_eventflow extract_weightflow_NOMINAL extract_jet_origins extract_taujet_origins extract_fakerates extract_pileup
	# DONE

extract_eventflow:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} eventflow_mu eventflow_el eventflow_elmu eventflow_elel eventflow_mumu

extract_weightflow_systematics:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_PU_DOWN weightflow_el_PU_DOWN weightflow_elmu_PU_DOWN weightflow_elel_PU_DOWN weightflow_mumu_PU_DOWN
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_PU_UP weightflow_el_PU_UP weightflow_elmu_PU_UP weightflow_elel_PU_UP weightflow_mumu_PU_UP
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_TOP_PT weightflow_el_TOP_PT weightflow_elmu_TOP_PT weightflow_elel_TOP_PT weightflow_mumu_TOP_PT
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_JES_UP weightflow_el_JES_UP weightflow_elmu_JES_UP weightflow_elel_JES_UP weightflow_mumu_JES_UP
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_JES_DOWN weightflow_el_JES_DOWN weightflow_elmu_JES_DOWN weightflow_elel_JES_DOWN weightflow_mumu_JES_DOWN
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_JER_UP weightflow_el_JER_UP weightflow_elmu_JER_UP weightflow_elel_JER_UP weightflow_mumu_JER_UP
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_JER_DOWN weightflow_el_JER_DOWN weightflow_elmu_JER_DOWN weightflow_elel_JER_DOWN weightflow_mumu_JER_DOWN
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_BTAG_UP weightflow_el_BTAG_UP weightflow_elmu_BTAG_UP weightflow_elel_BTAG_UP weightflow_mumu_BTAG_UP
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_BTAG_DOWN weightflow_el_BTAG_DOWN weightflow_elmu_BTAG_DOWN weightflow_elel_BTAG_DOWN weightflow_mumu_BTAG_DOWN
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_MET_UNCLUSTERED_EN_UP weightflow_el_MET_UNCLUSTERED_EN_UP weightflow_elmu_MET_UNCLUSTERED_EN_UP weightflow_elel_MET_UNCLUSTERED_EN_UP weightflow_mumu_MET_UNCLUSTERED_EN_UP
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_MET_UNCLUSTERED_EN_DOWN weightflow_el_MET_UNCLUSTERED_EN_DOWN weightflow_elmu_MET_UNCLUSTERED_EN_DOWN weightflow_elel_MET_UNCLUSTERED_EN_DOWN weightflow_mumu_MET_UNCLUSTERED_EN_DOWN

extract_weightflow_NOMINAL:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu_NOMINAL weightflow_el_NOMINAL weightflow_elmu_NOMINAL weightflow_elel_NOMINAL weightflow_mumu_NOMINAL

extract_weightflow:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu weightflow_el weightflow_elmu weightflow_elel weightflow_mumu

extract_jet_origins:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} singlemu_pretau_jet_partonFlavour singleel_pretau_jet_partonFlavour

extract_taujet_origins:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} singlemu_pretau_taujet_origins singleel_pretau_taujet_origins

extract_all_jet_origins:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} singlemu_pretau_jet_partonFlavour singleel_pretau_jet_partonFlavour \
	singlemu_passjet_jet_partonFlavour      \
	singleel_passjet_jet_partonFlavour      \
	dilep_pass2leps_jet_partonFlavour       \
	dilep_pass2leps_jet_hadronFlavour       \
	dilep_passjets_jet_partonFlavour        \
	dilep_passjets_jet_hadronFlavour        \
	dilep_passjetsNbtag_jet_partonFlavour   \
	dilep_passjetsNbtag_jet_hadronFlavour   \
	dilep_passbtagfinal_jet_partonFlavour   \
	dilep_passbtagfinal_jet_hadronFlavour   \
	elel_passjets_jet_partonFlavour         \
	elel_passjets_jet_hadronFlavour         \
	elel_passjetsNbtag_jet_partonFlavour    \
	elel_passjetsNbtag_jet_hadronFlavour    \
	elel_passbtagfinal_jet_partonFlavour    \
	elel_passbtagfinal_jet_hadronFlavour    \
	mumu_passjets_jet_partonFlavour         \
	mumu_passjets_jet_hadronFlavour         \
	mumu_passjetsNbtag_jet_partonFlavour    \
	mumu_passjetsNbtag_jet_hadronFlavour    \
	mumu_passbtagfinal_jet_partonFlavour    \
	mumu_passbtagfinal_jet_hadronFlavour    \
	elmu_passjets_jet_partonFlavour         \
	elmu_passjets_jet_hadronFlavour         \
	elmu_passjetsNbtag_jet_partonFlavour    \
	elmu_passjetsNbtag_jet_hadronFlavour    \
	elmu_passbtagfinal_jet_partonFlavour    \
	elmu_passbtagfinal_jet_hadronFlavour
	#elmu_passbtagfinal_jet_genParton_pdgId
	#dilep_pass2leps_jet_genParton_pdgId
	#elmu_passjetsNbtag_jet_genParton_pdgId
	#elmu_passjets_jet_genParton_pdgId
	#mumu_passbtagfinal_jet_genParton_pdgId
	#mumu_passjetsNbtag_jet_genParton_pdgId
	#mumu_passjets_jet_genParton_pdgId
	#elel_passbtagfinal_jet_genParton_pdgId
	#elel_passjetsNbtag_jet_genParton_pdgId
	#elel_passjets_jet_genParton_pdgId
	#dilep_passbtagfinal_jet_genParton_pdgId
	#dilep_passjets_jet_genParton_pdgId
	#dilep_passjetsNbtag_jet_genParton_pdgId


extract_fakerates:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} singleel_pretauselection_jettotaufakerates singlemu_pretauselection_jettotaufakerates

extract_pileup:
	extract_long_csv pileup_beforetrig_num_inters_rawWeight pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_beforetrig_num_inters_rawWeight.long.csv
	extract_long_csv pileup_beforetrig_num_inters_weight    pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_beforetrig_num_inters_weight.long.csv
	extract_long_csv pileup_passtrig_num_inters_rawWeight pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_passtrig_num_inters_rawWeight.long.csv
	extract_long_csv pileup_passtrig_num_inters_weight    pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_passtrig_num_inters_weight.long.csv
	extract_long_csv pileup_beforetrig_nGoodPV_rawWeight pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_beforetrig_nGoodPV_rawWeight.long.csv
	extract_long_csv pileup_beforetrig_nGoodPV_weight    pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_beforetrig_nGoodPV_weight.long.csv
	extract_long_csv pileup_passtrig_nGoodPV_rawWeight pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_passtrig_nGoodPV_rawWeight.long.csv
	extract_long_csv pileup_passtrig_nGoodPV_weight    pileup outdir/"${job}"/merged-sets/ `cat dtags_"${job}"` > outdir/"${job}"/merged-sets/jobsums/pileup_passtrig_nGoodPV_weight.long.csv
	#./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} pileup_beforetrig_num_inters_rawWeight pileup_beforetrig_num_inters_weight pileup_beforetrig_nGoodPV_rawWeight pileup_beforetrig_nGoodPV_weight pileup_passtrig_num_inters_rawWeight pileup_passtrig_num_inters_weight pileup_passtrig_nGoodPV_rawWeight pileup_passtrig_nGoodPV_weight


get_dtags:
	ls outdir/"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > dtags_${job}

get_job_dtags:
	# job output filenames:
	# outdir/v10.8/Data13TeV_SingleElectron2016B_03Feb2017_ver2_0._.root
	ls outdir/"${job}"/*root | sed -e 's,^.*'${job}'/,,' -e 's/_[0-9][0-9]*\.[a-z]*\.root//' | sort | uniq > dtags_job_${job}

jet_fakerate_measurements:
	jet_fakerate_measurements --add-header HLTjet_qcd outdir/${jjob}/merged-sets/ `cat dtags_${jjob}_JetHT ` 2> /dev/null >  ${jjob}_${vjob}_fakerate_measurement
	jet_fakerate_measurements HLTmu_wjets   outdir/${jjob}/merged-sets/ `cat dtags_${jjob}_SingleMuon ` 2> /dev/null >> ${jjob}_${vjob}_fakerate_measurement
	jet_fakerate_measurements mumu_passjets outdir/${vjob}/merged-sets/ `cat dtags_${vjob} ` 2> /dev/null >> ${jjob}_${vjob}_fakerate_measurement
	jet_fakerate_measurements elmu_passjets outdir/${vjob}/merged-sets/ `cat dtags_${vjob} ` 2> /dev/null >> ${jjob}_${vjob}_fakerate_measurement
	jet_fakerate_measurements elel_passjets outdir/${vjob}/merged-sets/ `cat dtags_${vjob} ` 2> /dev/null >> ${jjob}_${vjob}_fakerate_measurement

jet_fakerate_measurements_control:
	jet_fakerate_measurements --add-header singlemu_pretau outdir/${job}/merged-sets/ `cat dtags_${job}` 2> /dev/null >> ${job}_fakerate_measurement_control
	jet_fakerate_measurements singleel_pretau outdir/${job}/merged-sets/ `cat dtags_${job}` 2> /dev/null >> ${job}_fakerate_measurement_control


selection_control_histos: normalize=f
selection_control_histos: logy=f
selection_control_histos:
	# stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_vlooseTaus_tau_jets_distr_pt  f "" " vlooseTaus tau p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}`
	# looks like this now
	for h in `cat cur_th1d_histos_pt_eta_energy`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singleel: normalize=f
selection_control_histos_singleel: logy=f
selection_control_histos_singleel:
	for h in `cat cur_th1d_histos_pt_eta_energy_singleel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_singleel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_elel: normalize=f
selection_control_histos_elel: logy=f
selection_control_histos_elel:
	for h in `cat cur_th1d_histos_pt_eta_energy_elel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_elel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_elmu: normalize=f
selection_control_histos_elmu: logy=f
selection_control_histos_elmu:
	for h in `cat cur_th1d_histos_pt_eta_energy_elmu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_elmu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_mumu: normalize=f
selection_control_histos_mumu: logy=f
selection_control_histos_mumu:
	for h in `cat cur_th1d_histos_pt_eta_energy_mumu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_mumu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singleel_other: normalize=f
selection_control_histos_singleel_other: logy=f
selection_control_histos_singleel_other:
	for h in `cat cur_th1d_histos_other_singleel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singlemu: normalize=f
selection_control_histos_singlemu: logy=f
selection_control_histos_singlemu:
	for h in `cat cur_th1d_histos_pt_eta_energy_singlemu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_singlemu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singlemu_other: normalize=f
selection_control_histos_singlemu_other: logy=f
selection_control_histos_singlemu_other:
	for h in `cat cur_th1d_histos_other_singlemu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta: normalize=f
selection_control_histos_eta: logy=f
selection_control_histos_eta:
	for h in `cat cur_th1d_histos_eta`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta_singleel: normalize=f
selection_control_histos_eta_singleel: logy=f
selection_control_histos_eta_singleel:
	for h in `cat cur_th1d_histos_eta_singleel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta_singlemu: normalize=f
selection_control_histos_eta_singlemu: logy=f
selection_control_histos_eta_singlemu:
	for h in `cat cur_th1d_histos_eta_singlemu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 5 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_pileup: normalize=f
selection_control_histos_pileup: logy=f
selection_control_histos_pileup:
	for h in pileup_passtrig_nGoodPV_rawWeight pileup_passtrig_nGoodPV_weight pileup_singlelepton_nGoodPV_rawWeight pileup_singlelepton_nGoodPV_weight pileup_dilepton_nGoodPV_rawWeight pileup_dilepton_nGoodPV_weight pileup_sel_passos_nGoodPV_rawWeight pileup_sel_passos_nGoodPV_weight pileup_smu_passos_nGoodPV_rawWeight pileup_smu_passos_nGoodPV_weight ;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 0 50 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_rho: normalize=f
selection_control_histos_rho: logy=f
selection_control_histos_rho:
	for h in rho_beforetrig_rawWeight rho_beforetrig_weight rho_passtrig_rawWeight rho_passtrig_weight rho_singlelepton_rawWeight rho_singlelepton_weight  rho_dilepton_rawWeight rho_dilepton_weight \
		 rhoCentral_beforetrig_rawWeight rhoCentral_beforetrig_weight rhoCentral_passtrig_rawWeight rhoCentral_passtrig_weight rhoCentral_singlelepton_rawWeight rhoCentral_singlelepton_weight  rhoCentral_dilepton_rawWeight rhoCentral_dilepton_weight \
		 rhoCentralNeutral_beforetrig_rawWeight rhoCentralNeutral_beforetrig_weight rhoCentralNeutral_passtrig_rawWeight rhoCentralNeutral_passtrig_weight rhoCentralNeutral_singlelepton_rawWeight rhoCentralNeutral_singlelepton_weight  rhoCentralNeutral_dilepton_rawWeight rhoCentralNeutral_dilepton_weight \
		 rhoCentralChargedPileUp_beforetrig_rawWeight rhoCentralChargedPileUp_beforetrig_weight rhoCentralChargedPileUp_passtrig_rawWeight rhoCentralChargedPileUp_passtrig_weight rhoCentralChargedPileUp_singlelepton_rawWeight rhoCentralChargedPileUp_singlelepton_weight  rhoCentralChargedPileUp_dilepton_rawWeight rhoCentralChargedPileUp_dilepton_weight ;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h same $$h " Events/bin " 1 -1 -1 0.25 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

stacked_fakerate_histos:
	for h in HLTjetmu_wjets_jets_distr;  \
	do                                                           \
	stacked_histo_project3toX_distr ${lumi} $$h x 1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	stacked_histo_project3toX_distr ${lumi} $$h y 1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	stacked_histo_project3toX_distr ${lumi} $$h z 1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	#jobing/stacked_histo_project3toX_distr.cc

done_ck:
	ls outdir/"${job}"/*job_done | wc -l
	ls outdir/"${job}"/*root     | wc -l
	ls outdir/"${job}"/*cfg.py   | wc -l
	ls outdir/"${job}"/Data*job_done | wc -l
	ls outdir/"${job}"/Data*root     | wc -l
	ls outdir/"${job}"/Data*cfg.py   | wc -l
	ls outdir/"${job}"/MC*job_done | wc -l
	ls outdir/"${job}"/MC*root     | wc -l
	ls outdir/"${job}"/MC*cfg.py   | wc -l
	ls outdir/"${job}"/*TTJ*job_done | wc -l
	ls outdir/"${job}"/*TTJ*root     | wc -l
	ls outdir/"${job}"/*TTJ*cfg.py   | wc -l

ck_running:
	for j in `sed 's/^.*-J \([^ ]*\).*$$/\1/' outdir/${job}/FARM/inputs/*_bsub.sh`; \
	do                                 \
	st=$$(bjobs -J $$j 2>/dev/null);     \
	if ! [[ $$st == "" ]]; then echo $$j; fi;  \
	done | tee ${job}.running_jobs

ck_errors:
	wc -l outdir/"${job}"/FARM/outputs/* | sort


.PHONY:missing_jobs
missing_jobs:
	ls outdir/${job}/*cfg.py | sed 's/_cfg.py/.job_done/' > submitted_jobs
	for j in `cat submitted_jobs`;        \
	do                                    \
	if [ ! -e $$j ]; then echo $$j; fi;   \
	done | sed -e 's/job_done/sh/' -e "s;"${job}";"${job}"//FARM/inputs;" > missing_jobs
	wc -l missing_jobs

resubmit: missing_jobs
	for j in `cat missing_jobs`;                                    \
	do                                                              \
	  grep -h $$j outdir/"${job}"//FARM/inputs/*bsub.sh; \
	done > ${job}_missing_jobs_bsubs
	wc -l ${job}_missing_jobs*
	source ${job}_missing_jobs_bsubs

.PHONY:completeness
completeness:
	printf "dtag\tjob\tjson\troot\tcfg\n"
	printf "\n"
	for d in `echo ${dtags}`; \
	do                        \
	printf "%s\t%d\t%d\t%d\t%d\n" $$d `ls outdir/${job}/"$$d"*job_done | wc -l` `ls outdir/${job}/"$$d"*json | wc -l`  `ls outdir/${job}/"$$d"*root | wc -l` `ls outdir/${job}/"$$d"*cfg.py | wc -l`;  \
	done

#job=v9.15
#outdir=outdir/${job}
#echo $outdir
#dtags="SingleElectron2016B SingleElectron2016C SingleElectron2016D SingleElectron2016E SingleElectron2016F SingleElectron2016G SingleElectron2016H"


proxy:
	# for consistency
	make-x509-proxy

calculate_luminosity:
	cat outdir/"${job}"/*${runs}*json > "${job}"_${runs}_prejson
	sed -i 's/}{/,/g'     "${job}"_${runs}_prejson
	sed -i 's/,,,*/,/g'   "${job}"_${runs}_prejson
	sed -i 's/ ,"/,\n"/g' "${job}"_${runs}_prejson
	sed -i -e 's/"[0-9]/[&/' -e 's/]]/]]]/' -e 's/{/[/' -e 's/}/]/' -e 's/:/,/' "${job}"_${runs}_prejson
	#cd analysis/job-lumis
	./merge_json.py "${job}"_${runs}_prejson "${job}"_${runs}.json
	brilcalc lumi -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i "${job}"_${runs}.json | tee "${job}"_${runs}.brilcalc.out
	mv "${job}"_${runs}_prejson analysis/job-lumis
	mv "${job}"_${runs}.json analysis/job-lumis
	mv "${job}"_${runs}.brilcalc.out analysis/job-lumis


calculate_prejson:
	#cat outdir/"${job}"/*json > "${job}"_prejson
	sed -i 's/}{/,/g'     "${prejson}"
	sed -i 's/,,,*/,/g'   "${prejson}"
	sed -i 's/ ,"/,\n"/g' "${prejson}"
	sed -i -e 's/"[0-9]/[&/' -e 's/]]/]]]/' -e 's/{/[/' -e 's/}/]/' -e 's/:/,/' "${prejson}"
	#cd analysis/job-lumis
	./merge_json.py "${prejson}" "${prejson}".json
	brilcalc lumi -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i "${prejson}".json | tee "${prejson}".brilcalc.out
	#mv "${prejson}" analysis/job-lumis
	mv "${prejson}".json analysis/job-lumis
	mv "${prejson}".brilcalc.out analysis/job-lumis

cutncount_pretau_selection_control_mu:
	#control_met_singlemu_pretau_pt_eta
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_met_singlemu_pretau_pt_eta  f " E^{miss}_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_jet_singlemu_pretau_selJetsNoLep_pt_eta  f " jet p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_jet_singlemu_pretau_selJetsNoLep_pt_eta  f " jet #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_lep_singlemu_pretau_selLeptons_pt_eta  f " lepton p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_lep_singlemu_pretau_selLeptons_pt_eta  f " lepton #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_tau_singlemu_pretau_selTausNoLep_pt_eta  f " tau p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_tau_singlemu_pretau_selTausNoLep_pt_eta  f " tau #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} singlemu_pretauselection_njetsNoLep  f " N jets " " Events/bin " 1 3 10 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 

cutncount_pretau_selection_control_el:
	#control_met_singleel_pretau_pt_eta
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_met_singleel_pretau_pt_eta  f " E^{miss}_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_jet_singleel_pretau_selJetsNoLep_pt_eta  f " jet p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_jet_singleel_pretau_selJetsNoLep_pt_eta  f " jet #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_lep_singleel_pretau_selLeptons_pt_eta  f " lepton p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_lep_singleel_pretau_selLeptons_pt_eta  f " lepton #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} x control_tau_singleel_pretau_selTausNoLep_pt_eta  f " tau p_{T} " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo2d_distr_with_ratio f ${logy} ${lumi} y control_tau_singleel_pretau_selTausNoLep_pt_eta  f " tau #eta " " Events/bin " 2 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi}     singleel_pretauselection_njetsNoLep  f " N jets " " Events/bin " 1 3 10 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 


# elmu FAKERATES direct distrs
elmu_fake_rates_jets_distrs_medium:
	# scale_width inclusive_xsec normalize csv_out logy lumi distr mc_distr name_suffix distr_name yname rebin_factor xmin xmax ymin ymax ratio_range dir zero_bin zero_bin2 x_leg dtags
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f "" " mediumTaus jet p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta f "" " mediumTaus jet #eta "  " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 

elmu_fake_rates_jets_distrs_sys_errors:
	#stacked_histo_distr_with_ratio_elmuFakeRateSys f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f " mediumTaus jet p_{T} " " Events/bin "  1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta f " mediumTaus jet #eta "  " Events/bin "  1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 

elmu_fake_rates_jets_distrs_medium_systematics:
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_PU_UP_mediumTaus_jets_distr_pt    "  mediumTaus PU UP jet p_{T} " " Events/bin "    1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_PU_UP_mediumTaus_jets_distr_eta   "  mediumTaus PU UP jet #eta "  " Events/bin "    1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_PU_DOWN_mediumTaus_jets_distr_pt  "  mediumTaus PU DOWN jet p_{T} " " Events/bin "  1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_PU_DOWN_mediumTaus_jets_distr_eta "  mediumTaus PU DOWN jet #eta "  " Events/bin "  1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_JER_UP_mediumTaus_jets_distr_pt    " mediumTaus JER UP jet p_{T} " " Events/bin "   1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_JER_UP_mediumTaus_jets_distr_eta   " mediumTaus JER UP jet #eta "  " Events/bin "   1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_JER_DOWN_mediumTaus_jets_distr_pt  " mediumTaus JER DOWN jet p_{T} " " Events/bin " 1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_JER_DOWN_mediumTaus_jets_distr_eta " mediumTaus JER DOWN jet #eta "  " Events/bin " 1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_JES_UP_mediumTaus_jets_distr_pt    " mediumTaus JES UP jet p_{T} " " Events/bin "   1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_JES_UP_mediumTaus_jets_distr_eta   " mediumTaus JES UP jet #eta "  " Events/bin "   1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_JES_DOWN_mediumTaus_jets_distr_pt  " mediumTaus JES DOWN jet p_{T} " " Events/bin " 1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta elmu_passjets_JES_DOWN_mediumTaus_jets_distr_eta " mediumTaus JES DOWN jet #eta "  " Events/bin " 1 -1 -1 0.1 outdir/${job}/merged-sets/ `cat dtags_${job}` 

elmu_fake_rates_jets_distrs_all:
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_vlooseTaus_jets_distr_pt  f " vlooseTaus jet p_{T} " " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_vlooseTaus_jets_distr_eta f " vlooseTaus jet #eta "  " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_looseTaus_jets_distr_pt   f " looseTausp jet p_{T} " " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_looseTaus_jets_distr_eta  f " looseTaus jet #eta "   " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f " mediumTaus jet p_{T} " " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta f " mediumTaus jet #eta "  " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_tightTaus_jets_distr_pt   f " tightTaus jet p_{T} "  " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_tightTaus_jets_distr_eta  f " tightTaus jet #eta "   " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_vtightTaus_jets_distr_pt  f " vtightTaus jet p_{T} " " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_vtightTaus_jets_distr_eta f " vtightTaus jet #eta "  " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 

elmu_fake_rates_taus_distrs_medium:
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_pt  f " mediumTaus tau p_{T} " " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_eta f " mediumTaus tau #eta "  " Events/bin " 1 -1 -1 0.25 outdir/${job}/merged-sets/ `cat dtags_${job}` 

elmu_fake_rates_taus_distrs_all:
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_vlooseTaus_tau_jets_distr_pt  f "" " vlooseTaus tau p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_vlooseTaus_tau_jets_distr_eta f "" " vlooseTaus tau #eta "  " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_looseTaus_tau_jets_distr_pt   f "" " looseTausp tau p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_looseTaus_tau_jets_distr_eta  f "" " looseTaus tau #eta "   " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_pt  f "" " mediumTaus tau p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_eta f "" " mediumTaus tau #eta "  " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_tightTaus_tau_jets_distr_pt   f "" " tightTaus tau p_{T} "  " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_tightTaus_tau_jets_distr_eta  f "" " tightTaus tau #eta "   " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_vtightTaus_tau_jets_distr_pt  f "" " vtightTaus tau p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_vtightTaus_tau_jets_distr_eta f "" " vtightTaus tau #eta "  " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 

elmu_fake_rates_control_plots2:
	# NOMINAL VS TOP PT trend
	# leptons
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_pt_NOMINAL      f "" " NOMINAL lep p_{T} " " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_pt_TOP_PT       f "" " TOP PT lep p_{T} "  " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	# jets
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f "" " NOMINAL jet p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_TOP_PT_mediumTaus_jets_distr_pt  f "" " TOP PT jet p_{T} " " Events/bin " 1 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 

elmu_fake_rates_control_plots:
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} met_elmu_passjets                      f "" " Missing E_{T} " " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_met_pt                   f "" " Missing E_{T} " " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} pileup_elmu_passjets_nGoodPV_rawWeight f "" " N vertices "  " Raw Events "      1 0 50 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} pileup_elmu_passjets_nGoodPV_weight    f "" " N vertices "  " Weighted Events " 1 0 50 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} rho_elmu_passjets_rawWeight f ""    " #rho "  " Events/bin "          1 0 50 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} rho_elmu_passjets_weight    f ""    " #rho "  " Weighted Events/bin " 1 0 50 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_pt  f ""  " leading lepton p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep2_pt  f ""  " second  lepton p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_pt  f ""  " leading jet p_{T} "  " Events/bin "    5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet2_pt  f ""  " second  jet p_{T} "  " Events/bin "    5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_eta f ""  " leading lepton #eta "  " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep2_eta f ""  " second  lepton #eta "  " Events/bin "  5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_eta f ""  " leading jet #eta "  " Events/bin "     5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet2_eta f ""  " second  jet #eta "  " Events/bin "     5 -1 -1 -1 -1 0.5 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 

elmu_fake_rates_control_plots_SYS_distrs:
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_pt_PU_UP  f ""  " PU UP leading lepton p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_pt_PU_UP  f ""  " PU UP leading jet p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_lep1_pt_PU_DOWN  f ""  " PU DOWN leading lepton p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_pt_PU_DOWN  f ""  " PU DOWN leading jet p_{T} "  " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_pt_JES_UP   f ""  " JES UP leading jet p_{T} "   " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet1_pt_JES_DOWN f ""  " JES DOWN leading jet p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet2_pt_JES_UP   f ""  " JES UP second jet p_{T} "   " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_jet2_pt_JES_DOWN f ""  " JES DOWN second jet p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_JES_UP_mediumTaus_jets_distr_pt   f "" " JES UP jets p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_JES_DOWN_mediumTaus_jets_distr_pt f "" " JES DOWN jets p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_JER_UP_mediumTaus_jets_distr_pt   f "" " JER UP jets p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	stacked_histo_distr_with_ratio T f f f ${logy} ${lumi} elmu_passjets_JER_DOWN_mediumTaus_jets_distr_pt f "" " JER DOWN jets p_{T} " " Events/bin " 5 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 

# elmu FAKE RATE jet distr-s
elmu_fake_rates_control_plots_SYS:
	# stacked_histo_distr_with_ratio_elmuFakeRateSys f f 35000 elmu_passjets_mediumTaus_jets_distr_pt same elmu_passjets_mediumTaus_jets " Events " 1 -1 -1 0.25 outdir/t5.0/merged-sets/ `cat dtags_t5.0 `
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f   " jet p_{T} " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta f   " jet #eta "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_vlooseTaus_tau_jets_distr_pt  f   " vLoose tau jet p_{T} " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_vlooseTaus_tau_jets_distr_eta f   " vLoose tau jet #eta "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_looseTaus_tau_jets_distr_pt   f   " Loose tau jet p_{T} "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_looseTaus_tau_jets_distr_eta  f   " Loose tau jet #eta "   " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_pt  f   " Medium tau jet p_{T} " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_eta f   " Medium tau jet #eta "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_tightTaus_tau_jets_distr_pt   f   " Tight tau jet p_{T} "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_tightTaus_tau_jets_distr_eta  f   " Tight tau jet #eta "   " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_vtightTaus_tau_jets_distr_pt  f   " vTight tau jet p_{T} " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_vtightTaus_tau_jets_distr_eta f   " vTight tau jet #eta "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 

elmu_fake_rates_control_plots_SYS_medium:
	#stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_pt  f   " Medium tau jet p_{T} " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	#stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_tau_jets_distr_eta f   " Medium tau jet #eta "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_pt  f   " Medium tau jet p_{T} " " Events/bin " 1 -1 -1 10000 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f ${logy} ${lumi} elmu_passjets_mediumTaus_jets_distr_eta f   " Medium tau jet #eta "  " Events/bin " 1 -1 -1 10000 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	# and with top pt passjets_TOP_PT_mediumTaus
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f ${logy} ${lumi} elmu_passjets_TOP_PT_mediumTaus_jets_distr_pt  f   " (with Top p_{T} weight) Medium tau jet p_{T} " " Events/bin " 1 -1 -1 10000 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f ${logy} ${lumi} elmu_passjets_TOP_PT_mediumTaus_jets_distr_eta f   " (with Top p_{T} weight) Medium tau jet #eta "  " Events/bin " 1 -1 -1 10000 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}` 

# elmu ELMU FAKE RATE
elmu_fake_rates:
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_vlooseTaus _pt  1 -1 -1 "vLooseMVAIso" " vlooseTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_vlooseTaus _eta 1 -1 -1 "vLooseMVAIso" " vlooseTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_looseTaus  _pt  1 -1 -1 "LooseMVAIso"  " looseTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_looseTaus  _eta 1 -1 -1 "LooseMVAIso"  " looseTaus  jet #eta "   outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 "MediumMVAIso" " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _eta 1 -1 -1 "MediumMVAIso" " mediumTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_tightTaus  _pt  1 -1 -1 "TightMVAIso"  " tightTaus  jet p_{T} "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_tightTaus  _eta 1 -1 -1 "TightMVAIso"  " tightTaus  jet #eta "   outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_vtightTaus _pt  1 -1 -1 "vTightMVAIso" " vtightTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_vtightTaus _eta 1 -1 -1 "vTightMVAIso" " vtightTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fake_rates_TT_tunes_nominal:
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_nominal_fullMC " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}   `
	#_taujet_origins
	#stacked_histo_distr_with_ratio f T f f f ${logy} ${lumi} elmu_passjets_mediumTaus_taujet_origins  f ""  " jet origins "  " Events/bin " 1 -1 -1 -1 -1 0.25 outdir/${job}/merged-sets/ -1 -1 0.5 `cat dtags_${job}` 
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_powheg " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_powheg   `


elmu_fake_rates_TT_tunes:
	#jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_h " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_h   `
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_isrU " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_isrU`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_isrD " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_isrD`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_fsrU " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_fsrU`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_fsrD " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_fsrD`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_hdU " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_hdU `
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_hdD " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_hdD `
	#                                                                                                                                                     _hdD 
	#jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _eta 1 -1 -1 medium_taus " mediumTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`


elmu_fake_rates_TT_tunes_batch2:
	#jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_GluonMoveCR " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_GCR`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_QCDbasedCR  " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_QCR`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_TuneM2T4_U  " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_TnUp`
	jet_1D_fakerate_per_mc ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus_TuneM2T4_D  " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}_TnDw`

# MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune
# MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune
# MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down
# MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up


elmu_fakerate_tau_candidate_energies:
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_signalCands               f   " tau Signal cand E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_signalChargedHadrCands    f   " tau Sig Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_signalNeutrHadrCands      f   " tau Sig Nu Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_signalGammaCands          f   " tau Signal Gamma E " " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_isolationCands            f   " tau Iso cand E "     " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_isolationChargedHadrCands f   " tau Iso Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_isolationNeutrHadrCands   f   " tau Iso Nu Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_tau_isolationGammaCands       f   " tau Iso Gamma E "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fakerate_tau_candidate_energies_4vars:
	# 4 iso var-s
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_photonPtSumOutsideSignalCone  f   " tau Photon Pt Sum Outside Signal Cone "    " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_chargedIsoPtSum  f   " tau Charged Iso Pt Sum "    " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_neutralIsoPtSum  f   " tau Neutral Iso Pt Sum "    " Events/bin " 1 0 40 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_mediumTaus_puCorrPtSum      f   " tau Pu Corr Pt Sum "    " Events/bin " 2 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fakerate_tau_candidate_energies_NoIsoIDtaus:
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_signalCands               f " tau Signal cand E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_signalChargedHadrCands    f " tau Sig Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_signalNeutrHadrCands      f " tau Sig Nu Hadr E "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_signalGammaCands          f " tau Signal Gamma E " " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_isolationCands            f " tau Iso cand E "     " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_isolationChargedHadrCands f " tau Iso Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_isolationNeutrHadrCands   f " tau Iso Nu Hadr E "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_tau_isolationGammaCands       f " tau Iso Gamma E "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fakerate_tau_candidate_energies_NoIsoIDtaus_4vars:
	# 4 iso var-s
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_photonPtSumOutsideSignalCone f " tau Photon Pt Sum Outside Signal Cone "  " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_chargedIsoPtSum  f   " tau Charged Iso Pt Sum "    " Events/bin " 1 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_neutralIsoPtSum  f   " tau Neutral Iso Pt Sum "    " Events/bin " 1 0 40 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys f f f f f ${lumi} elmu_passjets_NoIsoIDtaus_puCorrPtSum      f   " tau Pu Corr Pt Sum "    " Events/bin " 2 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`



elmu_fakerate_tau_candidate_energies_looseTaus:
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_signalCands               f   " tau Signal cand E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_signalChargedHadrCands    f   " tau Sig Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_signalNeutrHadrCands      f   " tau Sig Nu Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_signalGammaCands          f   " tau Signal Gamma E " " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_isolationCands            f   " tau Iso cand E "     " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_isolationChargedHadrCands f   " tau Iso Ch Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_isolationNeutrHadrCands   f   " tau Iso Nu Hadr E "  " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_tau_isolationGammaCands       f   " tau Iso Gamma E "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
elmu_fakerate_tau_candidate_energies_looseTaus_4vars:
	# 4 iso var-s
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_photonPtSumOutsideSignalCone       f   " tau Photon Pt Sum Outside Signal Cone "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_chargedIsoPtSum  f   " tau Charged Iso Pt Sum "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_neutralIsoPtSum  f   " tau Neutral Iso Pt Sum "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`
	stacked_histo_distr_with_ratio_elmuFakeRateSys T T T f f ${lumi} elmu_passjets_looseTaus_puCorrPtSum      f   " tau Pu Corr Pt Sum "    " Events/bin " 10 -1 -1 0.5 outdir/${job}/merged-sets/ `cat dtags_${job}`


elmu_fake_rates_SYS_medium:
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_mediumTaus _eta 1 -1 -1 medium_taus " mediumTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fake_rates_SYS_medium_noSYS:
	jet_1D_fakerate_per_mc_SYS f f f ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 medium_taus " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS f f f ${lumi} elmu_passjets_mediumTaus _eta 1 -1 -1 medium_taus " mediumTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`

elmu_fake_rates_SYS_all:
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_vlooseTaus _pt  1 -1 -1 "vLooseMVAIso" " vlooseTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_vlooseTaus _eta 1 -1 -1 "vLooseMVAIso" " vlooseTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_looseTaus  _pt  1 -1 -1 "LooseMVAIso"  " looseTausp jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_looseTaus  _eta 1 -1 -1 "LooseMVAIso"  " looseTaus  jet #eta "   outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_mediumTaus _pt  1 -1 -1 "MediumMVAIso" " mediumTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_mediumTaus _eta 1 -1 -1 "MediumMVAIso" " mediumTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_tightTaus  _pt  1 -1 -1 "TightMVAIso"  " tightTaus  jet p_{T} "  outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_tightTaus  _eta 1 -1 -1 "TightMVAIso"  " tightTaus  jet #eta "   outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_vtightTaus _pt  1 -1 -1 "vTightMVAIso" " vtightTaus jet p_{T} " outdir/${job}/merged-sets/ `cat dtags_${job}`
	jet_1D_fakerate_per_mc_SYS T T T ${lumi} elmu_passjets_vtightTaus _eta 1 -1 -1 "vTightMVAIso" " vtightTaus jet #eta "  outdir/${job}/merged-sets/ `cat dtags_${job}`

jet_origins:
	# "HLTjet_qcd_jet_partonFlavour HLTmu_qcd_jet_partonFlavour HLTjetmu_qcd_jet_partonFlavour HLTjet_wjets_jet_partonFlavour HLTmu_wjets_jet_partonFlavour HLTjetmu_wjets_jet_partonFlavour"
	# "dilep_passjets_jet_partonFlavour dilep_passjetsNbtag_jet_partonFlavour elel_passjets_jet_partonFlavour"
	ls outdir/"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > test_dtags
	for d in `cat test_dtags | sed -n '/^MC/p' `; do jet_origins outdir/"${job}"/merged-sets/ $$d ${channels}; done

jet_weightflow:
	# "HLTjet_qcd_jet_partonFlavour HLTmu_qcd_jet_partonFlavour HLTjetmu_qcd_jet_partonFlavour HLTjet_wjets_jet_partonFlavour HLTmu_wjets_jet_partonFlavour HLTjetmu_wjets_jet_partonFlavour"
	# "dilep_passjets_jet_partonFlavour dilep_passjetsNbtag_jet_partonFlavour elel_passjets_jet_partonFlavour"
	ls outdir/"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > test_dtags
	for d in `cat test_dtags | sed -n '/^MC/p' `; do jet_weightflow outdir/"${job}"/merged-sets/ $$d; done

weightflow_table:
	mkdir -p sys_weightflows/${job}.${dtags}/
	weightflow_table ${lumi} ${distr} outdir/${job}/merged-sets/ `cat ${dtags}` > sys_weightflows/${job}.${dtags}/${distr}

.PHONY:fakerate_distr_plots
fakerate_distr_plots:
	test -s ${file}
	# make stops if a command in the recipe fails...
	histo_project3toN_distr x "pT (GeV)" ${content_name}    ${distr} 1  0 150   ${file}
	histo_project3toN_distr y eta        ${content_name}    ${distr} 1 -2.5 2.5 ${file}
	histo_project3toN_distr z radius     ${content_name}    ${distr} 1  0   0.3 ${file}
	histo_project3toN_distr_atslice x "pT (GeV)" ${content_name}    ${distr} 1  0 150   y 2 3 ${file}

jet_distr_plots_JetHT: file=outdir/${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
jet_distr_plots_JetHT: distr=HLTjet_qcd_jets_distr
jet_distr_plots_JetHT: content_name="N jets"
jet_distr_plots_JetHT: fakerate_distr_plots

jet_distr_plots_SingleMuon: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
jet_distr_plots_SingleMuon: distr=HLTmu_wjets_jets_distr
jet_distr_plots_SingleMuon: content_name="N jets"
jet_distr_plots_SingleMuon: fakerate_distr_plots

taujet_distr_plots_JetHT: file=outdir/${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
taujet_distr_plots_JetHT: distr=HLTjet_qcd_tau_jets_distr
taujet_distr_plots_JetHT: content_name="N tau jets"
taujet_distr_plots_JetHT: fakerate_distr_plots

taujet_distr_plots_SingleMuon: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
taujet_distr_plots_SingleMuon: distr=HLTmu_wjets_tau_jets_distr
taujet_distr_plots_SingleMuon: content_name="N tau jets"
taujet_distr_plots_SingleMuon: fakerate_distr_plots

.PHONY:fakerate_ratio_plots
fakerate_ratio_plots: csv=f
fakerate_ratio_plots:
	test -s ${file}
	# make stops if a command in the recipe fails...
	histo_project3toN_ratio_distr ${csv} x pT     "fake rate" ${selection}_tau_jets_distr${largebins} ${selection}_jets_distr${largebins} 1  0 150   ${file}
	histo_project3toN_ratio_distr ${csv} y eta    "fake rate" ${selection}_tau_jets_distr${largebins} ${selection}_jets_distr${largebins} 1 -2.5 2.5 ${file}
	histo_project3toN_ratio_distr ${csv} z radius "fake rate" ${selection}_tau_jets_distr${largebins} ${selection}_jets_distr${largebins} 1  0   0.3 ${file}

.PHONY:jet_fakerate_plots_JetHT
jet_fakerate_plots_JetHT: file=outdir/${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
jet_fakerate_plots_JetHT: selection=HLTjet_qcd
jet_fakerate_plots_JetHT: largebins=""
jet_fakerate_plots_JetHT: fakerate_ratio_plots

.PHONY:jet_fakerate_plots_SingleMuon
jet_fakerate_plots_SingleMuon: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
jet_fakerate_plots_SingleMuon: selection=HLTmu_wjets
jet_fakerate_plots_SingleMuon: largebins=""
jet_fakerate_plots_SingleMuon: fakerate_ratio_plots

jet_fakerate_plots: jet_fakerate_plots_JetHT jet_fakerate_plots_SingleMu

.PHONY:jet_fakerate_plots_control_el
jet_fakerate_plots_control_el: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
jet_fakerate_plots_control_el: selection=singleel_pretau
jet_fakerate_plots_control_el: largebins=""
jet_fakerate_plots_control_el: fakerate_ratio_plots

.PHONY:jet_fakerate_plots_control_mu
jet_fakerate_plots_control_mu: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
jet_fakerate_plots_control_mu: selection=singlemu_pretau
jet_fakerate_plots_control_mu: largebins=""
jet_fakerate_plots_control_mu: fakerate_ratio_plots

jet_fakerate_plots_control: jet_fakerate_plots_control_el jet_fakerate_plots_control_mu

mumu_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
mumu_fakerate_ratio_plots: selection=mumu_passjets
mumu_fakerate_ratio_plots: largebins=""
mumu_fakerate_ratio_plots: fakerate_ratio_plots

elmu_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elmu_fakerate_ratio_plots: selection=elmu_passjets
elmu_fakerate_ratio_plots: largebins=""
elmu_fakerate_ratio_plots: fakerate_ratio_plots

elel_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elel_fakerate_ratio_plots: selection=elel_passjets
elel_fakerate_ratio_plots: largebins=""
elel_fakerate_ratio_plots: fakerate_ratio_plots

dilep_fakerate_ratio_plots: elel_fakerate_ratio_plots mumu_fakerate_ratio_plots elmu_fakerate_ratio_plots

fakerate_study_extract_all:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} HLTjet_qcd_jet_partonFlavour HLTmu_wjets_jet_partonFlavour elmu_passjets_jet_partonFlavour mumu_passjets_jet_partonFlavour

jet_mc_distrs:
	stacked_histo_project3to1_distr --normalize 1.7 T HLTjet_qcd_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 f HLTjet_qcd_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 T HLTjet_qcd_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 T HLTmu_wjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 f HLTmu_wjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 T HLTmu_wjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

tau_mc_distrs:
	stacked_histo_project2to1_distr --normalize 1.7 HLTjet_qcd_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project2to1_distr --normalize 1.7 HLTjet_qcd_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project2to1_distr --normalize 1.7 HLTmu_wjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project2to1_distr --normalize 1.7 HLTmu_wjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`


jets_taus_mc_distrs_dileptons:
	stacked_histo_project3to1_distr --normalize 1.7 T mumu_passjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 f mumu_passjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T mumu_passjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T elmu_passjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 f elmu_passjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T elmu_passjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	#
	stacked_histo_project2to1_distr --normalize 1.7 mumu_passjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 mumu_passjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 elmu_passjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 elmu_passjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`

# elmu FAKERATES from 3d distrs
jets_mc_distrs_dileptons_lumi:
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_jets_distr${large_bins} x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_jets_distr${large_bins} y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_jets_distr${large_bins} z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	# tau
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_tau_jets_distr${large_bins} x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_tau_jets_distr${large_bins} y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr ${verbose} ${lumi} ${logy} elmu_passjets_tau_jets_distr${large_bins} z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`


jets_mc_distrs_dileptons_sf:
	# Usage : jet_fakerate_sf_mc_to_data[--verbose] [--normalize] lumi distr rebin_factor x_axis_min_range x_axis_max_range name_tag dir dtags
	jet_fakerate_sf_mc_to_data ${verbose} ${lumi} elmu_passjets 1 0 250 sf outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jets_taus_mc_distrs_control:
	stacked_histo_project3to1_distr --normalize 1.7 T singlemu_pretau_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 f singlemu_pretau_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T singlemu_pretau_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T singleel_pretau_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 f singleel_pretau_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 T singleel_pretau_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	#
	stacked_histo_project2to1_distr --normalize 1.7 singlemu_pretau_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 singlemu_pretau_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 singleel_pretau_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 singleel_pretau_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`


jet_tau_origin_distrs_stacked_jets_qcd:
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 

jet_tau_origin_distrs_stacked_jets:
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd_tau z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets_tau x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets_tau y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets_tau z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 


jet_tau_origin_distrs_stacked_dileptons:
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets_tau x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets_tau y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets_tau z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets_tau x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets_tau y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets_tau z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `


jet_origin_distrs_stacked_jets_qcd:
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 


jet_origin_distrs_stacked_jets:
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTjet_qcd z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_stacked_per_origin --normalize 0.25 HLTmu_wjets z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 

jet_origin_distrs_stacked_dileptons:
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 mumu_passjets z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_stacked_per_origin --normalize 0.25 elmu_passjets z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `



jet_origin_fractions_stacked_qcd: largebins=""
jet_origin_fractions_stacked_qcd:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTjet_qcd x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTjet_qcd y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTjet_qcd z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 

jet_origin_fractions_stacked_wjets: largebins=""
jet_origin_fractions_stacked_wjets:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTmu_wjets x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTmu_wjets y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 HLTmu_wjets z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 


jet_origin_fractions_stacked_dileptons: largebins=""
jet_origin_fractions_stacked_dileptons:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 mumu_passjets x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 mumu_passjets y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 mumu_passjets z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 elmu_passjets x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 elmu_passjets y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 elmu_passjets z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `

jet_origin_fractions_stacked_control_el: largebins=""
jet_origin_fractions_stacked_control_el:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singleel_pretau x 1 eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singleel_pretau y 1 eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singleel_pretau z 1 eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `

jet_origin_fractions_stacked_control_mu: largebins=""
jet_origin_fractions_stacked_control_mu:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_pretau x 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_pretau y 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_pretau z 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `

jet_origin_fractions_stacked_control2_mu: largebins=""
jet_origin_fractions_stacked_control2_mu:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_control_pretau x 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_control_pretau y 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize f ${largebins} ${csv} 0.25 singlemu_control_pretau z 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `


# unstacked fractions of jets

jet_origin_fractions_unstacked_control_mu: largebins=""
jet_origin_fractions_unstacked_control_mu:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau x 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau y 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau z 1 mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `

jet_origin_fractions_unstacked_control_mu_mistau: largebins=""
jet_origin_fractions_unstacked_control_mu_mistau:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau x 1 mutau_pretau_mistau outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_mistau`
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau y 1 mutau_pretau_mistau outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_mistau`
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 singlemu_pretau z 1 mutau_pretau_mistau outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_mistau`

jet_origin_fractions_unstacked_dileptons: largebins=""
jet_origin_fractions_unstacked_dileptons:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 mumu_passjets x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 mumu_passjets y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 mumu_passjets z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 elmu_passjets x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 elmu_passjets y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 elmu_passjets z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `

jet_origin_fractions_unstacked_jets: largebins=""
jet_origin_fractions_unstacked_jets:
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTjet_qcd  x 1 QCD   outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTjet_qcd  y 1 QCD   outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTjet_qcd  z 1 QCD   outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTmu_wjets x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTmu_wjets y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked ${verbose} --normalize T ${largebins} ${csv} 0.25 HLTmu_wjets z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 



jet_fakerate_fractions_stacked_qcd:
	jet_fakerate_fractions_per_origin  0.25 HLTjet_qcd x 1  0   200    QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_fakerate_fractions_per_origin  0.25 HLTjet_qcd y 1 -2.5   2.5  QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_fakerate_fractions_per_origin  0.25 HLTjet_qcd z 1  0     1    QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 

jet_fakerate_fractions_stacked_wjets:
	jet_fakerate_fractions_per_origin  0.25 HLTmu_wjets x 1  0   200    WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_fakerate_fractions_per_origin  0.25 HLTmu_wjets y 1 -2.5   2.5  WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_fakerate_fractions_per_origin  0.25 HLTmu_wjets z 1  0     1    WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 


jet_fakerate_fractions_stacked_dileptons:
	jet_fakerate_fractions_per_origin  0.25 mumu_passjets x 1  0   200    mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 mumu_passjets y 1 -2.5   2.5  mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 mumu_passjets z 1  0     1    mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 elmu_passjets x 1  0   200    elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 elmu_passjets y 1 -2.5   2.5  elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 elmu_passjets z 1  0     1    elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 

jet_fakerate_fractions_stacked_control:
	jet_fakerate_fractions_per_origin  0.25 singleel_pretau x 1  0   200    eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 singleel_pretau y 1 -2.5   2.5  eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 singleel_pretau z 1  0     1    eltau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 singlemu_pretau x 1  0   200    mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 singlemu_pretau y 1 -2.5   2.5  mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_fakerate_fractions_per_origin  0.25 singlemu_pretau z 1  0     1    mutau_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 





jet_tau_origin_distrs_overlayed:
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd_tau x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd_tau y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd_tau z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets_tau x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets_tau y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets_tau z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 

jet_origin_distrs_overlayed:
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTjet_qcd z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 HLTmu_wjets z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 

jet_origin_distrs_dilep_overlayed:
	jet_distrs_per_origin_overlayed --normalize 0.25 mumu_passjets x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 mumu_passjets y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 mumu_passjets z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 elmu_passjets x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 elmu_passjets y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 elmu_passjets z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 

jet_origin_distrs_control_overlayed:
	jet_distrs_per_origin_overlayed --normalize 0.25 singleel_pretau x 1 singleel_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 singleel_pretau y 1 singleel_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 singleel_pretau z 1 singleel_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 singlemu_pretau x 1 singlemu_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 singlemu_pretau y 1 singlemu_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 
	jet_distrs_per_origin_overlayed --normalize 0.25 singlemu_pretau z 1 singlemu_pretau outdir/"${job}"/merged-sets/ `cat dtags_"${job}" ` 

jet_mc_fakerates:
	jet_fakerate_per_mc 0.25 HLTjet_qcd  z 1  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTjet_qcd  y 1 -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTjet_qcd  x 1  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_mc 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_mc 0.25 HLTmu_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_mc_fakerates_dilep:
	jet_fakerate_per_mc 0.25 mumu_passjets  z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 mumu_passjets  y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 mumu_passjets  x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets z 1  0   0.4  elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets y 1 -2.5 2.5  elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets x 1  0 150    elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

# elmu FAKERATES
jet_mc_fakerates_dilep_lumi:
	jet_fakerate_per_mc ${lumi} elmu_passjets z 1  0   0.4  elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc ${lumi} elmu_passjets y 1 -2.5 2.5  elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc ${lumi} elmu_passjets x 1  0 150    elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_qcd:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_qcd z 1  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_qcd x 1  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`

jet_origin_fakerates_wjets:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTmu_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_origin_fakerates_wjets_HLTjet:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 HLTjet_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_origin_fakerates: set_logy=T
jet_origin_fakerates: large=f
jet_origin_fakerates: jet_origin_fakerates_qcd jet_origin_fakerates_wjets

jet_origin_fakerates_elmu:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_mumu:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_elel:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elel_passjets z 1  0   0.4 elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elel_passjets y 1 -2.5 2.5 elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 elel_passjets x 1  0 150   elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_control_mu:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_control_mu_notau: set_logy=f
jet_origin_fakerates_control_mu_notau:
	jet_fakerate_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 s-mu_notau outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 s-mu_notau outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   s-mu_notau outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_control_mu_onlyttlj:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`

jet_origin_fakerates_control_mu_noltau:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 s-mu_noltau outdir/"${job}"/merged-sets/ `cat dtags_${job}_noltau`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 s-mu_noltau outdir/"${job}"/merged-sets/ `cat dtags_${job}_noltau`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   s-mu_noltau outdir/"${job}"/merged-sets/ `cat dtags_${job}_noltau`

jet_origin_fakerates_control_mu_mistau:
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_per_origin f ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`

# distrs of fakes
jet_origin_fakes_jets:
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTmu_wjets z 1  0   0.4 wjets outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTmu_wjets y 1 -2.5 2.5 wjets outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTmu_wjets x 1  0 150   wjets outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTjet_qcd  z 1  0   0.4 qcd   outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTjet_qcd  y 1 -2.5 2.5 qcd   outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 HLTjet_qcd  x 1  0 150   qcd   outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakes_dilep:
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakes_control_mu:
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau z 1  0   0.4 smu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau y 1 -2.5 2.5 smu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakes_per_origin T ${large} ${set_logy} 0.25 singlemu_pretau x 1  0 150   smu outdir/"${job}"/merged-sets/ `cat dtags_${job}`





jet_origin_fakerates_dileptons: set_logy=T
jet_origin_fakerates_dileptons: large=f
jet_origin_fakerates_dileptons: jet_origin_fakerates_elmu jet_origin_fakerates_mumu jet_origin_fakerates_elel

jet_fakerate_separate_origins:
	jet_fakerate_separate_origins 0.25 HLTjet_qcd z 1  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTjet_qcd x 1  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_fakerate_separate_origins_overlayed: max_y=0.005
jet_fakerate_separate_origins_overlayed: set_logy=T
jet_fakerate_separate_origins_overlayed:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTjet_qcd z 1  ${max_y}  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTjet_qcd y 1  ${max_y} -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTjet_qcd x 1  ${max_y}  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTmu_wjets z 1 ${max_y}  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTmu_wjets y 1 ${max_y} -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 HLTmu_wjets x 1 ${max_y}  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	#jet_fakerate_separate_origins_overlayed.cc

jet_fakerate_separate_origins_dilep:
	jet_fakerate_separate_origins 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_dilep: max_y=0.005
jet_fakerate_separate_origins_overlayed_dilep: set_logy=T
jet_fakerate_separate_origins_overlayed_dilep:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 mumu_passjets z 1 ${max_y}  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 mumu_passjets y 1 ${max_y} -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 mumu_passjets x 1 ${max_y}  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 elmu_passjets z 1 ${max_y}  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 elmu_passjets y 1 ${max_y} -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 elmu_passjets x 1 ${max_y}  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_control_passjet: max_y=0.005
jet_fakerate_separate_origins_overlayed_control_passjet: set_logy=T
jet_fakerate_separate_origins_overlayed_control_passjet:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet z 1 ${max_y}  0   0.4 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet y 1 ${max_y} -2.5 2.5 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet x 1 ${max_y}  0 150   s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet z 1 ${max_y}  0   0.4 s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet y 1 ${max_y} -2.5 2.5 s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet x 1 ${max_y}  0 150   s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_control_passjet_onlyttlj: set_logy=T
jet_fakerate_separate_origins_overlayed_control_passjet_onlyttlj:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet z 1 0.005  0   0.4 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet y 1 0.005 -2.5 2.5 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_passjet x 1 0.005  0 150   s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet z 1 0.005  0   0.4 s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet y 1 0.005 -2.5 2.5 s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_passjet x 1 0.005  0 150   s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`

jet_fakerate_separate_origins_overlayed_control: max_y=0.005
jet_fakerate_separate_origins_overlayed_control: set_logy=T
jet_fakerate_separate_origins_overlayed_control:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau z 1 ${max_y}  0   0.4 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau y 1 ${max_y} -2.5 2.5 s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau x 1 ${max_y}  0 150   s-mu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau z 1 ${max_y}  0   0.4 s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau y 1 ${max_y} -2.5 2.5 s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau x 1 ${max_y}  0 150   s-el outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_control_onlyttlj:
	jet_fakerate_separate_origins_overlayed T 0.25 singlemu_control_pretau z 1 0.005  0   0.4 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed T 0.25 singlemu_control_pretau y 1 0.005 -2.5 2.5 s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed T 0.25 singlemu_control_pretau x 1 0.005  0 150   s-mu_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed T 0.25 singleel_control_pretau z 1 0.005  0   0.4 s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed T 0.25 singleel_control_pretau y 1 0.005 -2.5 2.5 s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`
	jet_fakerate_separate_origins_overlayed T 0.25 singleel_control_pretau x 1 0.005  0 150   s-el_onlyttlj outdir/"${job}"/merged-sets/ `cat dtags_${job}_onlyttlj`

jet_fakerate_separate_origins_overlayed_control_mistau: set_logy="f"
jet_fakerate_separate_origins_overlayed_control_mistau:
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau z 1 0.005  0   0.4 s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau y 1 0.005 -2.5 2.5 s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singlemu_control_pretau x 1 0.005  0 150   s-mu_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau z 1 0.005  0   0.4 s-el_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau y 1 0.005 -2.5 2.5 s-el_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`
	jet_fakerate_separate_origins_overlayed ${set_logy} 0.25 singleel_control_pretau x 1 0.005  0 150   s-el_mistau outdir/"${job}"/merged-sets/ `cat dtags_${job}_mistau`

jet_fakefactor_jets:
	jet_fakefactor_data 1 HLTjet_qcd x 1  0 100 qcd outdir/${job}/merged-sets `cat dtags_${job}_JetHT`
	jet_fakefactor_data 1 HLTjet_qcd x 1  0 500 qcd outdir/${job}/merged-sets `cat dtags_${job}_JetHT`
	jet_fakefactor_data 1 HLTjet_qcd y 1 -3 3   qcd outdir/${job}/merged-sets `cat dtags_${job}_JetHT`
	jet_fakefactor_data 1 HLTjet_qcd z 1  0 1   qcd outdir/${job}/merged-sets `cat dtags_${job}_JetHT`
	#
	jet_fakefactor_data 1 HLTmu_wjets x 1  0 100 wjets outdir/${job}/merged-sets `cat dtags_${job}_SingleMuon`
	jet_fakefactor_data 1 HLTmu_wjets x 1  0 500 wjets outdir/${job}/merged-sets `cat dtags_${job}_SingleMuon`
	jet_fakefactor_data 1 HLTmu_wjets y 1 -3 3   wjets outdir/${job}/merged-sets `cat dtags_${job}_SingleMuon`
	jet_fakefactor_data 1 HLTmu_wjets z 1  0 1   wjets outdir/${job}/merged-sets `cat dtags_${job}_SingleMuon`

jet_fakefactor_dileptons:
	jet_fakefactor_data 1 elel_passjets x 1  0 100 elel outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elel_passjets x 1  0 500 elel outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elel_passjets y 1 -3 3   elel outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elel_passjets z 1  0 1   elel outdir/${job}/merged-sets `cat dtags_${job}`
	#
	jet_fakefactor_data 1 elmu_passjets x 1  0 100 elmu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elmu_passjets x 1  0 500 elmu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elmu_passjets y 1 -3 3   elmu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 elmu_passjets z 1  0 1   elmu outdir/${job}/merged-sets `cat dtags_${job}`
	#
	jet_fakefactor_data 1 mumu_passjets x 1  0 100 mumu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 mumu_passjets x 1  0 500 mumu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 mumu_passjets y 1 -3 3   mumu outdir/${job}/merged-sets `cat dtags_${job}`
	jet_fakefactor_data 1 mumu_passjets z 1  0 1   mumu outdir/${job}/merged-sets `cat dtags_${job}`

overlay_std_fakerates:
	overlay_th1d_distrs Y ${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3
	overlay_th1d_distrs Y ${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3
	overlay_th1d_distrs Y ${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3

overlay_fakerates:
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3

overlay_fakerates_4:
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_${largebins}_fake-rates-per-x-pt_overlayed.png     qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr${largebins}_over_HLTjet_qcd_jets_distr${largebins}_x.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr${largebins}_over_HLTmu_wjets_jets_distr${largebins}_x.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr${largebins}_over_mumu_passjets_jets_distr${largebins}_x.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr${largebins}_over_elmu_passjets_jets_distr${largebins}_x.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_${largebins}_fake-rates-per-y-eta_overlayed.png    qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr${largebins}_over_HLTjet_qcd_jets_distr${largebins}_y.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr${largebins}_over_HLTmu_wjets_jets_distr${largebins}_y.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr${largebins}_over_mumu_passjets_jets_distr${largebins}_y.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr${largebins}_over_elmu_passjets_jets_distr${largebins}_y.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_${largebins}_fake-rates-per-z-radius_overlayed.png qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr${largebins}_over_HLTjet_qcd_jets_distr${largebins}_z.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr${largebins}_over_HLTmu_wjets_jets_distr${largebins}_z.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr${largebins}_over_mumu_passjets_jets_distr${largebins}_z.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr${largebins}_over_elmu_passjets_jets_distr${largebins}_z.root h3

overlay_fakerates_qcd_n_mumu:
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_jets_distr_x.root h3  mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3 
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_jets_distr_y.root h3  mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3 
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_jets_distr_z.root h3  mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3 

overlay_fakerates_sleps_control:
	overlay_th1d_distrs Y ${job}_fake-rates-per-x-pt_overlayed.png     ${job}-el outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singleel_pretau_tau_jets_distr_over_singleel_pretau_jets_distr_x.root h3 ${job}-mu outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singlemu_pretau_tau_jets_distr_over_singlemu_pretau_jets_distr_x.root h3
	overlay_th1d_distrs Y ${job}_fake-rates-per-y-eta_overlayed.png     ${job}-el outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singleel_pretau_tau_jets_distr_over_singleel_pretau_jets_distr_y.root h3 ${job}-mu outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singlemu_pretau_tau_jets_distr_over_singlemu_pretau_jets_distr_y.root h3
	overlay_th1d_distrs Y ${job}_fake-rates-per-z-radius_overlayed.png     ${job}-el outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singleel_pretau_tau_jets_distr_over_singleel_pretau_jets_distr_z.root h3 ${job}-mu outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_singlemu_pretau_tau_jets_distr_over_singlemu_pretau_jets_distr_z.root h3

overlay_dilep_fakerates:
	overlay_th1d_distrs Y ${ttjob}_dileps_fake-rates-per-x-pt_overlayed.png  elel-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_x.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs Y ${ttjob}_dileps_fake-rates-per-y-eta_overlayed.png  elel-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_y.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs Y ${ttjob}_dileps_fake-rates-per-z-redius_overlayed.png  elel-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_z.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_z.root h3




overlay_fakefactors_jets:
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-100.png  qcd outdir/${job}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_x_qcd_0-100.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${job}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_x_wjets_0-100.root HLTmu_wjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-500.png  qcd outdir/${job}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_x_qcd_0-500.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${job}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_x_wjets_0-500.root HLTmu_wjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_y.png  qcd outdir/${job}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_y_qcd_-3-3.root HLTjet_qcd_tau_jets_distr_y wjets outdir/${job}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_y_wjets_-3-3.root HLTmu_wjets_tau_jets_distr_y
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_z.png  qcd outdir/${job}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_z_qcd_0-1.root HLTjet_qcd_tau_jets_distr_z wjets outdir/${job}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_z_wjets_0-1.root HLTmu_wjets_tau_jets_distr_z

overlay_fakefactors_dileptons:
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-100.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_x_elel_0-100.root elel_passjets_tau_jets_distr_x elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-100.root elmu_passjets_tau_jets_distr_x mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-100.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-500.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_x_elel_0-500.root elel_passjets_tau_jets_distr_x elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-500.root elmu_passjets_tau_jets_distr_x mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-500.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_y_-3-3.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_y_elel_-3-3.root elel_passjets_tau_jets_distr_y elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_y_elmu_-3-3.root elmu_passjets_tau_jets_distr_y mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_y_mumu_-3-3.root mumu_passjets_tau_jets_distr_y
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_z_0-1.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_z_elel_0-1.root elel_passjets_tau_jets_distr_z elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_z_elmu_0-1.root elmu_passjets_tau_jets_distr_z mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_z_mumu_0-1.root mumu_passjets_tau_jets_distr_z

overlay_fakefactors_jets_n_dileptons:
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_x_0-100.png qcd   outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_x_qcd_0-100.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_x_wjets_0-100.root HLTmu_wjets_tau_jets_distr_x elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-100.root elmu_passjets_tau_jets_distr_x mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-100.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_x_0-500.png qcd   outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_x_qcd_0-500.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_x_wjets_0-500.root HLTmu_wjets_tau_jets_distr_x elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-500.root elmu_passjets_tau_jets_distr_x mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-500.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_y_-3-3.png qcd   outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_y_qcd_-3-3.root HLTjet_qcd_tau_jets_distr_y wjets outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_y_wjets_-3-3.root HLTmu_wjets_tau_jets_distr_y elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_y_elmu_-3-3.root elmu_passjets_tau_jets_distr_y mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_y_mumu_-3-3.root mumu_passjets_tau_jets_distr_y
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_z_0-1.png qcd   outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTjet_qcd_z_qcd_0-1.root HLTjet_qcd_tau_jets_distr_z wjets outdir/${jjob}/merged-sets/jobsums/FakeFactor_HLTmu_wjets_z_wjets_0-1.root HLTmu_wjets_tau_jets_distr_z elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_z_elmu_0-1.root elmu_passjets_tau_jets_distr_z mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_z_mumu_0-1.root mumu_passjets_tau_jets_distr_z




overlay_fakerates_test:
	overlay_th1d_distrs Y test_fake-rates-per-x-pt_overlayed.png     qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3
	overlay_th1d_distrs Y test_fake-rates-per-y-eta_overlayed.png    qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3
	overlay_th1d_distrs Y test_fake-rates-per-z-radius_overlayed.png qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3



ntuple_distrs_el:
	ntuple_distr --verbose el -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 " ${range} sel-el_sel1lep3jets_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose el -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && njets > 2 " ${range} sel-el_sel1lep3jets_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose el -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 " ${range} sel-el_sel1lep3jets40met_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose el -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 " ${range} sel-el_sel1lep3jets40met1b_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose el 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 " ${range} sel-el_sel1lep3jets40met1b1t_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose el 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && tau0_id * lep0_id < 0 " ${range} sel-el_sel1lep3jets40met1b1tOS_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`

ntuple_distrs_mu:
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 " ${range} sel-mu_sel1lep3jets_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 " ${range} sel-mu_sel1lep3jets_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 " ${range} sel-mu_sel1lep3jets40met_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 " ${range} sel-mu_sel1lep3jets40met1b_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 " ${range} sel-mu_sel1lep3jets40met1b1t_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && tau0_id * lep0_id < 0 " ${range} sel-mu_sel1lep3jets40met1b1tOS_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`


ntuple_distrs_mu_from3jets:
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 " ${range} sel-mu_sel1lep3jets_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 " ${range} sel-mu_sel1lep3jets40met_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu -1 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 " ${range} sel-mu_sel1lep3jets40met1b_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 " ${range} sel-mu_sel1lep3jets40met1b1t_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`
	ntuple_distr --verbose mu 0.97 ${bSF} T T T f f ${lumi_bcdef} ${lumi_gh} ${distr} "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && tau0_id * lep0_id < 0 " ${range} sel-mu_sel1lep3jets40met1b1tOS_${file_suffix} ${distr_X} outdir/${job}/merged-sets/ `cat dtags_${job}`

ntuple_likelihood_mu_ouLJ3j: bSF=T
ntuple_likelihood_mu_ouLJ3j: scale=1.036
ntuple_likelihood_mu_ouLJ3j: sys_shift=""
ntuple_likelihood_mu_ouLJ3j: tau_ES=
ntuple_likelihood_mu_ouLJ3j:
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets ==3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-3j-Scale-${scale}_WaMC_${sys_shift} " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`

ntuple_likelihood_mu_ouLJ4j: bSF=T
ntuple_likelihood_mu_ouLJ4j: scale=1.036
ntuple_likelihood_mu_ouLJ4j: sys_shift=""
ntuple_likelihood_mu_ouLJ4j: tau_ES=
ntuple_likelihood_mu_ouLJ4j:
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-4j-Scale-${scale}_WaMC_${sys_shift} " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`

ntuple_likelihood_mu: bSF=T
ntuple_likelihood_mu: scale=1.036
ntuple_likelihood_mu: sys_shift=""
ntuple_likelihood_mu: tau_ES=
ntuple_likelihood_mu:
	# in-lj
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance > 0 && lj_peak_distance < 1600)  ${tau_ES}" "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_inLJ-Scale-${scale}_WaMC " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`
	# out 3j
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets ==3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-3j-Scale-${scale}_WaMC " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`
	# out 4j
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-4j-Scale-${scale}_WaMC " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`
	# pre-tau
	# in-lj
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance > 0 && lj_peak_distance < 1600)  ${tau_ES}" "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_inLJ-Scale-${scale}_WaMC_pretau " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`
	# out 3j
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets ==3 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-3j-Scale-${scale}_WaMC_pretau " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`
	# out 4j
	ntuple_distr --verbose mu ${scale} -1 ${bSF} T T T ${sys_shift} f f 19094 16079 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 3 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-4j-Scale-${scale}_WaMC_pretau " #mu MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_mu_WaMC`

ntuple_likelihood_el: bSF=T
ntuple_likelihood_el: scale=0.983
ntuple_likelihood_el: sys_shift=""
ntuple_likelihood_el: tau_ES=
ntuple_likelihood_el:
	# in-lj
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance > 0 && lj_peak_distance < 1600) ${tau_ES}" "(25,0,250)"  likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_inLJ-Scale-${scale}_WaMC     " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`
	# out 3j
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets ==3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES}" "(25,0,250)" likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-3j-Scale-${scale}_WaMC " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`
	# out 4j
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 3 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES}" "(25,0,250)" likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-4j-Scale-${scale}_WaMC " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`
	# pre-tau
	# in-lj
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance > 0 && lj_peak_distance < 1600)  ${tau_ES}" "(25,0,250)"  likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_inLJ-Scale-${scale}_WaMC_pretau     " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`
	# out 3j
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets ==3 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-3j-Scale-${scale}_WaMC_pretau " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`
	# out 4j
	ntuple_distr --verbose el ${scale} -1 ${bSF} T T T ${sys_shift} f f 19102 12210 "sqrt(2*lep0_p4.pt()*met_corrected.pt() * (1 - cos(lep0_p4.phi() - met_corrected.phi())))" "HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 3 && met_corrected.pt() > 40 && nbjets > 0 && (lj_peak_distance < 0 || lj_peak_distance >= 1600) ${tau_ES} " "(25,0,250)" likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-4j-Scale-${scale}_WaMC_pretau " #el MET M_{T} " " Events " /eos/user/o/otoldaie/ttbar-leptons-80X_data/v11.8_merged-sets/ `cat dtags_v11.8_el_WaMC`

test_ntuple_channels:
	python channel_distrs.py outdir/v12.2/merged-sets/ el --lumi 28000
	python channel_distrs.py outdir/v12.2/merged-sets/ mu --lumi 30800
	python channel_distrs.py outdir/v12.2/merged-sets/ elmu_hltel --lumi 28000
	python channel_distrs.py outdir/v12.2/merged-sets/ elmu_hltmu --lumi 30800
	python channel_distrs.py outdir/v12.2/merged-sets/ elmu_hltelmu --lumi 30800
	python channel_distrs.py outdir/v12.2/merged-sets/ taumutauh --lumi 30800
	python channel_distrs.py outdir/v12.2/merged-sets/ mujets --lumi 30800


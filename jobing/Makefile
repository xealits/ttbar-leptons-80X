# maybe I just should make a bunch of scripts with these
# just..lets keep it together here for now

#just in case
SHELL:=/bin/bash

reduce:
	./jobing/reduce.zsh ${job} jobing/dtags/${dtags}

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
	ls outdir/"${job}"/*root | sed -e 's_^.*'${job}'/__' -e 's/_[0-9][0-9]*.root//' | sort | uniq > dtags_job_${job}

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
	for h in `cat cur_th1d_histos_pt_eta_energy`;     \
	do                                                \
	NEWstacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other`;             \
	do                                                \
	NEWstacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singleel: normalize=f
selection_control_histos_singleel: logy=f
selection_control_histos_singleel:
	for h in `cat cur_th1d_histos_pt_eta_energy_singleel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_singleel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_elel: normalize=f
selection_control_histos_elel: logy=f
selection_control_histos_elel:
	for h in `cat cur_th1d_histos_pt_eta_energy_elel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_elel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_elmu: normalize=f
selection_control_histos_elmu: logy=f
selection_control_histos_elmu:
	for h in `cat cur_th1d_histos_pt_eta_energy_elmu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_elmu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_mumu: normalize=f
selection_control_histos_mumu: logy=f
selection_control_histos_mumu:
	for h in `cat cur_th1d_histos_pt_eta_energy_mumu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_mumu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singleel_other: normalize=f
selection_control_histos_singleel_other: logy=f
selection_control_histos_singleel_other:
	for h in `cat cur_th1d_histos_other_singleel`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singlemu: normalize=f
selection_control_histos_singlemu: logy=f
selection_control_histos_singlemu:
	for h in `cat cur_th1d_histos_pt_eta_energy_singlemu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other_singlemu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_singlemu_other: normalize=f
selection_control_histos_singlemu_other: logy=f
selection_control_histos_singlemu_other:
	for h in `cat cur_th1d_histos_other_singlemu`;             \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta: normalize=f
selection_control_histos_eta: logy=f
selection_control_histos_eta:
	for h in `cat cur_th1d_histos_eta`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta_singleel: normalize=f
selection_control_histos_eta_singleel: logy=f
selection_control_histos_eta_singleel:
	for h in `cat cur_th1d_histos_eta_singleel`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_eta_singlemu: normalize=f
selection_control_histos_eta_singlemu: logy=f
selection_control_histos_eta_singlemu:
	for h in `cat cur_th1d_histos_eta_singlemu`;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 5 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_pileup: normalize=f
selection_control_histos_pileup: logy=f
selection_control_histos_pileup:
	for h in pileup_passtrig_nGoodPV_rawWeight pileup_passtrig_nGoodPV_weight pileup_singlelepton_nGoodPV_rawWeight pileup_singlelepton_nGoodPV_weight pileup_dilepton_nGoodPV_rawWeight pileup_dilepton_nGoodPV_weight pileup_sel_passos_nGoodPV_rawWeight pileup_sel_passos_nGoodPV_weight pileup_smu_passos_nGoodPV_rawWeight pileup_smu_passos_nGoodPV_weight ;     \
	do                                                \
	stacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 0 50 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

selection_control_histos_rho: normalize=f
selection_control_histos_rho: logy=f
selection_control_histos_rho:
	for h in rho_beforetrig_rawWeight rho_beforetrig_weight rho_passtrig_rawWeight rho_passtrig_weight rho_singlelepton_rawWeight rho_singlelepton_weight  rho_dilepton_rawWeight rho_dilepton_weight \
		 rhoCentral_beforetrig_rawWeight rhoCentral_beforetrig_weight rhoCentral_passtrig_rawWeight rhoCentral_passtrig_weight rhoCentral_singlelepton_rawWeight rhoCentral_singlelepton_weight  rhoCentral_dilepton_rawWeight rhoCentral_dilepton_weight \
		 rhoCentralNeutral_beforetrig_rawWeight rhoCentralNeutral_beforetrig_weight rhoCentralNeutral_passtrig_rawWeight rhoCentralNeutral_passtrig_weight rhoCentralNeutral_singlelepton_rawWeight rhoCentralNeutral_singlelepton_weight  rhoCentralNeutral_dilepton_rawWeight rhoCentralNeutral_dilepton_weight \
		 rhoCentralChargedPileUp_beforetrig_rawWeight rhoCentralChargedPileUp_beforetrig_weight rhoCentralChargedPileUp_passtrig_rawWeight rhoCentralChargedPileUp_passtrig_weight rhoCentralChargedPileUp_singlelepton_rawWeight rhoCentralChargedPileUp_singlelepton_weight  rhoCentralChargedPileUp_dilepton_rawWeight rhoCentralChargedPileUp_dilepton_weight ;     \
	do                                                \
	NEWstacked_histo_distr_with_ratio ${normalize} ${logy} ${lumi} $$h $$h 1 -1 -1 outdir/"${job}"/merged-sets/ `cat dtags_"${job}"`; \
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
fakerate_ratio_plots:
	test -s ${file}
	# make stops if a command in the recipe fails...
	histo_project3toN_ratio_distr x pT     "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1  0 150   ${file}
	histo_project3toN_ratio_distr y eta    "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1 -2.5 2.5 ${file}
	histo_project3toN_ratio_distr z radius "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1  0   0.3 ${file}

.PHONY:jet_fakerate_plots_JetHT
jet_fakerate_plots_JetHT: file=outdir/${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
jet_fakerate_plots_JetHT: selection=HLTjet_qcd
jet_fakerate_plots_JetHT: fakerate_ratio_plots

.PHONY:jet_fakerate_plots_SingleMuon
jet_fakerate_plots_SingleMuon: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
jet_fakerate_plots_SingleMuon: selection=HLTmu_wjets
jet_fakerate_plots_SingleMuon: fakerate_ratio_plots

jet_fakerate_plots: jet_fakerate_plots_JetHT jet_fakerate_plots_SingleMu

.PHONY:jet_fakerate_plots_control_el
jet_fakerate_plots_control_el: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
jet_fakerate_plots_control_el: selection=singleel_pretau
jet_fakerate_plots_control_el: fakerate_ratio_plots

.PHONY:jet_fakerate_plots_control_mu
jet_fakerate_plots_control_mu: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
jet_fakerate_plots_control_mu: selection=singlemu_pretau
jet_fakerate_plots_control_mu: fakerate_ratio_plots

jet_fakerate_plots_control: jet_fakerate_plots_control_el jet_fakerate_plots_control_mu

mumu_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
mumu_fakerate_ratio_plots: selection=mumu_passjets
mumu_fakerate_ratio_plots: fakerate_ratio_plots

elmu_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elmu_fakerate_ratio_plots: selection=elmu_passjets
elmu_fakerate_ratio_plots: fakerate_ratio_plots

elel_fakerate_ratio_plots: file=outdir/${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elel_fakerate_ratio_plots: selection=elel_passjets
elel_fakerate_ratio_plots: fakerate_ratio_plots

dilep_fakerate_ratio_plots: elel_fakerate_ratio_plots mumu_fakerate_ratio_plots elmu_fakerate_ratio_plots

fakerate_study_extract_all:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} HLTjet_qcd_jet_partonFlavour HLTmu_wjets_jet_partonFlavour elmu_passjets_jet_partonFlavour mumu_passjets_jet_partonFlavour

jet_mc_distrs:
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

tau_mc_distrs:
	stacked_histo_project2to1_distr --normalize 1.7 HLTjet_qcd_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project2to1_distr --normalize 1.7 HLTjet_qcd_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project2to1_distr --normalize 1.7 HLTmu_wjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project2to1_distr --normalize 1.7 HLTmu_wjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`


jets_taus_mc_distrs_dileptons:
	stacked_histo_project3to1_distr --normalize 1.7 mumu_passjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 mumu_passjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 mumu_passjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 elmu_passjets_jets_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 elmu_passjets_jets_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project3to1_distr --normalize 1.7 elmu_passjets_jets_distr z 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	#
	stacked_histo_project2to1_distr --normalize 1.7 mumu_passjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 mumu_passjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 elmu_passjets_taus_distr x 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	stacked_histo_project2to1_distr --normalize 1.7 elmu_passjets_taus_distr y 1 outdir/"${job}"/merged-sets/ `cat dtags_${job}`



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



jet_origin_fractions_stacked_qcd:
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTjet_qcd x 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTjet_qcd y 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTjet_qcd z 1 QCD outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_JetHT ` 

jet_origin_fractions_stacked_wjets:
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTmu_wjets x 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTmu_wjets y 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 
	jet_distrs_origin_fractions_stacked --normalize 0.25 HLTmu_wjets z 1 WJets outdir/"${job}"/merged-sets/ `cat dtags_"${job}"_SingleMuon ` 


jet_origin_fractions_stacked_dileptons:
	jet_distrs_origin_fractions_stacked --normalize 0.25 mumu_passjets x 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked --normalize 0.25 mumu_passjets y 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked --normalize 0.25 mumu_passjets z 1 mumu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked --normalize 0.25 elmu_passjets x 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked --normalize 0.25 elmu_passjets y 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `
	jet_distrs_origin_fractions_stacked --normalize 0.25 elmu_passjets z 1 elmu outdir/"${job}"/merged-sets/ `cat dtags_"${job}" `




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

jet_origin_fakerates_qcd:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_qcd z 1  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_qcd x 1  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`

jet_origin_fakerates_wjets:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTmu_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_origin_fakerates_wjets_HLTjet:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 HLTjet_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_origin_fakerates: set_logy=T
jet_origin_fakerates: large=f
jet_origin_fakerates: jet_origin_fakerates_qcd jet_origin_fakerates_wjets

jet_origin_fakerates_elmu:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_mumu:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_elel:
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elel_passjets z 1  0   0.4 elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elel_passjets y 1 -2.5 2.5 elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin ${large} ${set_logy} 0.25 elel_passjets x 1  0 150   elel outdir/"${job}"/merged-sets/ `cat dtags_${job}`

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

jet_fakerate_separate_origins_overlayed:
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd z 1  0   0.4 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd x 1  0 150   JetHT outdir/"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets x 1  0 150   SingleMuon outdir/"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	#jet_fakerate_separate_origins_overlayed.cc

jet_fakerate_separate_origins_dilep:
	jet_fakerate_separate_origins 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_dilep:
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets z 1  0   0.4 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets y 1 -2.5 2.5 mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets x 1  0 150   mumu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets z 1  0   0.4 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets y 1 -2.5 2.5 elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets x 1  0 150   elmu outdir/"${job}"/merged-sets/ `cat dtags_${job}`

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
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs Y ${ttjob}_${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel outdir/${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 mumu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3 elmu-sel outdir/${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_z.root h3

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
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-100.png  qcd outdir/${job}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_x_qcd_0-100.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${job}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_x_wjets_0-100.root HLTmu_wjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-500.png  qcd outdir/${job}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_x_qcd_0-500.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${job}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_x_wjets_0-500.root HLTmu_wjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_y.png  qcd outdir/${job}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_y_qcd_-3-3.root HLTjet_qcd_tau_jets_distr_y wjets outdir/${job}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_y_wjets_-3-3.root HLTmu_wjets_tau_jets_distr_y
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_z.png  qcd outdir/${job}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_z_qcd_0-1.root HLTjet_qcd_tau_jets_distr_z wjets outdir/${job}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_z_wjets_0-1.root HLTmu_wjets_tau_jets_distr_z

overlay_fakefactors_dileptons:
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-100.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_x_elel_0-100.root elel_passjets_tau_jets_distr_x elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-100.root elmu_passjets_tau_jets_distr_x mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-100.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_x_0-500.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_x_elel_0-500.root elel_passjets_tau_jets_distr_x elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-500.root elmu_passjets_tau_jets_distr_x mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-500.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_y_-3-3.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_y_elel_-3-3.root elel_passjets_tau_jets_distr_y elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_y_elmu_-3-3.root elmu_passjets_tau_jets_distr_y mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_y_mumu_-3-3.root mumu_passjets_tau_jets_distr_y
	overlay_th1d_distrs f ${job}_fakefactors_overlayed_z_0-1.png  elel outdir/${job}/merged-sets/jobsums/FakeFactor_elel_passjets_z_elel_0-1.root elel_passjets_tau_jets_distr_z elmu outdir/${job}/merged-sets/jobsums/FakeFactor_elmu_passjets_z_elmu_0-1.root elmu_passjets_tau_jets_distr_z mumu outdir/${job}/merged-sets/jobsums/FakeFactor_mumu_passjets_z_mumu_0-1.root mumu_passjets_tau_jets_distr_z

overlay_fakefactors_jets_n_dileptons:
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_x_0-100.png qcd   outdir/${jjob}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_x_qcd_0-100.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${jjob}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_x_wjets_0-100.root HLTmu_wjets_tau_jets_distr_x elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-100.root elmu_passjets_tau_jets_distr_x mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-100.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_x_0-500.png qcd   outdir/${jjob}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_x_qcd_0-500.root HLTjet_qcd_tau_jets_distr_x wjets outdir/${jjob}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_x_wjets_0-500.root HLTmu_wjets_tau_jets_distr_x elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_x_elmu_0-500.root elmu_passjets_tau_jets_distr_x mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_x_mumu_0-500.root mumu_passjets_tau_jets_distr_x
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_y_-3-3.png qcd   outdir/${jjob}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_y_qcd_-3-3.root HLTjet_qcd_tau_jets_distr_y wjets outdir/${jjob}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_y_wjets_-3-3.root HLTmu_wjets_tau_jets_distr_y elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_y_elmu_-3-3.root elmu_passjets_tau_jets_distr_y mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_y_mumu_-3-3.root mumu_passjets_tau_jets_distr_y
	overlay_th1d_distrs f ${jjob}_${vjob}_fakefactors_overlayed_z_0-1.png qcd   outdir/${jjob}/merged-sets/jobsums/HLTjet_qcd_FakeFactor_z_qcd_0-1.root HLTjet_qcd_tau_jets_distr_z wjets outdir/${jjob}/merged-sets/jobsums/HLTmu_wjets_FakeFactor_z_wjets_0-1.root HLTmu_wjets_tau_jets_distr_z elmu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_elmu_passjets_z_elmu_0-1.root elmu_passjets_tau_jets_distr_z mumu  outdir/${vjob}/merged-sets/jobsums/FakeFactor_mumu_passjets_z_mumu_0-1.root mumu_passjets_tau_jets_distr_z




overlay_fakerates_test:
	overlay_th1d_distrs Y test_fake-rates-per-x-pt_overlayed.png     qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3
	overlay_th1d_distrs Y test_fake-rates-per-y-eta_overlayed.png    qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3
	overlay_th1d_distrs Y test_fake-rates-per-z-radius_overlayed.png qcd-sel ${dir}/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 wjets-HLTjet ${dir}/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3





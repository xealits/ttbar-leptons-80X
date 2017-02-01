# maybe I just should make a bunch of scripts with these
# just..lets keep it together here for now

#just in case
SHELL:=/bin/bash

reduce:
	./jobing/reduce.zsh ${job} jobing/dtags/${dtags}

extract_all: extract_eventflow extract_weightflow extract_jet_origins extract_fakerates extract_pileup
	# DONE

extract_eventflow:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} eventflow_mu eventflow_el eventflow_elmu eventflow_elel eventflow_mumu

extract_weightflow:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} weightflow_mu weightflow_el weightflow_elmu weightflow_elel weightflow_mumu

extract_jet_origins:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} smu_pretau_jet_origin_ids sel_pretau_jet_origin_ids

extract_all_jet_origins:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} smu_pretau_jet_origin_ids sel_pretau_jet_origin_ids \
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
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} pileup_passtrig_rawWeight_pergoodpv pileup_passtrig_weight_pergoodpv

get_dtags:
	ls test/tests/outdir_"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > dtags_${job}

get_job_dtags:
	ls test/tests/outdir_"${job}"/*root | sed -e 's_^.*'${job}'/__' -e 's/_[0-9][0-9]*.root//' | sort | uniq > dtags_job_${job}

selection_control_histos:
	for h in `cat cur_th1d_histos_pt_eta_energy`;     \
	do                                                \
	stacked_histo_distr ${lumi} $$h 4 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other`;             \
	do                                                \
	stacked_histo_distr ${lumi} $$h 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done

stacked_fakerate_histos:
	for h in HLTjetmu_wjets_jets_distr;  \
	do                                                           \
	stacked_histo_project3toX_distr ${lumi} $$h x 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	stacked_histo_project3toX_distr ${lumi} $$h y 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	stacked_histo_project3toX_distr ${lumi} $$h z 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	#jobing/stacked_histo_project3toX_distr.cc

done_ck:
	ls test/tests/outdir_"${job}"/*job_done | wc -l
	ls test/tests/outdir_"${job}"/*root | wc -l
	ls test/tests/outdir_"${job}"/*cfg.py | wc -l

.PHONY:missing_jobs
missing_jobs:
	ls test/tests/outdir_${job}/*cfg.py | sed 's/_cfg.py/.job_done/' > submitted_jobs
	for j in `cat submitted_jobs`;        \
	do                                    \
	if [ ! -e $$j ]; then echo $$j; fi;   \
	done | sed -e 's/job_done/sh/' -e "s_"${job}"_"${job}"//FARM/inputs_" > missing_jobs
	wc -l missing_jobs

resubmit: missing_jobs
	for j in `cat missing_jobs`;                                    \
	do                                                              \
	  grep -h $$j test/tests/outdir_"${job}"//FARM/inputs/*bsub.sh; \
	done > missing_jobs_bsubs
	wc -l missing_jobs*
	source missing_jobs_bsubs

.PHONY:completeness
completeness:
	printf "dtag\tjob\tjson\troot\tcfg\n"
	printf "\n"
	for d in `echo ${dtags}`; \
	do                        \
	printf "%s\t%d\t%d\t%d\t%d\n" $$d `ls test/tests/outdir_${job}/"$$d"*job_done | wc -l` `ls test/tests/outdir_${job}/"$$d"*json | wc -l`  `ls test/tests/outdir_${job}/"$$d"*root | wc -l` `ls test/tests/outdir_${job}/"$$d"*cfg.py | wc -l`;  \
	done

#job=v9.15
#outdir=test/tests/outdir_${job}
#echo $outdir
#dtags="SingleElectron2016B SingleElectron2016C SingleElectron2016D SingleElectron2016E SingleElectron2016F SingleElectron2016G SingleElectron2016H"


proxy:
	# for consistency
	make-x509-proxy

calculate_luminosity:
	cat test/tests/outdir_"${job}"/*json > "${job}"_prejson
	sed -i 's/}{/,/g'     "${job}"_prejson
	sed -i 's/,,,*/,/g'   "${job}"_prejson
	sed -i 's/ ,"/,\n"/g' "${job}"_prejson
	sed -i -e 's/"[0-9]/[&/' -e 's/]]/]]]/' -e 's/{/[/' -e 's/}/]/' -e 's/:/,/' "${job}"_prejson
	#cd analysis/job-lumis
	./merge_json.py "${job}"_prejson "${job}".json
	brilcalc lumi -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i "${job}".json | tee "${job}".brilcalc.out
	mv "${job}"_prejson analysis/job-lumis
	mv "${job}".json analysis/job-lumis
	mv "${job}".brilcalc.out analysis/job-lumis


jet_origins:
	# "HLTjet_qcd_jet_partonFlavour HLTmu_qcd_jet_partonFlavour HLTjetmu_qcd_jet_partonFlavour HLTjet_wjets_jet_partonFlavour HLTmu_wjets_jet_partonFlavour HLTjetmu_wjets_jet_partonFlavour"
	# "dilep_passjets_jet_partonFlavour dilep_passjetsNbtag_jet_partonFlavour elel_passjets_jet_partonFlavour"
	ls test/tests/outdir_"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > test_dtags
	for d in `cat test_dtags | sed -n '/^MC/p' `; do jet_origins test/tests/outdir_"${job}"/merged-sets/ $$d ${channels}; done

jet_weightflow:
	# "HLTjet_qcd_jet_partonFlavour HLTmu_qcd_jet_partonFlavour HLTjetmu_qcd_jet_partonFlavour HLTjet_wjets_jet_partonFlavour HLTmu_wjets_jet_partonFlavour HLTjetmu_wjets_jet_partonFlavour"
	# "dilep_passjets_jet_partonFlavour dilep_passjetsNbtag_jet_partonFlavour elel_passjets_jet_partonFlavour"
	ls test/tests/outdir_"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > test_dtags
	for d in `cat test_dtags | sed -n '/^MC/p' `; do jet_weightflow test/tests/outdir_"${job}"/merged-sets/ $$d; done

.PHONY:fakerate_distr_plots
fakerate_distr_plots:
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root -a -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/${file}
	test -s ${file}
	# make stops if a command in the recipe fails...
	# || { echo foo; exit 1 }
	histo_project3toN_distr x "pT (GeV)" ${content_name}    ${distr} 1  0 150   ${file}
	histo_project3toN_distr y eta        ${content_name}    ${distr} 1 -2.5 2.5 ${file}
	histo_project3toN_distr z radius     ${content_name}    ${distr} 1  0   0.3 ${file}
	histo_project3toN_distr_atslice x "pT (GeV)" ${content_name}    ${distr} 1  0 150   y 2 3 ${file}

jet_distr_plots_JetHT: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
jet_distr_plots_JetHT: distr=HLTjet_qcd_jets_distr
jet_distr_plots_JetHT: content_name="N jets"
jet_distr_plots_JetHT: fakerate_distr_plots

jet_distr_plots_SingleMu: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
jet_distr_plots_SingleMu: distr=HLTmu_wjets_jets_distr
jet_distr_plots_SingleMu: content_name="N jets"
jet_distr_plots_SingleMu: fakerate_distr_plots

taujet_distr_plots_JetHT: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
taujet_distr_plots_JetHT: distr=HLTjet_qcd_tau_jets_distr
taujet_distr_plots_JetHT: content_name="N tau jets"
taujet_distr_plots_JetHT: fakerate_distr_plots

taujet_distr_plots_SingleMu: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
taujet_distr_plots_SingleMu: distr=HLTmu_wjets_tau_jets_distr
taujet_distr_plots_SingleMu: content_name="N tau jets"
taujet_distr_plots_SingleMu: fakerate_distr_plots

.PHONY:fakerate_ratio_plots
fakerate_ratio_plots:
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root -a -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/${file}
	test -s ${file}
	# make stops if a command in the recipe fails...
	# || { echo foo; exit 1 }
	histo_project3toN_ratio_distr x pT     "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1  0 150   ${file}
	histo_project3toN_ratio_distr y eta    "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1 -2.5 2.5 ${file}
	histo_project3toN_ratio_distr z radius "fake rate" ${selection}_tau_jets_distr ${selection}_jets_distr 1  0   0.3 ${file}

.PHONY:jet_fakerate_plots_JetHT
jet_fakerate_plots_JetHT: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_JetHT2016X.root
jet_fakerate_plots_JetHT: selection=HLTjet_qcd
jet_fakerate_plots_JetHT: fakerate_ratio_plots

.PHONY:jet_fakerate_plots_SingleMu
jet_fakerate_plots_SingleMu: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleMuon2016X.root
jet_fakerate_plots_SingleMu: selection=HLTmu_wjets
jet_fakerate_plots_SingleMu: fakerate_ratio_plots

jet_fakerate_plots: jet_fakerate_plots_JetHT jet_fakerate_plots_SingleMu

mumu_fakerate_ratio_plots: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
mumu_fakerate_ratio_plots: selection=mumu_passjets
mumu_fakerate_ratio_plots: fakerate_ratio_plots

elmu_fakerate_ratio_plots: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elmu_fakerate_ratio_plots: selection=elmu_passjets
elmu_fakerate_ratio_plots: fakerate_ratio_plots

elel_fakerate_ratio_plots: file=test/tests/outdir_${job}/merged-sets/jobsums/Data13TeV_SingleLepton2016X.root
elel_fakerate_ratio_plots: selection=elel_passjets
elel_fakerate_ratio_plots: fakerate_ratio_plots

dilep_fakerate_ratio_plots: elel_fakerate_ratio_plots mumu_fakerate_ratio_plots elmu_fakerate_ratio_plots

jet_mc_distrs:
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr x 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr y 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTjet_qcd_jets_distr z 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr x 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr y 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	stacked_histo_project3to1_distr --normalize 1.7 HLTmu_wjets_jets_distr z 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_mc_fakerates:
	jet_fakerate_per_mc 0.25 HLTjet_qcd  z 1  0   0.4 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTjet_qcd  y 1 -2.5 2.5 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTjet_qcd  x 1  0 150   JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_mc 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_mc 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_mc 0.25 HLTmu_wjets x 1  0 150   SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_mc_fakerates_dilep:
	jet_fakerate_per_mc 0.25 mumu_passjets  z 1  0   0.4 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 mumu_passjets  y 1 -2.5 2.5 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 mumu_passjets  x 1  0 150   mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets z 1  0   0.4  elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets y 1 -2.5 2.5  elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_mc 0.25 elmu_passjets x 1  0 150    elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_qcd:
	jet_fakerate_per_origin 0.25 HLTjet_qcd z 1  0   0.4 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_per_origin 0.25 HLTjet_qcd x 1  0 150   JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`

jet_origin_fakerates_wjets:
	jet_fakerate_per_origin 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_per_origin 0.25 HLTmu_wjets x 1  0 150   SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_origin_fakerates: jet_origin_fakerates_qcd jet_origin_fakerates_wjets

jet_origin_fakerates_elmu:
	jet_fakerate_per_origin 0.25 elmu_passjets z 1  0   0.4 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 elmu_passjets y 1 -2.5 2.5 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 elmu_passjets x 1  0 150   elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_mumu:
	jet_fakerate_per_origin 0.25 mumu_passjets z 1  0   0.4 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 mumu_passjets y 1 -2.5 2.5 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 mumu_passjets x 1  0 150   mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

jet_origin_fakerates_elel:
	jet_fakerate_per_origin 0.25 elel_passjets z 1  0   0.4 elel test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 elel_passjets y 1 -2.5 2.5 elel test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_per_origin 0.25 elel_passjets x 1  0 150   elel test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins:
	jet_fakerate_separate_origins 0.25 HLTjet_qcd z 1  0   0.4 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTjet_qcd x 1  0 150   JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins 0.25 HLTmu_wjets x 1  0 150   SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`

jet_fakerate_separate_origins_overlayed:
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd z 1  0   0.4 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd y 1 -2.5 2.5 JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTjet_qcd x 1  0 150   JetHT test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_JetHT`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets z 1  0   0.4 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets y 1 -2.5 2.5 SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	jet_fakerate_separate_origins_overlayed 0.25 HLTmu_wjets x 1  0 150   SingleMuon test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}_SingleMuon`
	#jet_fakerate_separate_origins_overlayed.cc

jet_fakerate_separate_origins_dilep:
	jet_fakerate_separate_origins 0.25 mumu_passjets z 1  0   0.4 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets y 1 -2.5 2.5 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 mumu_passjets x 1  0 150   mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets z 1  0   0.4 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets y 1 -2.5 2.5 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins 0.25 elmu_passjets x 1  0 150   elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

jet_fakerate_separate_origins_overlayed_dilep:
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets z 1  0   0.4 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets y 1 -2.5 2.5 mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 mumu_passjets x 1  0 150   mumu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets z 1  0   0.4 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets y 1 -2.5 2.5 elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`
	jet_fakerate_separate_origins_overlayed 0.25 elmu_passjets x 1  0 150   elmu test/tests/outdir_"${job}"/merged-sets/ `cat dtags_${job}`

overlay_std_fakerates:
	overlay_th1d_distrs ${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3
	overlay_th1d_distrs ${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3
	overlay_th1d_distrs ${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3

overlay_fakerates:
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3

overlay_fakerates_4:
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-x-pt_overlayed.png     qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_x.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_x.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-y-eta_overlayed.png    qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_y.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_y.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs ${ttjob}_${jetjob}_fake-rates-per-z-radius_overlayed.png qcd-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_JetHT2016X_HLTjet_qcd_tau_jets_distr_over_HLTjet_qcd_jets_distr_z.root h3 wjets-sel test/tests/outdir_${jetjob}/merged-sets/jobsums/Data13TeV_SingleMuon2016X_HLTmu_wjets_tau_jets_distr_over_HLTmu_wjets_jets_distr_z.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_z.root h3

overlay_dilep_fakerates:
	overlay_th1d_distrs ${ttjob}_dileps_fake-rates-per-x-pt_overlayed.png  elel-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_x.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_x.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_x.root h3
	overlay_th1d_distrs ${ttjob}_dileps_fake-rates-per-y-eta_overlayed.png  elel-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_y.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_y.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_y.root h3
	overlay_th1d_distrs ${ttjob}_dileps_fake-rates-per-z-redius_overlayed.png  elel-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elel_passjets_tau_jets_distr_over_elel_passjets_jets_distr_z.root h3 mumu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_mumu_passjets_tau_jets_distr_over_mumu_passjets_jets_distr_z.root h3 elmu-sel test/tests/outdir_${ttjob}/merged-sets/jobsums/Data13TeV_SingleLepton2016X_elmu_passjets_tau_jets_distr_over_elmu_passjets_jets_distr_z.root h3



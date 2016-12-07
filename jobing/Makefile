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

extract_fakerates:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} singleel_pretauselection_jettotaufakerates singlemu_pretauselection_jettotaufakerates

extract_pileup:
	./jobing/extract_csv.zsh ${job} jobing/dtags/${dtags} pileup_passtrig_rawWeight_pergoodpv pileup_passtrig_weight_pergoodpv

get_dtags:
	ls test/tests/outdir_"${job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > dtags_${job}

selection_control_histos:
	for h in `cat cur_th1d_histos_pt_eta_energy`;     \
	do                                                \
	stacked_histo_distr ${lumi} $$h 4 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done
	for h in `cat cur_th1d_histos_other`;             \
	do                                                \
	stacked_histo_distr ${lumi} $$h 1 test/tests/outdir_"${job}"/merged-sets/ `cat dtags_"${job}"`; \
	done


done_ck:
	ls test/tests/outdir_"${job}"/*job_done | wc -l
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

make_proxy:
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
	ls test/tests/outdir_"${jet_job}"/merged-sets/*root | sed -e 's_.*sets/__' -e 's/.root$$//' > test_dtags
	for d in `cat test_dtags | sed -n '/^MC/p' `; do jet_origins test/tests/outdir_"${jet_job}"/merged-sets/ $$d; done

.PHONY:fakerate_ratio_plots
fakerate_ratio_plots:
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root -a -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	#test -s test/tests/outdir_"${jet_job}"/merged-sets/${file}
	test -s ${file}
	# make stops if a command in the recipe fails...
	# || { echo foo; exit 1 }
	histo_project3toN_ratio_distr x pT     ${selection}_tau_jets_distr ${selection}_jets_distr 1  0 150   ${file}
	histo_project3toN_ratio_distr y eta    ${selection}_tau_jets_distr ${selection}_jets_distr 1 -2.5 2.5 ${file}
	histo_project3toN_ratio_distr z radius ${selection}_tau_jets_distr ${selection}_jets_distr 1  0   0.3 ${file}

jet_fakerate_plots_JetHT: file=test/tests/outdir_${job}/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root
jet_fakerate_plots_JetHT: selection=HLTjet_qcd
jet_fakerate_plots_JetHT: fakerate_ratio_plots

jet_fakerate_plots_SingleMu: file=test/tests/outdir_${job}/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
jet_fakerate_plots_SingleMu: selection=HLTmu_wjets
jet_fakerate_plots_SingleMu: fakerate_ratio_plots

jet_fakerate_plots: jet_fakerate_plots_JetHT jet_fakerate_plots_SingleMu

jet_fakerate_plots_old:
	test -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root -a -s test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	# make stops if a command in the recipe fails...
	# || { echo foo; exit 1 }
	histo_project3toN_ratio_distr x HLTjet_qcd_tau_jets_distr HLTjet_qcd_jets_distr 1  0 150   test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root
	histo_project3toN_ratio_distr y HLTjet_qcd_tau_jets_distr HLTjet_qcd_jets_distr 1 -2.5 2.5 test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root
	histo_project3toN_ratio_distr z HLTjet_qcd_tau_jets_distr HLTjet_qcd_jets_distr 1  0   0.3 test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_JetHT2016X_23Sep2016_v1.root
	#
	histo_project3toN_ratio_distr x HLTmu_wjets_tau_jets_distr HLTmu_wjets_jets_distr 1  0 150   test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	histo_project3toN_ratio_distr y HLTmu_wjets_tau_jets_distr HLTmu_wjets_jets_distr 1 -2.5 2.5 test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root
	histo_project3toN_ratio_distr z HLTmu_wjets_tau_jets_distr HLTmu_wjets_jets_distr 1  0   0.3 test/tests/outdir_"${jet_job}"/merged-sets/Data13TeV_SingleMuon2016X_23Sep2016_v1.root



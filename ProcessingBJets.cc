
struct bTaggingEfficiencyHistograms {
	TH2F* b_alljet   ;
	TH2F* b_tagged   ;
	TH2F* c_alljet   ;
	TH2F* c_tagged   ;
	TH2F* udsg_alljet;
	TH2F* udsg_tagged;
	};

double bTagging_b_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.b_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.b_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.b_tagged->GetBinContent(global_bin_id);
	return N_tagged/N_alljets;
	}

double bTagging_c_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& eta, double& pt)
	{
	Int_t global_bin_id = bEffs.c_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.c_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.c_tagged->GetBinContent(global_bin_id);
	return N_tagged/N_alljets;
	}

double bTagging_udsg_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& eta, double& pt)
	{
	Int_t global_bin_id = bEffs.udsg_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.udsg_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.udsg_tagged->GetBinContent(global_bin_id);
	return N_tagged/N_alljets;
	}

int processBJets_BTag(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
	BTagCalibrationReader& btagCal, BTagSFUtil& btsfutil,
	struct bTaggingEfficiencyHistograms& bEffs,
	string& b_tagger_label, float b_tag_WP,
	pat::JetCollection& selBJets,                          // output
	bool record, bool debug) // more output

{
// for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)

for (size_t ijet = 0; ijet < jets.size(); ++ijet)
	{
	pat::Jet& jet = jets[ijet];

	double eta=jet.eta();
	double pt=jet.pt();

	bool hasCSVtag(jet.bDiscriminator(b_tagger_label) > b_tag_WP);
	bool raw_CSV_tag = hasCSVtag;
	//bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.935);
	bool hasCSVtag_BTagUp(false), hasCSVtag_BTagDown(false);

	if (record)
		{
		fill_2d(string("btag_b_tagging_candidate_jets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		fill_1d(string("btag_b_tagging_candidate_jets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		if (hasCSVtag) fill_2d(string("btag_b_tagging_candidate_jets_pt_eta_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		}

	//update according to the SF measured by BTV
	// new fency procedure with CSV files
	// 80X calibrators in btagCal
	// usage:
	// btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, b_jet.eta(), b_jet.pt());
	// -- one calibrator for central value, and up/down
	// thus the specification of the value to callibrate,
	// instead of different callibrators
	if(isMC){
		// int flavId=jet.partonFlavour();
		int flavId=jet.hadronFlavour();
		// also: patJet->genParton().pdgId()
		// fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

		double sf = 1.0, eff = 1.0;
		/* If the jet is tagged -- weight *= SF of the jet
		 * if not weight *= (1 - eff*SF)/(1 - eff)
		 */

		if (abs(flavId)==5) {
			if (record)
				{
				/* recording the jets, per-flavour, candidates and tagged --- for the efficiency measurement and for control
				 * NOTICE: the recorded jets are weighted with the current genWeight, pile-up etc
				 */
				fill_2d(string("btag_b_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				if (hasCSVtag)
					fill_2d(string("btag_b_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			// get SF for the jet
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_b"), 200, 0., 2.,   sf, weight);

			// get eff for the jet
			eff = bTagging_b_jet_efficiency(bEffs, eta, pt);
			if (record) fill_1d(string("btag_eff_flavour_b"), 200, 0., 2.,   sf, weight);
			}
		else if(abs(flavId)==4) {

			if (record)
				{ 
				fill_2d(string("btag_c_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				if (hasCSVtag)
					fill_2d(string("btag_c_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_c"), 200, 0., 2.,   sf, weight);

			eff = bTagging_c_jet_efficiency(bEffs, eta, pt);
			if (record) fill_1d(string("btag_eff_flavour_c"), 200, 0., 2.,   sf, weight);
			}
		else {

			if (record)
				{ 
				fill_2d(string("btag_udsg_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				if (hasCSVtag)
					fill_2d(string("btag_udsg_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_udsg"), 200, 0., 2.,   sf, weight);

			eff = bTagging_udsg_jet_efficiency(bEffs, eta, pt);
			if (record) fill_1d(string("btag_eff_flavour_udsg"), 200, 0., 2.,   sf, weight);
			}

	if (hasCSVtag) // a tagged jet
		bTaggingSF_eventWeight *= sf;
	else // not tagged
		bTaggingSF_eventWeight *= (1 - sf*eff) / (1 - eff);
	}

	if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
		{
		selBJets.push_back(jet);
		if (record)
			{
			fill_2d(string("btag_control_jet_selBJets_pt_eta"), 250, 0., 500., 200, -4., 4., pt, jet.eta(), weight);
			fill_1d(string("btag_control_jet_selBJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		}
	}


return 0;
}






int processBJets_BTag(pat::JetCollection& jets, bool isMC, double weight, // input
	BTagCalibrationReader& btagCal, BTagSFUtil& btsfutil,
	string& b_tagger_label, float b_tag_WP,
	float beff, float leff,
	pat::JetCollection& selBJets,                          // output
	bool record, bool debug) // more output
{

// for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
for (size_t ijet = 0; ijet < jets.size(); ++ijet)
	{
	pat::Jet& jet = jets[ijet];

	double eta=jet.eta();

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

		double sf;
		if (abs(flavId)==5) {
			// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
			// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_b_jets_pt_eta_beforesf"), jet.pt(), eta, weight);
			// int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight)

			if (record)
				{ 
				fill_2d(string("btag_b_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				if (hasCSVtag)
					{
					fill_2d(string("btag_b_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
					//fill_2d(string("btag_b_tagged_b_jets_pt_eta_beforesf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
					// jet_radius(pat::Jet& jet)
					//fill_1d(string("btag_b_tagged_b_jets_radius"), 400, 0., 4.,   jet_radius(jet), weight);
					}
				}

			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, jet.pt(), 0.);
			if (record) fill_1d(string("btag_sf_flavour_b"), 200, 0., 2.,   sf, weight);

			btsfutil.modifyBTagsWithSF(hasCSVtag, sf, beff);

			// recording bulk effect of the SF on N b jets
			fill_1d(string("btag_b_falvour_SF_effect"), 4, 0., 5.,   (raw_CSV_tag ? 0 : 2) + (hasCSVtag ? 0 : 1), weight);
			// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
			// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
			// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_b_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
			if (record && hasCSVtag) fill_2d(string("btag_b_tagged_b_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		} else if(abs(flavId)==4) {

			if (record)
				{ 
				fill_2d(string("btag_c_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				if (hasCSVtag)
					fill_2d(string("btag_c_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				}

			// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, jet.pt(), 0.);
			// fill_btag_sf(string("btag_sf_flavour_c"), sf, weight);
			if (record) fill_1d(string("btag_sf_flavour_c"), 200, 0., 2.,   sf, weight);

			btsfutil.modifyBTagsWithSF(hasCSVtag, sf, beff);
			fill_1d(string("btag_c_falvour_SF_effect"), 4, 0., 5.,   (raw_CSV_tag ? 0 : 2) + (hasCSVtag ? 0 : 1), weight);
			// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
			// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
			// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_c_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
			if (record && hasCSVtag) fill_2d(string("btag_b_tagged_c_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		} else {

			if (record)
				{ 
				fill_2d(string("btag_udsg_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				if (hasCSVtag)
					fill_2d(string("btag_udsg_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				}

			// hasCSVtag = false; // FIXME: 8.26-27 no b-tagged light jets in MC <<<<<-----------------------------------

			// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCalL.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, jet.pt(), 0.);
			// fill_btag_sf(string("btag_sf_flavour_udsg"), sf, weight);
			if (record) fill_1d(string("btag_sf_flavour_udsg"), 200, 0., 2.,   sf, weight);

			btsfutil.modifyBTagsWithSF(hasCSVtag, sf, leff);
			fill_1d(string("btag_udsg_falvour_SF_effect"), 4, 0., 5.,   (raw_CSV_tag ? 0 : 2) + (hasCSVtag ? 0 : 1), weight);
			// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
			// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
			// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_udsg_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
			if (record && hasCSVtag) fill_2d(string("btag_b_tagged_udsg_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		}
	if (record && hasCSVtag) fill_2d(string("btag_b_tagged_candidate_jets_pt_eta_tagged_after_mc_sfs"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
	}

	if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
		{
		selBJets.push_back(jet);
		if (record)
			{
			fill_2d(string("btag_control_jet_selBJets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("btag_control_jet_selBJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		}
	}


return 0;
}



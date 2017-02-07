


int processMuons_ID_ISO_Kinematics(pat::MuonCollection& muons, reco::Vertex goodPV, double weight, // input
	patUtils::llvvMuonId::MuonId mu_ID, patUtils::llvvMuonId::MuonId veto_mu_ID,               // config/cuts
	patUtils::llvvMuonIso::MuonIso mu_ISO, patUtils::llvvMuonIso::MuonIso veto_mu_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	pat::MuonCollection& selMuons, LorentzVector& muDiff, unsigned int& nVetoMu,               // output
	bool record, bool debug) // more output
{
//LorentzVector muDiff(0., 0., 0., 0.);
// std::vector<patUtils::GenericLepton> selLeptons;
//pat::MuonCollection selMuons;
//unsigned int nVetoMu(0);
// unsigned int count_idiso_muons = 0;

for(unsigned int count_idiso_muons = 0, n=0; n<muons.size (); ++n)
	{
	// patUtils::GenericLepton& lepton = leptons[n];
	pat::Muon& muon = muons[n];

	if (record)
		{
		fill_2d(string("control_mu_slimmedmions_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
		fill_1d(string("control_mu_slimmedmions_phi"), 128, -3.2, 3.2, muon.phi(), weight);
		}

	bool
		passKin(true),     passId(true),     passIso(true),
		passVetoKin(true), passVetoId(true), passVetoIso(true);

	int lid(muon.pdgId()); // should always be 13

	//apply muon corrections
	/* no lepcorrs in 13.6
	if(abs(lid) == 13 && muCor)
		{
		float qter;
		TLorentzVector p4(muon.px(), muon.py(), muon.pz(), muon.energy());
		// old corrections:
		// muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
		// if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
		// roch-cor (rochcor) corrections:
		if (isMC) muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
		else muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
		muDiff -= muon.p4();
		muon.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
		muDiff += muon.p4();
		}
	*/


	//no need for charge info any longer
	//lid = abs(lid);
	//TString lepStr(lid == 13 ? "mu" : "e");
			
	// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
	// no need to mess with photon ID // double minDRlg(9999.);
	// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
	// no need to mess with photon ID // if(minDRlg<0.1) continue;

	//Cut based identification

	// ------------------------- muon IDs
	passId     = patUtils::passId(muon, goodPV, mu_ID,      patUtils::CutVersion::ICHEP16Cut);
	passVetoId = patUtils::passId(muon, goodPV, veto_mu_ID, patUtils::CutVersion::ICHEP16Cut);

	// ------------------------- muon isolation
	passIso     = patUtils::passIso(muon, mu_ISO,      patUtils::CutVersion::ICHEP16Cut);
	passVetoIso = patUtils::passIso(muon, veto_mu_ISO, patUtils::CutVersion::ICHEP16Cut);

	if (record && passId && passIso)
		{
		fill_2d(string("control_mu_idiso_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
		//fill_pt_e( string("all_muons_idiso_pt"), muon.pt(), weight);
		if (count_idiso_muons < 2)
			fill_2d(string("control_mu_idiso_top2pt_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
			//fill_pt_e( string("top2pt_muons_idiso_pt"), muon.pt(), weight);
		count_idiso_muons += 1;
		}


	// ---------------------------- Muon Kinematics
	double leta( muon.eta());

	// ---------------------- Main muon kin
	if(muon.pt() < pt_cut)   passKin = false;
	// if(leta > 2.1)        passKin = false;
	//if(muon.pt() < 26.)   passKin = false;
	if(leta > eta_cut)        passKin = false;

	// ---------------------- Veto muon kin
	//if (muon.pt () < 20)  passVetoKin = false;
	// if (leta > 2.1)       passVetoKin = false;
	if (muon.pt() < veto_pt_cut)   passVetoKin = false;
	if (leta > veto_eta_cut)       passVetoKin = false;

	if     (passKin     && passId     && passIso)
		{
		selMuons.push_back(muon);
		if (record)
			{
			fill_2d(string("control_mu_selMuons_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
			fill_1d(string("control_mu_selMuons_phi"), 128, -3.2, 3.2, muon.phi(), weight);
			}
		}
	else if(passVetoKin && passVetoId && passVetoIso) nVetoMu++;

	// TODO: ------------ ID/Iso scale factors:
	// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Muon_reconstruction_identificati
	// if it is MC and ID tight muon is found -- multiply weight by factor,
	// corresponding to its' tight scale factor distr
	// if loose -- loose scale factor

	/* TODO: add ID/Iso and run
	if (isMC) {
		if (passId) {
			// NOTE: scale factors are given for absolute eta
			weight *= muID_Tight_sf.GetBinContent( muID_Tight_sf.FindBin(leta, muon.pt()) );
		} else if (passVetoId) {
			weight *= muID_Loose_sf.GetBinContent( muID_Loose_sf.FindBin(leta, muon.pt()) );
		}
		if (passIso) {
			// NOTE: scale factors are given for absolute eta
			weight *= muIso_Tight_sf.GetBinContent( muIso_Tight_sf.FindBin(leta, muon.pt()) );
		} else if (passVetoIso) {
			weight *= muIso_Loose_sf.GetBinContent( muIso_Loose_sf.FindBin(leta, muon.pt()) );
		}
	}
	*/
	}


// TODO: there should be no need to sort selected muons here again -- they are in order of Pt
std::sort(selMuons.begin(),   selMuons.end(),   utils::sort_CandidatesByPt);

return 0;
}



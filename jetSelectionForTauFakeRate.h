#ifndef JETSEL_H
#define JETSEL_H

/*
 * couple useful functions for jet-tau fakerates
 */

#include <iostream>
#include <map>
#include <string>

#include "TH3F.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

using namespace std;



extern string mc_decay;

extern std::map<std::pair <string,string>, TH1D> th1d_distr_control;
extern std::map<string, TH1D> th1d_distr_control_headers;

extern std::map<std::pair <string,string>, TH2D> th2d_distr_control;
extern std::map<string, TH2D> th2d_distr_control_headers;

extern std::map<string, std::map<string, TH3D>> th3d_distr_maps_control;
extern std::map<string, TH3D> th3d_distr_maps_control_headers;

extern int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight);
extern int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight);


// TODO: I wonder where Float_t is defined?
// but probably TH3F depends on it and pulls it in anyway

// good bins 1, 2
// Float_t bins_pt[11] = { 0, 29, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
//Float_t bins_pt[11] = { 0, 30, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
//static Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 80, 150, 500 }; // 11 bins, 12 edges
//static int n_bins_pt = 11;
// bins 3 (exactly from AN489)
//Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 500 }; // 13 bins, 14 edges
//int n_bins_pt = 13;
// bins 4 (based on j6.5 distributions)
//static Float_t bins_pt[15] = { 0, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 150, 500 }; // 14 bins 15 edges
//static int n_bins_pt = 14;

// 0.1 pt bins
static Float_t bins_pt[49] = { 0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 65, 70, 75, 80, 100, 150, 500 }; // 48 bins 49 edges
static int n_bins_pt = 48;




// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
//Float_t bins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
//int n_bins_eta = 7;
//static Float_t bins_eta[10] = { -3, -2.4, -2.3, -1.5, -0.45, 0.45, 1.5, 2.3, 2.4, 3 }; // 9 bins 10 edges
//static int n_bins_eta = 9;
// AN bins go in 0.5, and are done for abs eta
// I will use 0.5, but for both signs of eta:
//static Float_t bins_eta[15] = { -3, -2.4, -2.3, -2., -1.5, -1, -0.5, 0, 0.5, 1., 1.5, 2., 2.3, 2.4, 3 }; // 14 bins 15 edges
//static int n_bins_eta = 14;

// now (2016) we've got more statistics -- thus smaller eta bins
// 0.1 eta bins
static Float_t bins_eta[51] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 }; // 50 bins 51 edges
static int n_bins_eta = 50;


//Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	//0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges
//int n_bins_rad = 15;
// exactly from AN489:
//static Float_t bins_rad[12] = { 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 1, 2 }; // 11 bins 12 edges
//static int n_bins_rad = 11;

// 0.1 radius bins:
static Float_t bins_rad[36] = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.4, 1, 2 }; // 35 bins 36 edges
static int n_bins_rad = 35;


//float tau_fake_distance = 0.1; // first-try distance
float tau_fake_distance = 0.3; // the distance to tau for a jet to be considered tau's origin


double jet_radius(pat::Jet& jet)
	{
	//return sqrt(jet.EtaPhiMoments::etaEtaMoment + jet.EtaPhiMoments::phiPhiMoment);
	//return sqrt(jet.etaEtaMoment() + jet.phiPhiMoment());
	return sqrt(jet.etaetaMoment() + jet.phiphiMoment());
	}

double jet_radius(pat::Tau& jet)
	{
	//return sqrt(jet.EtaPhiMoments::etaEtaMoment + jet.EtaPhiMoments::phiPhiMoment);
	//return sqrt(jet.etaEtaMoment() + jet.phiPhiMoment());
	return sqrt(jet.etaetaMoment() + jet.phiphiMoment());
	}



int fill_jet_distr(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th3d_distr_maps_control
	if (th3d_distr_maps_control.find(mc_decay) == th3d_distr_maps_control.end() )
		{
		//
		th3d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH3D>()));
		}

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad);
		// th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}


	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, radius, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad)));
		}

	// return success:
	return 0;
	}



int record_jets_fakerate_distrs(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, vector<LorentzVector>& visible_gen_taus, double event_weight, bool isMC)
	{

	for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
		{
		pat::Jet& jet = selJets[ijet];
		// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
		// selection jets
		fill_jet_distr(channel + selection + ("_jets_distr"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
		//fill_3d(channel + selection + ("_jets_distr"), 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad, 300, 0, 300,   20, event_weight);

		// study jet_R ~ jet_eta -- does calculation of R depend on eta of jet?
		// (maybe it would be nice to use bins of the 3d jet distr?)
		fill_2d(channel + selection + ("_jet_radius_vs_eta"), 100, 0., 2., 120, -3., 3., jet_radius(jet), jet.eta(), event_weight);

		// const reco::GenParticle* genParton()
		// jet parton origin
		if (isMC)
			{
			//const reco::GenParticle* jet_origin = jet.genParton();
			// the ID should be in:
			// jet_origin->pdgId();
			int jet_origin = abs(jet.partonFlavour());

			// match to visible_gen_taus
			// if partonFlavour == 0 and matches to a tau (dR < 1)
			// add it as 15 --- tau

			double minDRtj (9999.);
			unsigned int closest_tau = -1;
			for (size_t i = 0; i < visible_gen_taus.size(); i++)
				{
				double jet_tau_distance = TMath::Min(minDRtj, reco::deltaR (jet, visible_gen_taus[i]));
				if (jet_tau_distance<minDRtj)
					{
					closest_tau = i;
					minDRtj = jet_tau_distance;
					}
				}
			if (jet_origin == 0 && (minDRtj < 1)) jet_origin = 15;

			fill_1d(channel + selection + ("_jet_partonFlavour"), 100, 0, 100,   jet_origin, event_weight);
			fill_1d(channel + selection + ("_jet_hadronFlavour"), 100, 0, 100,   abs(jet.hadronFlavour()), event_weight);

			// per-jet-origin distrs
			if (jet_origin == 21) // gluons
				fill_jet_distr(channel + selection + ("_jets_distr_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 5) // b-quarks
				fill_jet_distr(channel + selection + ("_jets_distr_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 0) // other stuff (taus probably)
				fill_jet_distr(channel + selection + ("_jets_distr_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 15) // the matched to gen_taus (dR < 1) jets
				fill_jet_distr(channel + selection + ("_jets_distr_t"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else // other (light) quarks
				fill_jet_distr(channel + selection + ("_jets_distr_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			}

		for(size_t itau=0; itau < selTaus.size(); ++itau)
			{
			// selection taus (fake taus)
			//for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			//{
			//double fake_distance = reco::deltaR(selJets[ijet], selTaus[itau]);
			double fake_distance = reco::deltaR(jet, selTaus[itau]);
			//qcd_taujet_distance->Fill(fake_distance);
			fill_1d(channel + selection + ("_taujet_distance"), 100, 0, 2,   fake_distance, event_weight);

			if (fake_distance <= tau_fake_distance)
				{
				// the tau is fake by this jet -- save distr
				//qcd_tau_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
				fill_jet_distr(channel + selection + ("_tau_jets_distr"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
				// fill_pt_pt(string control_point_name, double pt1, double pt2, double event_weight)
				//fill_pt_pt(channel + selection + ("_tau_taujet_pts"), jet.pt(), selTaus[itau].pt(), event_weight);
				fill_2d(channel + selection + ("_tau_taujet_pts"), 400, 0., 400., 400, 0., 400, jet.pt(), selTaus[itau].pt(), event_weight);

				// N tau-jets
				//increment( channel + selection + ("_selection_ntaujet"), event_weight );

				// const reco::GenParticle* genParton()
				// jet parton origin for faking jet (just in case)
				if (isMC)
					{
					int jet_origin = abs(jet.partonFlavour());
					//const reco::GenParticle* jet_origin = selJets[ijet].genParton();
					// the ID should be in:
					// jet_origin->pdgId();
					//qcd_taujet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(channel + selection + ("_taujet_origins"), 100, 0, 100,   jet_origin, event_weight);
					//qcd_taujet_origin->Fill(abs( jet_origin->pdgId() ));
					// qcd_jet_origin

					// per-jet-origin distrs
					if (jet_origin == 21) // gluons
						fill_jet_distr(channel + selection + ("_tau_jets_distr_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 5) // b-quarks
						fill_jet_distr(channel + selection + ("_tau_jets_distr_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 0) // other stuff (taus probably)
						fill_jet_distr(channel + selection + ("_tau_jets_distr_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else // other (light) quarks
						fill_jet_distr(channel + selection + ("_tau_jets_distr_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					}
				continue;
				}
			//}
			}
		}

	// return success
	return 0;
	}



#endif /* JETSEL_H */


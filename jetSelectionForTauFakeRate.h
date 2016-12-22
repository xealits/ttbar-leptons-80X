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
static Float_t bins_pt[15] = { 0, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 150, 500 }; // 14 bins 15 edges
static int n_bins_pt = 14;


// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
//Float_t bins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
//int n_bins_eta = 7;
static Float_t bins_eta[10] = { -3, -2.4, -2.3, -1.5, -0.45, 0.45, 1.5, 2.3, 2.4, 3 }; // 9 bins 10 edges
static int n_bins_eta = 9;

//Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	//0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges
//int n_bins_rad = 15;
// exactly from AN489:
static Float_t bins_rad[12] = { 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 1, 2 }; // 11 bins 12 edges
static int n_bins_rad = 11;


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



int record_jets_fakerate_distrs(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, double event_weight, bool isMC)
	{

	for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
		{
		pat::Jet& jet = selJets[ijet];
		// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
		// selection jets
		fill_jet_distr(channel + selection + ("_jets_distr"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
		//fill_3d(channel + selection + ("_jets_distr"), 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad, 300, 0, 300,   20, event_weight);

		// const reco::GenParticle* genParton()
		// jet parton origin
		if (isMC)
			{
			//const reco::GenParticle* jet_origin = jet.genParton();
			// the ID should be in:
			// jet_origin->pdgId();
			//qcd_jet_origin->Fill(abs( jet.partonFlavour() ));
			fill_1d(channel + selection + ("_jet_partonFlavour"), 100, 0, 100,   abs(jet.partonFlavour()), event_weight);
			fill_1d(channel + selection + ("_jet_hadronFlavour"), 100, 0, 100,   abs(jet.hadronFlavour()), event_weight);
			//fill_1d(channel + selection + ("_jet_genParton_pdgId"), 100, 0, 100,   abs(jet.genParton()->pdgId()), event_weight);
			//qcd_jet_origin->Fill(abs( jet_origin->pdgId() ));
			// qcd_taujet_origin
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
					//const reco::GenParticle* jet_origin = selJets[ijet].genParton();
					// the ID should be in:
					// jet_origin->pdgId();
					//qcd_taujet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(channel + selection + ("_taujet_origins"), 100, 0, 100,   abs(jet.partonFlavour()), event_weight);
					//qcd_taujet_origin->Fill(abs( jet_origin->pdgId() ));
					// qcd_jet_origin
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


//
// Oleksii Toldaiev, <oleksii.toldaiev@gmail.com>
//
// ttbar to leptons event selection
// CMSSW version 80X (8_0_X, 8_0_5 for example)
//

//#define STANDALONE 1 // it seems to grab all the memory

#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
//#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

//#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"

// should work in CMSSW_8_0_12 and CMSSW_8_1_0
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
// #include "CondFormats/BTauObjects/interface/BTagCalibration.h"
// #include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
// this one is for 80X -> #include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h" 
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibrator.h"  
//#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibrator.h" 

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

#include <map>
#include <string>

#include "UserCode/ttbar-leptons-80X/interface/recordFuncs.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingBJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingMuons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingElectrons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingTaus.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingDRCleaning.h"

#include "UserCode/ttbar-leptons-80X/interface/SystematicShifts.h"

using namespace std;

namespace utils
	{
	namespace cmssw
		{
		// TODO: it is the same jetCorrector as in MacroUtils, only Fall_ prefix is set
		// Fall15_25nsV2_
		FactorizedJetCorrector* getJetCorrector(TString baseDir, TString pf, bool isMC)
			{
			gSystem->ExpandPathName(baseDir);
			//TString pf(isMC ? "MC" : "DATA");
			// TString pf("Fall15_25nsV2_");
			//pf += (isMC ? "MC" : "DATA");

			//order matters: L1 -> L2 -> L3 (-> Residuals)
			std::vector<std::string> jetCorFiles;
			std::cout << baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt" << std::endl;
			jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
			if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// now there is a practically empty file Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt
			// adding the run on it anyway
			//jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// it is dummy/empty file for MC and apparently is is not used
			// but in v13.1 it seemed to influence selection a bit
			// adding it for v13.4 -- will test later without it
			// and removing in 13.7 test -- compare with 13.4 & 13.4_repeat

			//init the parameters for correction
			std::vector<JetCorrectorParameters> corSteps;
			for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
			//return the corrector
			return new FactorizedJetCorrector(corSteps);
			}

		enum METvariations { NOMINAL, JERUP, JERDOWN, JESUP, JESDOWN, UMETUP, UMETDOWN, LESUP, LESDOWN };

		std::vector<LorentzVector> getMETvariations(LorentzVector &rawMETP4, pat::JetCollection &jets, std::vector<patUtils::GenericLepton> &leptons,bool isMC)
			{
			std::vector<LorentzVector> newMetsP4(9,rawMETP4);
			if(!isMC) return newMetsP4;

			LorentzVector nullP4(0,0,0,0);
			//recompute the clustered and unclustered fluxes with energy variations
			for(size_t ivar=1; ivar<=8; ivar++)
				{
				//leptonic flux
				LorentzVector leptonFlux(nullP4), lepDiff(nullP4);
				for(size_t ilep=0; ilep<leptons.size(); ilep++)
					{
					LorentzVector lepton = leptons[ilep].p4();
					double varSign( (ivar==LESUP ? 1.0 : (ivar==LESDOWN ? -1.0 : 0.0) ) );
					int id( abs(leptons[ilep].pdgId()) );
					double sf(1.0);
					if(id==13) sf=(1.0+varSign*0.01);
					if(id==11)
						{
						if(fabs(leptons[ilep].eta())<1.442) sf=(1.0+varSign*0.02);
						else                                sf=(1.0-varSign*0.05);
						}
					leptonFlux += lepton;
					lepDiff += (sf-1)*lepton;
					}

				//clustered flux
				LorentzVector jetDiff(nullP4), clusteredFlux(nullP4);
				for(size_t ijet=0; ijet<jets.size(); ijet++)
					{
					if(jets[ijet].pt()==0) continue;
					double jetsf(1.0);
					// FIXME: change the way this is stored (to not storing it)
					/// if(ivar==JERUP)   jetsf=jets[ijet].getVal("jerup")/jets[ijet].pt();
					/// if(ivar==JERDOWN) jetsf=jets[ijet].getVal("jerdown")/jets[ijet].pt();
					/// if(ivar==JESUP)   jetsf=jets[ijet].getVal("jesup")/jets[ijet].pt();
					/// if(ivar==JESDOWN) jetsf=jets[ijet].getVal("jesdown")/jets[ijet].pt();
					//LorentzVector newJet( jets[ijet] ); newJet *= jetsf;
					LorentzVector newJet = jets[ijet].p4(); newJet *= jetsf;
					jetDiff       += (newJet-jets[ijet].p4());
					clusteredFlux += jets[ijet].p4();
					}
				LorentzVector iMet=rawMETP4-jetDiff-lepDiff;

				//unclustered flux
				if(ivar==UMETUP || ivar==UMETDOWN)
					{
					LorentzVector unclusteredFlux=-(iMet+clusteredFlux+leptonFlux);
					unclusteredFlux *= (ivar==UMETUP ? 1.1 : 0.9); 
					iMet = -clusteredFlux -leptonFlux - unclusteredFlux;
					}

				//save new met
				newMetsP4[ivar]=iMet;
				}

			//all done here
			return newMetsP4;
			}
		}
	}


bool hasLeptonAsDaughter(const reco::GenParticle p)
	{
	bool foundL(false);
	if(p.numberOfDaughters()==0) return foundL;

	// cout << "Particle " << p.pdgId() << " with status " << p.status() << " and " << p.numberOfDaughters() << endl;
	const reco::Candidate *part = &p;
	// loop on the daughter particles to check if it has an e/mu as daughter
	while ((part->numberOfDaughters()>0))
		{
		const reco::Candidate* DaughterPart = part->daughter(0);
		// cout << "\t\t Daughter: " << DaughterPart->pdgId() << " with status " << DaughterPart->status() << endl;
		if (fabs(DaughterPart->pdgId()) == 11 || fabs(DaughterPart->pdgId() == 13))
			{
			foundL = true;
			break;
			}
		part=DaughterPart;
		}
	return foundL;
	}


bool hasWasMother(const reco::GenParticle  p)
	{
	bool foundW(false);
	if(p.numberOfMothers()==0) return foundW;
	const reco::Candidate* part =&p; // (p.mother());
	// loop on the mother particles to check if it has a W as mother
	while ((part->numberOfMothers()>0))
		{
		const reco::Candidate* MomPart =part->mother();
		if (fabs(MomPart->pdgId())==24)
			{
			foundW = true;
			break;
			}
		part = MomPart;
		}
	return foundW;
	}

bool hasTauAsMother(const reco::GenParticle  p)
	{
	bool foundTau(false);
	if (p.numberOfMothers()==0) return foundTau;
	const reco::Candidate* part = &p; //(p.mother());
	// loop on the mother particles to check if it has a tau as mother
	while ((part->numberOfMothers()>0))
		{
		const reco::Candidate* MomPart =part->mother();
		if (fabs(MomPart->pdgId())==15)// && MomPart->status() == 2) // Not sure the status check is needed.
			{
			foundTau = true;
			break;
			}
		part = MomPart;
		}
	return foundTau;
	}


/*
 * return ID of lepton if the tau decays leptonically
 * if it's hadronic return 0
 */
int tau_leptonic_ID(const reco::Candidate* t, bool debug)
{
const reco::Candidate* tau = t;
while(true)
	{
	int n = tau->numberOfDaughters();
	//if (debug) cout << "tau n daughters = " << n << endl;
	switch(n)
		{
		case 0: return 0;
		case 1: // it should be (unheard of) tau->tau transition
			tau = tau->daughter(0);
			break; // there should not be no infinite loop here!
		default: /* so, it should be the decay
			  * check daughters of tau, if leptons 11/13 is among them -- return, if tau -- loop, if no tau and no leptons -- return 0, it's hadronic
			  */
			{
			// find lepton among daughters
			//if (debug) cout << "tau decay, daughters";
			int tau_daughter_i = -1;
			for (int i=0; i<n; i++)
				{
				//if (debug) cout << " (" << i;
				int d_id = fabs(tau->daughter(i)->pdgId());
				//if (debug) cout << ", " << d_id << ")";
				if (d_id == 11 || d_id == 13)
					return d_id; // if lepton is found -- return immediately
				else if (d_id == 15) // if tau is found -- store it's index for possible loop
					tau_daughter_i = i;
				}
			//if (debug) cout << endl;
			if (tau_daughter_i != -1)
				tau = tau->daughter(tau_daughter_i); // loop on tau
			else
				// if no 11/13 lepton was found and no tau -- it's hadronic decay
				return 0;
			}
		}
	}
}


/*
 * string find_W_decay(const reco::Candidate * W) {
 *
 * returns a substring {el, mu, tau, q} identifing which decay happend
 * or returns "" string if something weird happend
 *
 * Usage:
 * const reco::Candidate * W = p.daughter( W_num );
 * mc_decay += find_W_decay(W);
*/
string find_W_decay(const reco::Candidate * W, bool debug) {
	const reco::Candidate * p = W; // eeh.. how C passes const arguments? are they local to function
	int d0_id, d1_id; // ids of decay daughters
	//bool found_decay=false;
	// follow W, whatever happens to it (who knows! defensive programming here)
	// until leptonic/quarkonic decay is found
	// 
	// assume there can be only W->W transitions inbetween (TODO: to check actually)
	//while (!found_decay) {
	while (true) {
		int n = p->numberOfDaughters();
		switch(n) {
		case 0: return string("");
		case 1: // it should be another W->W transition
			p = p->daughter(0);
			break; // there should not be no infinite loop here!
		case 2: // so, it should be the decay
			//found_decay = true;
			d0_id = fabs(p->daughter(0)->pdgId());
			d1_id = fabs(p->daughter(1)->pdgId());
			if (d0_id == 15 || d1_id == 15 )
				{
				//if (debug) cout << "found tau decay of W" << endl;
				// find tau->lepton or tau->hadronic?
				const reco::Candidate* tau = (d0_id == 15? p->daughter(0) : p->daughter(1));
				//if (debug) cout << "got tau" << endl;
				int tau_decay_id = tau_leptonic_ID(tau, debug);
				//if (debug) cout << "found ID of tau decay = " << tau_decay_id << endl;
				if (tau_decay_id == 11)// leptonic return tauel or taumu
					return string("tauel");
				else if (tau_decay_id == 13)
					return string("taumu");
				else // hadronic tau -- just tau
					return string("tauh");
				}
			if (d0_id == 11 || d1_id == 11 ) return string("el");
			if (d0_id == 13 || d1_id == 13 ) return string("mu");
			return string("q"); // FiXME: quite dangerous control-flow!
		default: // and this is just crazy
			return string("");
		}
	}
}

// int n = p.numberOfDaughters();
// if (n == 2) { // it is a decay vertes of t to something
// int d0_id = p.daughter(0)->pdgId();
// int d1_id = p.daughter(1)->pdgId();






#define MULTISEL_SIZE 256
#define NPARTICLES_SIZE 50
#define NORMKINO_DISTR_SIZE 400


// JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

struct JobDef
{
	string type;
	string dtag;
	string job_num;
};


// inline CONTROL:

// map for TH1D and map for double (weight counters)

// inline control functions usage:
//   fill_pt_e( "control_point_name", value, weight)
//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
//   increment( "control_point_name", weight )
//   printout_distrs(FILE * out)
//   printout_counters(FILE * out)

// TODO: and the multiselect weight-flow?

//std::map<string, TH1D*> th1d_distr_control;


//jetToTauFakeRate(TH3D * tau_fake_rate_jets_histo, TH3D * tau_fake_rate_taus_histo, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]))
//jetToTauFakeRate(tau_fake_rate_jets_histo1, tau_fake_rate_taus_histo1, tau_fake_rate_jets_histo2, tau_fake_rate_taus_histo2, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]));
//jetToTauFakeRate(TH3D * tau_fake_rate_jets_histo1, TH3D * tau_fake_rate_taus_histo1, TH3D * tau_fake_rate_jets_histo2, TH3D * tau_fake_rate_taus_histo2, Double_t tau_fake_rate_histo1_fraction, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius)
// later tau_fake_rate_histo1_fraction can be a TH3D histogram with fractions per pt, eta, radius
double jetToTauFakeRate_vanila(TH3D * tau_fake_rate_jets_histo1, TH3D * tau_fake_rate_taus_histo1, TH3D * tau_fake_rate_jets_histo2, TH3D * tau_fake_rate_taus_histo2, Double_t tau_fake_rate_histo1_fraction, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	// the tau_fake_rate_jets_histo and tau_fake_rate_taus_histo
	// are identical TH3D histograms
	// Int_t TH1::FindBin 	( 	Double_t  	x,
	//	Double_t  	y = 0,
	//	Double_t  	z = 0 
	// )
	// virtual Double_t TH3::GetBinContent 	( 	Int_t  	bin	) 	const

	Int_t global_bin_id = tau_fake_rate_jets_histo1->FindBin(jet_pt, jet_eta, jet_radius);

	Double_t jets_rate1 = tau_fake_rate_jets_histo1->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate_taus_histo1->GetBinContent(global_bin_id);

	Double_t jets_rate2 = tau_fake_rate_jets_histo2->GetBinContent(global_bin_id);
	Double_t taus_rate2 = tau_fake_rate_taus_histo2->GetBinContent(global_bin_id);

	// now linear mix of the two fakerates, according to the fraction:
	//    tau_fake_rate_histo1_fraction * taus_rate1/jets_rate1 +
	//    (1 - tau_fake_rate_histo1_fraction) * taus_rate2/jets_rate2
	// and also for small rates:
	//    (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1)
	
	Double_t fakerate = tau_fake_rate_histo1_fraction * (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1) + (1 - tau_fake_rate_histo1_fraction) * (jets_rate2 < 1 ? 0 : taus_rate2/jets_rate2);

	if (debug)
		{
		cout << jet_pt << " " << jet_eta << " " << jet_radius << " : " << global_bin_id << " : ";
		cout << taus_rate1 << "/" << jets_rate1 << ", " << taus_rate2 << "/" << jets_rate2 << "; "<< fakerate << "\n";
		}

	return fakerate;
	}



double jetToTauFakeRate(TH3D * tau_fake_rate_jets_histo1, TH3D * tau_fake_rate_taus_histo1, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	// the tau_fake_rate_jets_histo and tau_fake_rate_taus_histo
	// are identical TH3D histograms
	// Int_t TH1::FindBin 	( 	Double_t  	x,
	//	Double_t  	y = 0,
	//	Double_t  	z = 0 
	// )
	// virtual Double_t TH3::GetBinContent 	( 	Int_t  	bin	) 	const

	Int_t global_bin_id = tau_fake_rate_jets_histo1->FindBin(jet_pt, jet_eta, jet_radius);

	Double_t jets_rate1 = tau_fake_rate_jets_histo1->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate_taus_histo1->GetBinContent(global_bin_id);

	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);

	if (debug)
		{
		cout << jet_pt << " " << jet_eta << " " << jet_radius << " : " << global_bin_id << " : ";
		cout << taus_rate1 << "/" << jets_rate1 << endl;
		}

	return fakerate;
	}






/* New fake rates
 * via extension of existing distributions
 * assuming pt, eta and radius distributions are independent.
 */
typedef struct {
	TH1D* x;
	TH1D* y;
	TH1D* z;
	double integral;
} FakeRateProjections;

double jetToTauFakeRate_Projections(
		FakeRateProjections & tau_fake_rate_jets_histo,
		FakeRateProjections & tau_fake_rate_taus_histo,
		Double_t tau_fake_rate_histo1_fraction,
		Double_t jet_pt, Double_t jet_eta, Double_t jet_radius,
		bool debug)

	{
	// this method has to be normalized to average fake rates
	// since the product multiplies fake rates 3 time the average becomes cubed
	double average_fake_rate = tau_fake_rate_taus_histo.integral / tau_fake_rate_jets_histo.integral;

	Int_t bin_x_id = tau_fake_rate_jets_histo.x->FindBin(jet_pt);
	Int_t bin_y_id = tau_fake_rate_jets_histo.y->FindBin(jet_eta);
	Int_t bin_z_id = tau_fake_rate_jets_histo.z->FindBin(jet_radius);

	Double_t jets_ratex = tau_fake_rate_jets_histo.x->GetBinContent(bin_x_id);
	Double_t taus_ratex = tau_fake_rate_taus_histo.x->GetBinContent(bin_x_id);
	Double_t jets_ratey = tau_fake_rate_jets_histo.y->GetBinContent(bin_y_id);
	Double_t taus_ratey = tau_fake_rate_taus_histo.y->GetBinContent(bin_y_id);
	Double_t jets_ratez = tau_fake_rate_jets_histo.z->GetBinContent(bin_z_id);
	Double_t taus_ratez = tau_fake_rate_taus_histo.z->GetBinContent(bin_z_id);

	// just multiply separate fake rates
	double fakerate = (taus_ratex/jets_ratex) * (taus_ratey/jets_ratey) * (taus_ratez/jets_ratez) / (average_fake_rate * average_fake_rate);

	/*
	if (debug)
		{
		cout << jet_pt << " " << jet_eta << " " << jet_radius << " : " << global_bin_id << " : ";
		cout << taus_rate1 << "/" << jets_rate1 << ", " << taus_rate2 << "/" << jets_rate2 << "; "<< fakerate << "\n";
		}
	*/

	return fakerate;
	}






// Top pT reweighting
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
// -- the study is done for 7-8 Tev
//    but still recommended for 13TeV

double top_pT_SF(double x)
	{
	// the SF function is SF(x)=exp(a+bx)
	// where x is pT of the top quark (at generation?)
	// sqrt(s) 	channel     	a     	b
	// 7 TeV 	all combined 	0.199 	-0.00166
	// 7 TeV 	l+jets      	0.174 	-0.00137
	// 7 TeV 	dilepton    	0.222 	-0.00197
	// 8 TeV 	all combined 	0.156 	-0.00137
	// 8 TeV 	l+jets       	0.159 	-0.00141
	// 8 TeV 	dilepton     	0.148 	-0.00129
	// 13 TeV	all combined	0.0615	-0.0005
	// -- taking all combined 13 TeV
	double a = 0.0615;
	double b = -0.0005;
	return TMath::Exp(a + b*x);
	}




//string mc_decay("");
// Some MC datasets are inclusive, but analysis needs info on separate channels from them
// thus the events are traversed and the channel is found
// currently it is done only for TTbar channel (isTTbarMC)
//
// the sub-channel of MC is paired together with the distr_name string into the key of control distrs
// it is then printed out with dtag

std::map<std::pair <string,string>, double> weight_flow_control;


std::map<std::pair <string,string>, TH1D> th1d_distr_control;
std::map<string, TH1D> th1d_distr_control_headers;

std::map<std::pair <string,string>, TH1I> th1i_distr_control;
std::map<string, TH1I> th1i_distr_control_headers;

std::map<std::pair <string,string>, TH2D> th2d_distr_control;
std::map<string, TH2D> th2d_distr_control_headers;






// TODO: rearrange the code into particle selection and channel selection
// TODO organize the code with new general record functions and remove these

int fill_btag_eff(string control_point_name, double pt, double eta, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th2d_distr_control.find(key) == th2d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th2d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		// th2d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 400.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th2d_distr_control.insert( std::make_pair(key, TH2D((mc_decay + control_point_name).c_str(), ";;Pt-Eta", 250, 0., 500., 200, -4., 4.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th2d_distr_control[key].Fill(pt, eta, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th2d_distr_control[control_point_name].Integral() << endl;

	if (th2d_distr_control_headers.find(string("btag_eff")) == th2d_distr_control_headers.end() )
		{
		// th2d_distr_control_headers[string("btag_eff")] = TH2D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th2d_distr_control_headers.insert( std::make_pair(string("btag_eff"), TH2D("Header of B-tag efficiency distributions", ";;Pt-Eta", 250, 0., 500., 200, -4., 4.)));
		}

	// return success:
	return 0;
	}


int fill_n(string control_point_name, unsigned int value, double weight)
	{
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;N", 100, 0., 100.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;N", 100, 0., 100.)) );
		}

	// fill the distribution:
	// th1i_distr_control[key].Fill(value, weight);
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	//if (th1i_distr_control_headers.find(string("n")) == th1i_distr_control_headers.end() )
	//	{
	//	th1i_distr_control_headers[string("n")] = TH1I("Header of particle counter distributions", ";;N", 100, 0., 100.);
	//	}
	if (th1d_distr_control_headers.find(string("n")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("n")] = TH1D("Header of particle counter distributions", ";;N", 100, 0., 100.);
		th1d_distr_control_headers.insert( std::make_pair(string("n"), TH1D("Header of particle counter distributions", ";;N", 100, 0., 100.)));
		}

	// return success:
	return 0;
	}


int fill_particle_ids(string control_point_name, int value, double weight)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1i_distr_control.find(key) == th1i_distr_control.end() )
	//if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;ID", 600, -300., 300.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th1i_distr_control.insert( std::make_pair(key, TH1I((mc_decay + control_point_name).c_str(), ";;ID", 600, -300., 300.)));
		//th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;ID", 600, -300., 300.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1i_distr_control[key].Fill(value, weight);
	//th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	if (th1i_distr_control_headers.find(string("p_id")) == th1i_distr_control_headers.end() )
	//if (th1d_distr_control_headers.find(string("p_id")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("p_id")] = TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.);
		th1i_distr_control_headers.insert( std::make_pair(string("p_id"), TH1I("Header of particle ID distributions", ";;ID", 600, -300., 300.)));
		//th1d_distr_control_headers.insert( std::make_pair(string("p_id"), TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.)));
		}

	// return success:
	return 0;
	}


int fill_pu(string control_point_name, double value, double weight)
	{
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;nVtx", 100, 0., 100.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;nVtx", 100, 0., 100.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1d_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("pu")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("pu")] = TH1D("Header of pile up distributions", ";;nVtx", 100, 0., 100.);
		th1d_distr_control_headers.insert( std::make_pair(string("pu"), TH1D("Header of pile up distributions", ";;nVtx", 100, 0., 100.)));
		}

	// return success:
	return 0;
	}


int fill_pt_e(string control_point_name, double value, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 400.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;Pt/E(GeV)", 400, 0., 400.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1d_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("pt_e")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("pt_e")] = TH1D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th1d_distr_control_headers.insert( std::make_pair(string("pt_e"), TH1D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.)));
		}

	// return success:
	return 0;
	}


int fill_eta(string control_point_name, double value, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;Eta", 200, -4., 4.)));
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);

	if (th1d_distr_control_headers.find(string("eta")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("eta")] = TH1D("Header of Eta distributions", ";;Eta", 200, -4., 4.);
			th1d_distr_control_headers.insert( std::make_pair(string("eta"), TH1D("Header of Eta distributions", ";;Eta", 200, -4., 4.)));
		}

	// return success:
	return 0;
	// TODO: return the return value of Fill call?
	}


int fill_btag_sf(string control_point_name, double value, double weight)
	{
	//
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;B_sf", 200, 0., 2.)));
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);

	if (th1d_distr_control_headers.find(string("b_sf")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("eta")] = TH1D("Header of Eta distributions", ";;Eta", 200, -4., 4.);
			th1d_distr_control_headers.insert( std::make_pair(string("b_sf"), TH1D("Header of B SF distributions", ";;Eta", 200, 0., 2.)));
		}

	// return success:
	return 0;
	}


int increment(string control_point_name, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (weight_flow_control.find(key) == weight_flow_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		// weight_flow_control[key] = 0.;
		weight_flow_control.insert( std::make_pair(key, 0.));
		}

	// fill the distribution:
	weight_flow_control[key] += weight;

	// return success:
	return 0;
	}




// Plain text output printing:

int printout_distrs(FILE * out, JobDef JD)
	{
	//th1d_distr_control
	fprintf(out, "distr_header:distr_name:header/content,data_type,dtag,job_num,distr_values\n");
	//th1d_distr_control_headers
	//th1i_distr_control_headers

	/*
	for(std::map<string, TH1D>::iterator it = th1d_distr_control_headers.begin(); it != th1d_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1D * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}
	*/

	for(std::map<string, TH1I>::iterator it = th1i_distr_control_headers.begin(); it != th1i_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1I * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}

	/*
	for(std::map<std::pair <string,string>, TH1D>::iterator it = th1d_distr_control.begin(); it != th1d_distr_control.end(); ++it)
		{
		// iterator->first = key
		// iterator->second = value

		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		TH1D * distr = & it->second;

		// // Content:
		fprintf(out, "%s:content, %s,%s,%s", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str());
		for (int i=0; i < distr->GetSize(); i++) fprintf(out, ",%g", distr->GetBinContent(i));
		fprintf(out, "\n");

		// New output:
		// fprintf(out, "%s:content", name.c_str());

		// JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str()
		//for (int i=0; i < distr->GetSize(); i++)
		//	{
		//	fprintf(out, "%s,%s,%g,%g\n", name.c_str(), prefix.c_str(), distr->GetBinCenter(i), distr->GetBinContent(i));
		//	}
		//fprintf(out, "\n");
		}
	*/

	for(std::map<std::pair <string,string>, TH1I>::iterator it = th1i_distr_control.begin(); it != th1i_distr_control.end(); ++it)
		{
		// iterator->first = key
		// iterator->second = value

		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		TH1I * distr = & it->second;

		// // Content:
		fprintf(out, "%s:content, %s,%s,%s", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str());
		for (int i=0; i < distr->GetSize(); i++) fprintf(out, ",%g", distr->GetBinContent(i));
		fprintf(out, "\n");

		// New output:
		// fprintf(out, "%s:content", name.c_str());

		// JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str()
		//for (int i=0; i < distr->GetSize(); i++)
		//	{
		//	fprintf(out, "%s,%s,%g,%g\n", name.c_str(), prefix.c_str(), distr->GetBinCenter(i), distr->GetBinContent(i));
		//	}
		//fprintf(out, "\n");
		}

	return 0;
	}


int printout_counters(FILE * out, JobDef JD)
	{
	// weight flow control
	//struct JobDef
	//{
	//	string type;
	//	string dtag;
	//	string job_num;
	//};

	// fprintf(out, "weight_flow_header:%s,%s,%s\n", JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str());
	fprintf(out, "counter_header:val_name,data_type,dtag,job_num,value\n");

	for(std::map<std::pair <string,string>, double>::iterator it = weight_flow_control.begin(); it != weight_flow_control.end(); ++it)
		{
		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		double weight_sum = it->second;
		fprintf(out, "%s,%s,%s,%s,%g\n", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str(), weight_sum);
		}
	fprintf(out, "weight_flow end\n");

	return 0;
	}









int main (int argc, char *argv[])
{
//##############################################
//########2    GLOBAL INITIALIZATION     ########
//##############################################

// check arguments
if (argc < 2)
	{
	std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
	exit (0);
	}
	
// load framework libraries
gSystem->Load ("libFWCoreFWLite");
AutoLibraryLoader::enable ();

// random numbers for corrections & uncertainties
TRandom3 *r3 = new TRandom3();

// configure the process
const edm::ParameterSet & runProcess = edm::readPSetsFrom (argv[1])->getParameter < edm::ParameterSet > ("runProcess");

bool debug           = runProcess.getParameter<bool>  ("debug");
int debug_len        = runProcess.getParameter<int>  ("debug_len");
bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
bool isMC            = runProcess.getParameter<bool>  ("isMC");
double xsec          = runProcess.getParameter<double>("xsec");
int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
TString dtag         = runProcess.getParameter<std::string>("dtag");
string dtag_s        = runProcess.getParameter<std::string>("dtag");
string job_num       = runProcess.getParameter<std::string>("job_num");

cout << "Job:" << endl << dtag_s << '_' << job_num << endl;
JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

TString outUrl = runProcess.getParameter<std::string>("outfile");
TString outdir = runProcess.getParameter<std::string>("outdir");

string  muHLT_MC1   = runProcess.getParameter<std::string>("muHLT_MC1"),   muHLT_MC2   = runProcess.getParameter<std::string>("muHLT_MC2"),
	muHLT_Data1 = runProcess.getParameter<std::string>("muHLT_Data1"), muHLT_Data2 = runProcess.getParameter<std::string>("muHLT_Data2");

cout << "Triggers:" << endl;
cout << muHLT_MC1 << '\t' << muHLT_MC2 << '\t' << muHLT_Data1 << '\t' << muHLT_Data2 << endl;

// lepton efficiency scale factors:
TString muon_effs_dirname     = runProcess.getParameter < std::string > ("muon_effs");
TString electron_effs_dirname = runProcess.getParameter < std::string > ("electron_effs");
gSystem->ExpandPathName(muon_effs_dirname    );
gSystem->ExpandPathName(electron_effs_dirname);

cout << "dirs with lepton efficiencies:" << endl;
cout << muon_effs_dirname << endl;
cout << electron_effs_dirname << endl;

// now, muons have
//      track(reconstruction) efficiency, which is recommended per eta of muon now (however there should be something about N vertices too..
//      trk->id eff (eff identify reconstructed muon)
//      id->iso
//      iso->trig
// -- in 2016 all this stuff is per run of datataking (also they call it "era" of datataking) before HIP fix (BCDEF) and after (GH)
//    I use full dataset, thus MC should be weighted to some average of full dataset
//    (randomly split MC events in portions corresponding to luminosity of each part)
//    thus I take average of two SF, weighted by luminosity of two eras
//
// the trig eff for dilepton case is: apply negative of it for both leptons
cout << "unpacking muon eff SFs" << endl;

TFile * muon_effs_tracking_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_tracking_more_BCDEF_fits.root").Data() );
TFile * muon_effs_tracking_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_tracking_more_GH_fits.root").Data() );
TGraphAsymmErrors* muon_effs_tracking_BCDEF_graph = (TGraphAsymmErrors*) muon_effs_tracking_BCDEF_file->Get("ratio_eff_aeta_dr030e030_corr");
TGraphAsymmErrors* muon_effs_tracking_GH_graph    = (TGraphAsymmErrors*) muon_effs_tracking_GH_file->Get("ratio_eff_aeta_dr030e030_corr");
cout << "Y tracking (reconstruction)" << endl;

TFile* muon_effs_id_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_BCDEF.root").Data() );
TFile* muon_effs_id_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_GH.root").Data() );
TH2D* muon_effs_id_BCDEF_histo = (TH2D*) ((TDirectoryFile*) muon_effs_id_BCDEF_file->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get("pt_abseta_ratio");
TH2D* muon_effs_id_GH_histo    = (TH2D*) ((TDirectoryFile*) muon_effs_id_GH_file   ->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get("pt_abseta_ratio");
cout << "Y id" << endl;

TFile* muon_effs_iso_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_BCDEF.root").Data() );
TFile* muon_effs_iso_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_GH.root").Data() );
TH2D* muon_effs_iso_BCDEF_histo = (TH2D*) ((TDirectoryFile*) muon_effs_iso_BCDEF_file->Get("TightISO_TightID_pt_eta"))->Get("pt_abseta_ratio");
TH2D* muon_effs_iso_GH_histo    = (TH2D*) ((TDirectoryFile*) muon_effs_iso_GH_file->   Get("TightISO_TightID_pt_eta"))->Get("pt_abseta_ratio");

// --- yep, everywhere here Tight ID and ISO is used, since that's the leptons I use

TFile* muon_effs_trg_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_RunBtoF.root").Data() );
TFile* muon_effs_trg_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_Period4.root").Data() );
TH2D* muon_effs_trg_BCDEF_histo = (TH2D*) ((TDirectoryFile*) muon_effs_trg_BCDEF_file->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins"))->Get("pt_abseta_ratio");
TH2D* muon_effs_trg_GH_histo    = (TH2D*) ((TDirectoryFile*) muon_effs_trg_GH_file->   Get("IsoMu24_OR_IsoTkMu24_PtEtaBins"))->Get("pt_abseta_ratio");
cout << "Y trigger" << endl;

// From run v9.45 (reduced TopTrig-recommended LumiCert, 32 fb^-1, missing the 2nd version of H?) the luminosity ranges are:

double SingleMuon_data_bcdef_fraction = 19716.274 / (19716.274 + 15931.028);
double SingleMuon_data_gh_fraction    = 15931.028 / (19716.274 + 15931.028);


// now, electrons have
//      track(reconstruction) efficiency, which is recommended per eta of muon now (however there should be something about N vertices too..
//      and ID sf
//      also trigger
//
// the trig eff for dilepton case is: apply negative of it for both leptons
cout << "unpacking electron eff SFs" << endl;

TFile* electron_effs_tracking_all_file = TFile::Open((electron_effs_dirname + "/2016_Sept23_ElectronReconstructionSF_egammaEffi.txt_EGM2D.root").Data() );
TH2D* electron_effs_tracking_all_histo = (TH2D*) electron_effs_tracking_all_file->Get("EGamma_SF2D");
cout << "Y tracking (reconstruction)" << endl;

// for the selected electrons, Tight ID
// not for Veto
TFile* electron_effs_id_all_file = TFile::Open((electron_effs_dirname + "/2016_Sept23_ElectronID_TightCutBased_egammaEffi.txt_EGM2D.root").Data() );
TH2D* electron_effs_id_all_histo = (TH2D*) electron_effs_id_all_file->Get("EGamma_SF2D");
cout << "Y id" << endl;

//analysis/electron-effs/2016_03Feb_TriggerSF_Run2016All_v1.root
TFile* electron_effs_trg_all_file = TFile::Open((electron_effs_dirname + "/2016_03Feb_TriggerSF_Run2016All_v1.root").Data() );
TH2D* electron_effs_trg_all_histo = (TH2D*) electron_effs_trg_all_file->Get("Ele27_WPTight_Gsf");
cout << "Y trigger" << endl;

// --- these SFs will be applied to the selected leptons independently




// Kino cuts
double kino_cuts_muon_pt          = runProcess.getParameter<double>("kino_cuts_muon_pt");
double kino_cuts_muon_eta         = runProcess.getParameter<double>("kino_cuts_muon_eta");
double kino_cuts_muon_veto_pt     = runProcess.getParameter<double>("kino_cuts_muon_veto_pt");
double kino_cuts_muon_veto_eta    = runProcess.getParameter<double>("kino_cuts_muon_veto_eta");

double kino_cuts_electron_pt          = runProcess.getParameter<double>("kino_cuts_electron_pt");
double kino_cuts_electron_eta         = runProcess.getParameter<double>("kino_cuts_electron_eta");
double kino_cuts_electron_veto_pt     = runProcess.getParameter<double>("kino_cuts_electron_veto_pt");
double kino_cuts_electron_veto_eta    = runProcess.getParameter<double>("kino_cuts_electron_veto_eta");

double jet_kino_cuts_pt          = runProcess.getParameter<double>("jet_kino_cuts_pt");
double jet_kino_cuts_eta         = runProcess.getParameter<double>("jet_kino_cuts_eta");
double tau_kino_cuts_pt          = runProcess.getParameter<double>("tau_kino_cuts_pt");
double tau_kino_cuts_eta         = runProcess.getParameter<double>("tau_kino_cuts_eta");

double jettaufr_jet_kino_cuts_pt          = runProcess.getParameter<double>("jettaufr_jet_kino_cuts_pt");
double jettaufr_jet_kino_cuts_eta         = runProcess.getParameter<double>("jettaufr_jet_kino_cuts_eta");
double jettaufr_tau_kino_cuts_pt          = runProcess.getParameter<double>("jettaufr_tau_kino_cuts_pt");
double jettaufr_tau_kino_cuts_eta         = runProcess.getParameter<double>("jettaufr_tau_kino_cuts_eta");

cout << "Kino cuts" << endl;
cout << "jets: (pt)\t" << jet_kino_cuts_pt << "\t(eta)" << jet_kino_cuts_eta << endl;
cout << "taus: (pt)\t" << tau_kino_cuts_pt << "\t(eta)" << tau_kino_cuts_eta << endl;

// Tau IDs:
string tau_decayMode = runProcess.getParameter<std::string>("tau_decayMode"),
	tau_ID       = runProcess.getParameter<std::string>("tau_ID"),
	tau_againstMuon     = runProcess.getParameter<std::string>("tau_againstMuon"),
	tau_againstElectron = runProcess.getParameter<std::string>("tau_againstElectron");

double tauIDsf          = runProcess.getParameter<double>("tauIDsf");
double tauIDsf_shift    = runProcess.getParameter<double>("tauIDsf_shift");

//("byLooseCombinedIsolationDeltaBetaCorr3Hits");
/* all tau MVA IDs:
 * byVLooseIsolationMVArun2v1DBoldDMwLT
 * byLooseIsolationMVArun2v1DBoldDMwLT
 * byMediumIsolationMVArun2v1DBoldDMwLT
 * byTightIsolationMVArun2v1DBoldDMwLT
 * byVTightIsolationMVArun2v1DBoldDMwLT
 */
string tau_VLoose_ID("byVLooseIsolationMVArun2v1DBoldDMwLT");
string tau_Loose_ID ("byLooseIsolationMVArun2v1DBoldDMwLT");
string tau_Medium_ID("byMediumIsolationMVArun2v1DBoldDMwLT");
string tau_Tight_ID ("byTightIsolationMVArun2v1DBoldDMwLT");
string tau_VTight_ID("byVTightIsolationMVArun2v1DBoldDMwLT");

cout << "Tau IDs:" << tau_decayMode << '\t' << tau_ID << '\t' << tau_againstMuon << '\t' << tau_againstElectron << "\t| " << tau_Loose_ID << "\t| " << tau_Tight_ID << endl;


string jetID_s = runProcess.getParameter<std::string>("jetID"),
	jetPUID_s = runProcess.getParameter<std::string>("jetPUID");
jet_id    jetID;
pu_jet_id jetPUID;

if (jetID_s == string("Loose"))
	jetID = LooseJET;
else if (jetID_s == string("Tight"))
	jetID = TightJET;

if      (jetPUID_s == string("LoosePU"))
	jetPUID = LoosePUJET;
else if (jetPUID_s == string("MediumPU"))
	jetPUID = MediumPUJET;
else if (jetPUID_s == string("TightPU"))
	jetPUID = TightPUJET;

bool with_PU            = runProcess.getParameter<bool>  ("with_PU");

cout << "Jet IDs: (main) " << jetID << '\t' << "(PU)" << jetPUID << "\t(with PU)" << with_PU << endl;

cout << "Output directory: " << outdir << endl;

const edm::ParameterSet& myVidElectronIdConf = runProcess.getParameterSet("electronidparas");
const edm::ParameterSet& myVidElectronMainIdWPConf = myVidElectronIdConf.getParameterSet("tight");
const edm::ParameterSet& myVidElectronVetoIdWPConf = myVidElectronIdConf.getParameterSet("loose");
	
VersionedPatElectronSelector electronVidMainId(myVidElectronMainIdWPConf);
VersionedPatElectronSelector electronVidVetoId(myVidElectronVetoIdWPConf);
	
TString suffix = runProcess.getParameter < std::string > ("suffix");
std::vector < std::string > urls = runProcess.getUntrackedParameter < std::vector < std::string > >("input");

// Good lumi mask
// v2
lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

// for new orthogonality [TESTING]
bool isSingleElectronDataset = !isMC && dtag.Contains ("SingleEle");
// old trigger dataset orthogonality:
bool
	filterOnlySINGLEE  (false),
	filterOnlySINGLEMU (false);
if (!isMC)
	{
	if (dtag.Contains ("SingleMuon")) filterOnlySINGLEMU = true;
	if (dtag.Contains ("SingleEle"))  filterOnlySINGLEE  = true;
	}

// it is not used now
bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));
// Reactivate for diboson shapes  
// bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
// bool isMC_WZ      (isMC && (string (dtag.Data ()).find ("TeV_WZ") != string::npos));

// adding W1Jets, W2Jets, W3Jets, W4Jets
// thus, the events with 1, 2, 3, 4 of WJets whould be excluded
// (what are the x-sections for each of these?)
bool isW0JetsSet   (isMC && (dtag.Contains ("W0Jets")));
// so, WJets provides only W0Jets events
// all the rest is in WNJets (W4Jets = W>=4Jets)
// the cross-section should probably be = WJets - sum(WNJets)
// and how to do it:
// < LHEEventProduct > lheEPHandle; lheEPHandle->hepeup().NUP  <-- this parameter
// it = 6 for W1Jets, 7 for W2Jets, 8 for W3Jets, 9 for W4Jets (why not >=9?)
// WJets have NUP = all those numbers
// for W0Jets one takes NUP = 5

//bool isNoHLT = runProcess.getParameter<bool>  ("isNoHLT");
bool isNoHLT = dtag.Contains("noHLT");

bool withTauIDSFs = runProcess.getParameter <bool> ("withTauIDSFs");

bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );

// Summer16_23Sep2016BCDV4_DATA_        Summer16_23Sep2016EFV4_DATA_        Summer16_23Sep2016GV4_DATA_        Summer16_23Sep2016HV4_DATA_
bool period_BCD = !isMC && (dtag.Contains("2016B") || dtag.Contains("2016C") || dtag.Contains("2016D"));
bool period_EF  = !isMC && (dtag.Contains("2016E") || dtag.Contains("2016F"));
bool period_G   = !isMC && (dtag.Contains("2016G"));
bool period_H   = !isMC && (dtag.Contains("2016H"));




	
//##############################################
//######## GET READY FOR THE EVENT LOOP ########
//##############################################
size_t totalEntries(0);


// Data-driven tau fakerate background FAKERATES

//TString bTaggingEfficiencies_filename   = runProcess.getParameter<std::string>("bTaggingEfficiencies");
//gSystem->ExpandPathName(bTaggingEfficiencies_filename);
//TFile* bTaggingEfficiencies_file = TFile::Open(bTaggingEfficiencies_filename.Data());

TString fakerate1_filename = runProcess.getParameter < std::string > ("dataDriven_tauFakeRates1");
TString fakerate2_filename = runProcess.getParameter < std::string > ("dataDriven_tauFakeRates2");
TString fakerate_dileptons_filename = runProcess.getParameter < std::string > ("dataDriven_tauFakeRates_dileptons");
gSystem->ExpandPathName(fakerate1_filename);
gSystem->ExpandPathName(fakerate2_filename);
gSystem->ExpandPathName(fakerate_dileptons_filename);

cout << "files with fake rates:" << endl;
cout << fakerate1_filename << endl;
cout << fakerate2_filename << endl;
cout << fakerate_dileptons_filename << endl;

TFile * tau_fake_rate1_file = TFile::Open(fakerate1_filename.Data() );
TFile * tau_fake_rate2_file = TFile::Open(fakerate2_filename.Data() );
TFile * tau_fake_rate_file_dileptons = TFile::Open(fakerate_dileptons_filename.Data());
cout << "fake rate QCD file: " << tau_fake_rate1_file << endl;

Double_t tau_fake_rate_histo1_fraction = runProcess.getParameter < Double_t > ("tau_fake_rate_histo1_fraction");

// LARGE BINS for normal/vanila rates
// TODO: so the two fakerates are done 2ce -- two files and each file has qcd and wjets jistos
// rate1 = file1 = JetHT data file
TH3D * tau_fake_rate1_jets_histo_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_jets_distr_large_bins");
TH3D * tau_fake_rate1_taus_histo_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_tau_jets_distr_large_bins");

cout << "fake rate QCD histogram: " << tau_fake_rate1_jets_histo_q << endl;
tau_fake_rate1_jets_histo_q->Print();

// rate2 = file2 = SingleMuon data file
TH3D * tau_fake_rate2_jets_histo_w = (TH3D *) tau_fake_rate2_file->Get("HLTmu_wjets_jets_distr_large_bins");
TH3D * tau_fake_rate2_taus_histo_w = (TH3D *) tau_fake_rate2_file->Get("HLTmu_wjets_tau_jets_distr_large_bins");

// dilepton fake rates
TH3D * tau_fake_rate_jets_histo_elmu = (TH3D *) tau_fake_rate_file_dileptons->Get("elmu_passjets_jets_distr_large_bins");
TH3D * tau_fake_rate_taus_histo_elmu = (TH3D *) tau_fake_rate_file_dileptons->Get("elmu_passjets_tau_jets_distr_large_bins");
TH3D * tau_fake_rate_jets_histo_mumu = (TH3D *) tau_fake_rate_file_dileptons->Get("mumu_passjets_jets_distr_large_bins");
TH3D * tau_fake_rate_taus_histo_mumu = (TH3D *) tau_fake_rate_file_dileptons->Get("mumu_passjets_tau_jets_distr_large_bins");
TH3D * tau_fake_rate_jets_histo_elel = (TH3D *) tau_fake_rate_file_dileptons->Get("elel_passjets_jets_distr_large_bins");
TH3D * tau_fake_rate_taus_histo_elel = (TH3D *) tau_fake_rate_file_dileptons->Get("elel_passjets_tau_jets_distr_large_bins");




// THE SMOOTH 3D PROJECTIONS METHOD

// TODO: so the two fakerates are done 2ce -- two files and each file has qcd and wjets jistos
// rate1 = file1 = JetHT data file
TH3D * tau_fake_rate1_jets_histoPROJECTIONS_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_jets_distr");
TH3D * tau_fake_rate1_taus_histoPROJECTIONS_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_tau_jets_distr");
TH3D * tau_fake_rate1_jets_histoPROJECTIONS_w = (TH3D *) tau_fake_rate1_file->Get("HLTmu_qcd_jets_distr");
TH3D * tau_fake_rate1_taus_histoPROJECTIONS_w = (TH3D *) tau_fake_rate1_file->Get("HLTmu_qcd_tau_jets_distr");


// rate2 = file2 = SingleMuon data file
TH3D * tau_fake_rate2_jets_histoPROJECTIONS_q = (TH3D *) tau_fake_rate2_file->Get("HLTmu_qcd_jets_distr");
TH3D * tau_fake_rate2_taus_histoPROJECTIONS_q = (TH3D *) tau_fake_rate2_file->Get("HLTmu_qcd_tau_jets_distr");
TH3D * tau_fake_rate2_jets_histoPROJECTIONS_w = (TH3D *) tau_fake_rate2_file->Get("HLTmu_wjets_jets_distr");
TH3D * tau_fake_rate2_taus_histoPROJECTIONS_w = (TH3D *) tau_fake_rate2_file->Get("HLTmu_wjets_tau_jets_distr");




// dilepton fake rates
TH3D * tau_fake_rate_jets_histoPROJECTIONS_elmu = (TH3D *) tau_fake_rate_file_dileptons->Get("elmu_passjets_jets_distr");
TH3D * tau_fake_rate_taus_histoPROJECTIONS_elmu = (TH3D *) tau_fake_rate_file_dileptons->Get("elmu_passjets_tau_jets_distr");
TH3D * tau_fake_rate_jets_histoPROJECTIONS_mumu = (TH3D *) tau_fake_rate_file_dileptons->Get("mumu_passjets_jets_distr");
TH3D * tau_fake_rate_taus_histoPROJECTIONS_mumu = (TH3D *) tau_fake_rate_file_dileptons->Get("mumu_passjets_tau_jets_distr");
TH3D * tau_fake_rate_jets_histoPROJECTIONS_elel = (TH3D *) tau_fake_rate_file_dileptons->Get("elel_passjets_jets_distr");
TH3D * tau_fake_rate_taus_histoPROJECTIONS_elel = (TH3D *) tau_fake_rate_file_dileptons->Get("elel_passjets_tau_jets_distr");

cout << "averages of fake rate distrs" << endl << "X,Y,Z projections" << endl;

cout << "qcd:\t";
cout << tau_fake_rate1_taus_histoPROJECTIONS_q->ProjectionX()->Integral() / tau_fake_rate1_jets_histoPROJECTIONS_q->ProjectionX()->Integral() << ',';
cout << tau_fake_rate1_taus_histoPROJECTIONS_q->ProjectionY()->Integral() / tau_fake_rate1_jets_histoPROJECTIONS_q->ProjectionY()->Integral() << ',';
cout << tau_fake_rate1_taus_histoPROJECTIONS_q->ProjectionZ()->Integral() / tau_fake_rate1_jets_histoPROJECTIONS_q->ProjectionZ()->Integral() << endl;
cout << "wjets:\t";
cout << tau_fake_rate2_taus_histoPROJECTIONS_w->ProjectionX()->Integral() / tau_fake_rate2_jets_histoPROJECTIONS_w->ProjectionX()->Integral() << ',';
cout << tau_fake_rate2_taus_histoPROJECTIONS_w->ProjectionY()->Integral() / tau_fake_rate2_jets_histoPROJECTIONS_w->ProjectionY()->Integral() << ',';
cout << tau_fake_rate2_taus_histoPROJECTIONS_w->ProjectionZ()->Integral() / tau_fake_rate2_jets_histoPROJECTIONS_w->ProjectionZ()->Integral() << endl;
cout << "mumu:\t";
cout << tau_fake_rate_taus_histoPROJECTIONS_mumu->ProjectionX()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_mumu->ProjectionX()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_mumu->ProjectionY()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_mumu->ProjectionY()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_mumu->ProjectionZ()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_mumu->ProjectionZ()->Integral() << endl;
cout << "elmu:\t";
cout << tau_fake_rate_taus_histoPROJECTIONS_elmu->ProjectionX()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elmu->ProjectionX()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_elmu->ProjectionY()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elmu->ProjectionY()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_elmu->ProjectionZ()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elmu->ProjectionZ()->Integral() << endl;
cout << "elel:\t";
cout << tau_fake_rate_taus_histoPROJECTIONS_elel->ProjectionX()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elel->ProjectionX()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_elel->ProjectionY()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elel->ProjectionY()->Integral() << ',';
cout << tau_fake_rate_taus_histoPROJECTIONS_elel->ProjectionZ()->Integral() / tau_fake_rate_jets_histoPROJECTIONS_elel->ProjectionZ()->Integral() << endl;

// and use? the Y integral -- all bins are in range

// Make fake rate projections for smooth procedure
//TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);

// QCD
TH1D* tau_fake_rate1_jets_histoPROJECTIONS_q_x = (TH1D*) tau_fake_rate1_jets_histoPROJECTIONS_q->Project3D("x");
TH1D* tau_fake_rate1_jets_histoPROJECTIONS_q_y = (TH1D*) tau_fake_rate1_jets_histoPROJECTIONS_q->Project3D("y");
TH1D* tau_fake_rate1_jets_histoPROJECTIONS_q_z = (TH1D*) tau_fake_rate1_jets_histoPROJECTIONS_q->Project3D("z");
TH1D* tau_fake_rate1_taus_histoPROJECTIONS_q_x = (TH1D*) tau_fake_rate1_taus_histoPROJECTIONS_q->Project3D("x");
TH1D* tau_fake_rate1_taus_histoPROJECTIONS_q_y = (TH1D*) tau_fake_rate1_taus_histoPROJECTIONS_q->Project3D("y");
TH1D* tau_fake_rate1_taus_histoPROJECTIONS_q_z = (TH1D*) tau_fake_rate1_taus_histoPROJECTIONS_q->Project3D("z");
FakeRateProjections frates_qcd_jets_proj;
frates_qcd_jets_proj.x = tau_fake_rate1_jets_histoPROJECTIONS_q_x;
frates_qcd_jets_proj.y = tau_fake_rate1_jets_histoPROJECTIONS_q_y;
frates_qcd_jets_proj.z = tau_fake_rate1_jets_histoPROJECTIONS_q_z;
frates_qcd_jets_proj.integral = tau_fake_rate1_jets_histoPROJECTIONS_q->Integral();
FakeRateProjections frates_qcd_taus_proj;
frates_qcd_taus_proj.x = tau_fake_rate1_taus_histoPROJECTIONS_q_x;
frates_qcd_taus_proj.y = tau_fake_rate1_taus_histoPROJECTIONS_q_y;
frates_qcd_taus_proj.z = tau_fake_rate1_taus_histoPROJECTIONS_q_z;
frates_qcd_taus_proj.integral = tau_fake_rate1_taus_histoPROJECTIONS_q->Integral();


// WJets
TH1D* tau_fake_rate2_jets_histoPROJECTIONS_w_x = (TH1D*) tau_fake_rate2_jets_histoPROJECTIONS_w->Project3D("x");
TH1D* tau_fake_rate2_jets_histoPROJECTIONS_w_y = (TH1D*) tau_fake_rate2_jets_histoPROJECTIONS_w->Project3D("y");
TH1D* tau_fake_rate2_jets_histoPROJECTIONS_w_z = (TH1D*) tau_fake_rate2_jets_histoPROJECTIONS_w->Project3D("z");
TH1D* tau_fake_rate2_taus_histoPROJECTIONS_w_x = (TH1D*) tau_fake_rate2_taus_histoPROJECTIONS_w->Project3D("x");
TH1D* tau_fake_rate2_taus_histoPROJECTIONS_w_y = (TH1D*) tau_fake_rate2_taus_histoPROJECTIONS_w->Project3D("y");
TH1D* tau_fake_rate2_taus_histoPROJECTIONS_w_z = (TH1D*) tau_fake_rate2_taus_histoPROJECTIONS_w->Project3D("z");
FakeRateProjections frates_wjets_jets_proj;
frates_wjets_jets_proj.x = tau_fake_rate2_jets_histoPROJECTIONS_w_x;
frates_wjets_jets_proj.y = tau_fake_rate2_jets_histoPROJECTIONS_w_y;
frates_wjets_jets_proj.z = tau_fake_rate2_jets_histoPROJECTIONS_w_z;
frates_wjets_jets_proj.integral = tau_fake_rate2_jets_histoPROJECTIONS_w->Integral();
FakeRateProjections frates_wjets_taus_proj;
frates_wjets_taus_proj.x = tau_fake_rate2_taus_histoPROJECTIONS_w_x;
frates_wjets_taus_proj.y = tau_fake_rate2_taus_histoPROJECTIONS_w_y;
frates_wjets_taus_proj.z = tau_fake_rate2_taus_histoPROJECTIONS_w_z;
frates_wjets_taus_proj.integral = tau_fake_rate2_taus_histoPROJECTIONS_w->Integral();

// ELMU
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elmu_x = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elmu->Project3D("x");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elmu_y = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elmu->Project3D("y");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elmu_z = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elmu->Project3D("z");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elmu_x = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elmu->Project3D("x");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elmu_y = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elmu->Project3D("y");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elmu_z = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elmu->Project3D("z");
FakeRateProjections frates_elmu_jets_proj;
frates_elmu_jets_proj.x = tau_fake_rate_jets_histoPROJECTIONS_elmu_x;
frates_elmu_jets_proj.y = tau_fake_rate_jets_histoPROJECTIONS_elmu_y;
frates_elmu_jets_proj.z = tau_fake_rate_jets_histoPROJECTIONS_elmu_z;
frates_elmu_jets_proj.integral = tau_fake_rate_jets_histoPROJECTIONS_elmu->Integral();
FakeRateProjections frates_elmu_taus_proj;
frates_elmu_taus_proj.x = tau_fake_rate_taus_histoPROJECTIONS_elmu_x;
frates_elmu_taus_proj.y = tau_fake_rate_taus_histoPROJECTIONS_elmu_y;
frates_elmu_taus_proj.z = tau_fake_rate_taus_histoPROJECTIONS_elmu_z;
frates_elmu_taus_proj.integral = tau_fake_rate_taus_histoPROJECTIONS_elmu->Integral();


// MUMU
TH1D* tau_fake_rate_jets_histoPROJECTIONS_mumu_x = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_mumu->Project3D("x");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_mumu_y = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_mumu->Project3D("y");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_mumu_z = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_mumu->Project3D("z");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_mumu_x = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_mumu->Project3D("x");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_mumu_y = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_mumu->Project3D("y");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_mumu_z = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_mumu->Project3D("z");
FakeRateProjections frates_mumu_jets_proj;
frates_mumu_jets_proj.x = tau_fake_rate_jets_histoPROJECTIONS_mumu_x;
frates_mumu_jets_proj.y = tau_fake_rate_jets_histoPROJECTIONS_mumu_y;
frates_mumu_jets_proj.z = tau_fake_rate_jets_histoPROJECTIONS_mumu_z;
frates_mumu_jets_proj.integral = tau_fake_rate_jets_histoPROJECTIONS_mumu->Integral();
FakeRateProjections frates_mumu_taus_proj;
frates_mumu_taus_proj.x = tau_fake_rate_taus_histoPROJECTIONS_mumu_x;
frates_mumu_taus_proj.y = tau_fake_rate_taus_histoPROJECTIONS_mumu_y;
frates_mumu_taus_proj.z = tau_fake_rate_taus_histoPROJECTIONS_mumu_z;
frates_mumu_taus_proj.integral = tau_fake_rate_taus_histoPROJECTIONS_mumu->Integral();

// ELEL
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elel_x = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elel->Project3D("x");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elel_y = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elel->Project3D("y");
TH1D* tau_fake_rate_jets_histoPROJECTIONS_elel_z = (TH1D*) tau_fake_rate_jets_histoPROJECTIONS_elel->Project3D("z");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elel_x = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elel->Project3D("x");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elel_y = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elel->Project3D("y");
TH1D* tau_fake_rate_taus_histoPROJECTIONS_elel_z = (TH1D*) tau_fake_rate_taus_histoPROJECTIONS_elel->Project3D("z");
FakeRateProjections frates_elel_jets_proj;
frates_elel_jets_proj.x = tau_fake_rate_jets_histoPROJECTIONS_elel_x;
frates_elel_jets_proj.y = tau_fake_rate_jets_histoPROJECTIONS_elel_y;
frates_elel_jets_proj.z = tau_fake_rate_jets_histoPROJECTIONS_elel_z;
frates_elel_jets_proj.integral = tau_fake_rate_jets_histoPROJECTIONS_elel->Integral();
FakeRateProjections frates_elel_taus_proj;
frates_elel_taus_proj.x = tau_fake_rate_taus_histoPROJECTIONS_elel_x;
frates_elel_taus_proj.y = tau_fake_rate_taus_histoPROJECTIONS_elel_y;
frates_elel_taus_proj.z = tau_fake_rate_taus_histoPROJECTIONS_elel_z;
frates_elel_taus_proj.integral = tau_fake_rate_taus_histoPROJECTIONS_elel->Integral();

//MC normalization (to 1/pb)
if(debug) cout << "DEBUG: xsec: " << xsec << endl;

// ------------------------------------- jet energy scale and uncertainties 
TString jecDir = runProcess.getParameter < std::string > ("jecDir");
gSystem->ExpandPathName (jecDir);
// v1
// getJetCorrector(TString baseDir, TString pf, bool isMC)
//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite

// in 2016 the corrections for data are per-period:
// Summer16_23Sep2016BCDV4_DATA_        Summer16_23Sep2016EFV4_DATA_        Summer16_23Sep2016GV4_DATA_        Summer16_23Sep2016HV4_DATA_
TString jet_corr_files;
if (isMC)
	jet_corr_files = "/Summer16_23Sep2016V4_MC";
else if (period_BCD)
	jet_corr_files = "/Summer16_23Sep2016BCDV4_DATA";
else if (period_EF)
	jet_corr_files = "/Summer16_23Sep2016EFV4_DATA";
else if (period_G)
	jet_corr_files = "/Summer16_23Sep2016GV4_DATA";
else if (period_H)
	jet_corr_files = "/Summer16_23Sep2016HV4_DATA";
FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector (jecDir, jet_corr_files, isMC);

JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + jet_corr_files + "_Uncertainty_AK4PFchs.txt").Data ());
// <-- used by smearJES
// (slimmedJets in MINIAOD are AK4PFchs jets)

//string jetResolutionFileName   = runProcess.getParameter<edm::FileInPath>("resolutionFile").fullPath();
//string jetResolutionSFFileName = runProcess.getParameter<edm::FileInPath>("scaleFactorFile").fullPath();

TString TjetResolutionFileName   = runProcess.getParameter<std::string>("resolutionFile");
TString TjetResolutionSFFileName = runProcess.getParameter<std::string>("scaleFactorFile");
gSystem->ExpandPathName(TjetResolutionFileName);
gSystem->ExpandPathName(TjetResolutionSFFileName);

string jetResolutionFileName   (TjetResolutionFileName);
string jetResolutionSFFileName (TjetResolutionSFFileName);
// <---- ROOT & CMSSW are best friends
JME::JetResolution jet_resolution_in_pt = JME::JetResolution(jetResolutionFileName);
JME::JetResolutionScaleFactor jet_resolution_sf_per_eta = JME::JetResolutionScaleFactor(jetResolutionSFFileName);

/*
 * from CondFormats/JetMETObjects/interface/JetResolutionObject.h
 *
 * enum class Variation {
 *   NOMINAL = 0,
 *   DOWN = 1,
 *   UP = 2
 *};
*/


// ------------------------------------- muon energy scale and uncertainties
// MuScleFitCorrector *muCor = NULL;
// FIXME: MuScle fit corrections for 13 TeV not available yet (more Zs are needed) getMuonCorrector (jecDir, dtag);
// TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
// gSystem->ExpandPathName(muscleDir);
// v1
//rochcor2015* muCor = new rochcor2015(); // This replaces the old MusScleFitCorrector that was used at RunI
// comes from:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
// last muon POG talk:
// https://indico.cern.ch/event/533054/contributions/2171540/attachments/1274536/1891597/rochcor_run2_MuonPOG_051716.pdf

// ------------------------------------- electron energy scale
// Electron energy scale, based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer and adapted to this framework
// v1
//string EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015";
//EpCombinationTool theEpCombinationTool;
//theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC 
//ElectronEnergyCalibratorRun2 ElectronEnCorrector(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
//ElectronEnCorrector.initPrivateRng(new TRandom(1234));


// --------------------------------------- b-tagging 

/* https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
 * https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
 * --- doing the simple 1) method, with event weight = P(Data)/P(MC)
 *     where P(Data) = Pow(SFi * effi) * Pow(1 - SFj * effj)
 *           p(MC)   = Pow(effi) * Pow(1 - effj)
 *     (so, the Pow(effi) cancels out)
 * the efficiencies are measured as in twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#b_tagging_efficiency_in_MC_sampl
 * per hadronFlavour
 * and passed to the job in the bTaggingEfficiencies parameter as the filename
 */

TString bTaggingEfficiencies_filename   = runProcess.getParameter<std::string>("bTaggingEfficiencies");
gSystem->ExpandPathName(bTaggingEfficiencies_filename);
TFile* bTaggingEfficiencies_file = TFile::Open(bTaggingEfficiencies_filename.Data());

cout << "b-tagging eff-s, filename: " << bTaggingEfficiencies_filename << endl;

TH2D* bTaggingEfficiencies_b_alljet   ;
TH2D* bTaggingEfficiencies_b_tagged   ;
TH2D* bTaggingEfficiencies_c_alljet   ;
TH2D* bTaggingEfficiencies_c_tagged   ;
TH2D* bTaggingEfficiencies_udsg_alljet;
TH2D* bTaggingEfficiencies_udsg_tagged;

TString backup_b_eff_distr("MC2016_Summer16_DYJetsToLL_50toInf_madgraph");

// ( ? : ) would look too much here
if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_b_hadronFlavour_candidates_tagged"))
	{
	bTaggingEfficiencies_b_alljet    = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_b_hadronFlavour_candidates");
	bTaggingEfficiencies_b_tagged    = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_b_hadronFlavour_candidates_tagged");
	}
else
	{
	bTaggingEfficiencies_b_alljet    = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates");
	bTaggingEfficiencies_b_tagged    = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates_tagged");
	}
if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_c_hadronFlavour_candidates_tagged"))
	{
	bTaggingEfficiencies_c_alljet    = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_c_hadronFlavour_candidates");
	bTaggingEfficiencies_c_tagged    = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_c_hadronFlavour_candidates_tagged");
	}
else
	{
	bTaggingEfficiencies_c_alljet    = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates");
	bTaggingEfficiencies_c_tagged    = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates_tagged");
	}
if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_udsg_hadronFlavour_candidates_tagged"))
	{
	bTaggingEfficiencies_udsg_alljet = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_udsg_hadronFlavour_candidates");
	bTaggingEfficiencies_udsg_tagged = (TH2D*) bTaggingEfficiencies_file->Get(dtag + "_btag_udsg_hadronFlavour_candidates_tagged");
	}
else
	{
	bTaggingEfficiencies_udsg_alljet = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates");
	bTaggingEfficiencies_udsg_tagged = (TH2D*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates_tagged");
	}

struct bTaggingEfficiencyHistograms bEffs;

bEffs.b_alljet    = bTaggingEfficiencies_b_alljet   ;
bEffs.b_tagged    = bTaggingEfficiencies_b_tagged   ;
bEffs.c_alljet    = bTaggingEfficiencies_c_alljet   ;
bEffs.c_tagged    = bTaggingEfficiencies_c_tagged   ;
bEffs.udsg_alljet = bTaggingEfficiencies_udsg_alljet;
bEffs.udsg_tagged = bTaggingEfficiencies_udsg_tagged;

// Prescriptions taken from:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X

// new btag calibration
// TODO: check callibration readers in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
// and latest standalone callibrator:
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h

// Setup calibration readers
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
// BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_76X.csv");
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_ichep_80X.csv");
BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");

// and there:

// in *80X
/*
The name of the measurements is

  "incl" for light jets,
  "mujets" (from QCD methods only) or “comb” (combination of QCD and ttbar methods)
    for b and c jets for what concerns the pT/eta dependence for the different WP for CSVv2.
    The precision of the “comb” measurement is better than for the “mujets”,
    however for precision measurements on top physics done in the 2lepton channel, it is recommended to use the "mujets" one.
  "ttbar" for b and c jets for what concerns the pT/eta dependence for the different WP for cMVAv2,
  but only to be used for jets with a pT spectrum similar to that in ttbar.
  The measurement "iterativefit" provides the SF as a function of the discriminator shape. 
*/

// only 1 reader:

BTagCalibrationReader btagCal(BTagEntry::OP_MEDIUM,  // operating point
// BTagCalibrationReader btagCal(BTagEntry::OP_TIGHT,  // operating point
                             "central",             // central sys type
                             {"up", "down"});      // other sys types
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_B,      // btag flavour
//            "comb");              // they say "comb" is better precision, but mujets are independent from ttbar dilepton channels
          "mujets");                //
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_C,      // btag flavour
          "mujets");               // measurement type
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_UDSG,   // btag flavour
            "incl");                // measurement type

/* usage in 80X:
double jet_scalefactor    = reader.eval_auto_bounds(
          "central", 
          BTagEntry::FLAV_B, 
          b_jet.eta(), 
          b_jet.pt()
      );
double jet_scalefactor_up = reader.eval_auto_bounds(
          "up", BTagEntry::FLAV_B, b_jet.eta(), b_jet.pt());
double jet_scalefactor_do = reader.eval_auto_bounds(
          "down", BTagEntry::FLAV_B, b_jet.eta(), b_jet.pt()); 

*/

// in 76X

//The name of the measurements is
//
//  * "incl" for light jets,
//  * "mujets" for b and c jets for what concerns the pT/eta dependence for the different WP for JP and CSVv2 and
//  * "ttbar" for b and c jets for what concerns the pT/eta dependence for the different WP for cMVAv2, but only to be used for jets with a pT spectrum similar to that in ttbar.
//  * The measurement "iterativefit" provides the SF as a function of the discriminator shape. 
// --- so "incl" instead of "comb" for light-quarks

// TODO: update btag CSVv2
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X#Data_MC_Scale_Factors
// OP_MEDIUM instead of OP_LOOSE?
// v1
//BTagCalibrationReader btagCal   (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "central");  // calibration instance, operating point, measurement type, systematics type
//BTagCalibrationReader btagCalUp (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "up"     );  // sys up
//BTagCalibrationReader btagCalDn (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "down"   );  // sys down
//BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_LOOSE, "comb", "central");  // calibration instance, operating point, measurement type, systematics type
//BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_LOOSE, "comb", "up"     );  // sys up
//BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_LOOSE, "comb", "down"   );  // sys down
//BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_MEDIUM, "incl", "central");  // calibration instance, operating point, measurement type, systematics type
//BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_MEDIUM, "incl", "up"     );  // sys up
//BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_MEDIUM, "incl", "down"   );  // sys down


/* TODO: CMSSW calibration:
The twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X#Supported_Algorithms_and_Operati says

The name of the measurements is
 * "incl" for light jets,
 * "mujets" for b and c jets for what concerns the pT/eta dependence for the different WP for JP and CSVv2 and 
 * "ttbar" for b and c jets for what concerns the pT/eta dependence for the different WP for cMVAv2, but only to be used for jets with a pT spectrum similar to that in ttbar.
 * The measurement "iterativefit" provides the SF as a function of the discriminator shape. 


// setup calibration readers
// Before CMSSW_8_0_11 and CMSSW_8_1_0, one reader was needed 
// for every OperatingPoint/MeasurementType/SysType combination.
BTagCalibration calib("csvv1", "CSVV1.csv");
BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
                             "central");           // systematics type
BTagCalibrationReader reader_up(BTagEntry::OP_LOOSE, "up");  // sys up
BTagCalibrationReader reader_do(BTagEntry::OP_LOOSE, "down");  // sys down

reader.load(&calib,               // calibration instance
            BTagEntry::FLAV_B,    // btag flavour
            "comb")               // measurement type
// reader_up.load(...)
// reader_down.load(...)



// Usage:

float JetPt = b_jet.pt(); bool DoubleUncertainty = false;
if (JetPt>MaxBJetPt)  { // use MaxLJetPt for  light jets
  JetPt = MaxBJetPt; 
  DoubleUncertainty = true;
}  

// Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
double jet_scalefactor = reader.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 
double jet_scalefactor_up =  reader_up.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 
double jet_scalefactor_do =  reader_do.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 

if (DoubleUncertainty) {
   jet_scalefactor_up = 2*(jet_scalefactor_up - jet_scalefactor) + jet_scalefactor; 
   jet_scalefactor_do = 2*(jet_scalefactor_do - jet_scalefactor) + jet_scalefactor; 
}

*/




// --------------------------------------- pileup weighting
// pile-up is done directly with direct_pileup_reweight
// v7-10 pile-up is back
std::vector<double> direct_pileup_reweight = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct");
std::vector<double> direct_pileup_reweight_up = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_up");
std::vector<double> direct_pileup_reweight_down = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_down");
	
gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

// generator token for event Gen weight call:
//edm::EDGetTokenT<GenEventInfoProduct> generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator")));


// ----------------------------
// So here we got all the parameters from the config


cout << "Some input parameters\n";
cout << "isMC = " << isMC << "\n";
cout << "isW0JetsSet = " << isW0JetsSet << "\n";
cout << "isTTbarMC = "    << isTTbarMC << "\n";
cout << "isNLOMC = "      << isNLOMC << "\n";
cout << "period BCD EF G H = " << period_BCD << " " << period_EF << " " << period_G << " " << period_H << endl;
cout << "jecDir = "      << jecDir << "\n";











//##############################################
//########    INITIATING HISTOGRAMS     ########
//##############################################

/*
 * now all output histograms are handled with recordFuncs -- fill_1d etc
 */
//TH1D* singlelep_ttbar_initialevents  = (TH1D*) new TH1D("singlelep_ttbar_init",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
//TH1D* singlelep_ttbar_preselectedevents = (TH1D*) new TH1D("singlelep_ttbar_presele",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 




//##############################################
//########           EVENT LOOP         ########
//##############################################
//loop on all the events

// this progress bar is idle now -- at times jobs were crashing on printing a dot in the progress bar
printf ("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");

int nMultiChannel(0);

unsigned int iev = 0;
double weightflow_control_el_selection = 0;
double weightflow_control_mu_selection = 0;
double weightflow_control_elel_selection = 0;
double weightflow_control_mumu_selection = 0;
double weightflow_control_elmu_selection = 0;

for(size_t f=0; f<urls.size();++f)
	{
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());
	cout << "Processing file: " << urls[f].c_str() << endl;

	iev++;

	int treeStep (ev.size()/50);

	for (ev.toBegin(); !ev.atEnd(); ++ev)
		{
		if(debug)
			{
			cout << "Processing event " << iev << "\n\n" ;
			if(iev > debug_len)
				{
				cout << "Got to the event " << iev << " in the file, exiting" << endl;
				//return 0;
				break;
				}
			}

		edm::EventBase const & iEvent = ev;

		// mc_decay = string("");
		mc_decay.clear();

		reco::GenParticleCollection gen;
		fwlite::Handle<reco::GenParticleCollection> genHandle;
		genHandle.getByLabel(ev, "prunedGenParticles");
		if(genHandle.isValid() ) gen = *genHandle;

                vector<LorentzVector> visible_gen_taus;
		if (isMC)
			for(size_t i = 0; i < gen.size(); ++ i)
				{
				const reco::GenParticle & p = gen[i];
				unsigned int id = fabs(p.pdgId());
				int st = p.status();
				int n_daughters = p.numberOfDaughters();

				// if it is a final state tau
				//  the status is 1 or 2
				//  1. (final state, not decays)
				//  2. (decayed or fragmented -- the case for tau)
				if (id == 15 && st == 1)
					visible_gen_taus.push_back(p.p4()); 
				else if (id == 15 && st == 2)
					{
					// it's a final state tau
					// select its' daughters, skipping neutrinos
					// add their momenta -- use the sum as a visible_gen_tau
					LorentzVector vis_ds(0,0,0,0);
					for (int j = 0; j < n_daughters; ++j)
						{
						const reco::Candidate * d = p.daughter(j);
						unsigned int d_id = fabs(d->pdgId());
						if (d_id == 12 || d_id == 14 || d_id == 16) continue;
						vis_ds += d->p4();
						}
					visible_gen_taus.push_back(vis_ds); 
					}
				}

		fwlite::Handle < LHEEventProduct > lheEPHandle;
		lheEPHandle.getByLabel (ev, "externalLHEProducer");

		if (debug && isMC)
			{
			cout << "number of gen particles = " << gen.size() << "\n";

			//mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
			cout << "nup = " << lheEPHandle->hepeup().NUP << "\n";
			}

		// take only W0Jets events from WJets set: (W0Jets have NUP == 5)
		if (isW0JetsSet && (lheEPHandle->hepeup().NUP != 5))
			continue;

		// -------------------------- trying to extract what decay was generated here
		// iEvent.getByLabel("genParticles", genParticles);
		// for(size_t i = 0; i < genParticles->size(); ++ i) {
		// }
		// listings of particles' mother-daughters  show
		// W status 44 can turn into the same W status 44 several times
		// W status 62 decays into leton-neutrino
		//
		// otozh
		// there are numerous ways to do it:
		//  - find motherless particles (protons in our case) and move from them
		//  - just check out two first particles in the collection -- it seems these are the mothers always
		//  - find t-quarks among particles and transverse their branches,
		//    hoping that t and tbar cannot happen in not a ttbar pair
		//
		// let's do the naive 3 version
		// targeting:
		//   find decay vertex t-W
		//   follow along W and find decay vertex to leptons or quarks
		// -- thus recursion is needed

		// traversing separate t quarks
		//if (debug) {
		if (isTTbarMC && genHandle.isValid()) {
			//string mc_decay(""); // move to all job parameters
			// every found t-quark decay branch will be added as a substring to this string (el, mu, tau, q)
			// TODO: what to do if there are > t-decays than 1 ttbar pair in an event? (naive traverse approach also doesn't work?)

			// For reference, some PDG IDs:
			// QUARKS
			// d  1
			// u  2
			// s  3
			// c  4
			// b  5
			// t  6
			// b' 7
			// t' 8
			// g 21
			// gamma 22
			// Z     23
			// W     24
			// h     25
			// e, ve     11, 12
			// mu, vmu   13, 14
			// tau, vtau 15, 16

			for(size_t i = 0; i < gen.size(); ++ i) {
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)

				if (id == 6) { // if it is a t quark
					// looking for W+ and its' leptonic decay
					// without checks, just checking 2 daughters
					int n = p.numberOfDaughters();
					if (n == 2) { // it is a decay vertes of t to something
						int d0_id = p.daughter(0)->pdgId();
						int d1_id = p.daughter(1)->pdgId();
						int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						mc_decay += find_W_decay(W, debug);
					}
				}

				if (id == -6) { // if it is a tbar quark
					// looking for W- and its' leptonic decay
					int n = p.numberOfDaughters();
					if (n == 2) { // it is a decay vertes of t to something
						int d0_id = p.daughter(0)->pdgId();
						int d1_id = p.daughter(1)->pdgId();
						int W_num = d0_id == -24 ? 0 : (d1_id == -24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						mc_decay += find_W_decay(W, debug);
						// mc_decay += string("bar");
						// removing the difference between two branches
						// to reduce the jobs produced for TTbar
					}
				}
			}
			// so, mc_decay will be populated with strings matching t decays
			// hopefully, in TTbar sample only ttbar decays are present and it is not ambigous
		}

		std::sort(mc_decay.begin(), mc_decay.end()); // no dif between elmu and muel etc

		if (debug) {
			cout << "MC suffix " << mc_decay << " is found\n";
			}

		//if (!mc_decay.empty()) mc_decay = string("_") + mc_decay;
		//mc_decay = string("_") + mc_decay; // so we'll have "_" or "_mcdecay"
		//mc_decay = mc_decay; // so we'll have "" or "mcdecay"

		//* List of mother-daughters for all particles
		//* TODO: make it into a separate function

		/*
		if (debug) {
			for(size_t i = 0; i < gen.size(); ++ i) {
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status();
				int n = p.numberOfDaughters();
				cout << i << ": " << id << " " << st;
				if (p.numberOfMothers() != 0) cout << " <- " ;
				for (int j = 0 ; j < p.numberOfMothers(); ++j) {
					const reco::Candidate * mom = p.mother(j);
					cout << " " << mom->pdgId() << " " << mom->status() << ";";
					}
				cout << "\n";
				if (n>0) {
					cout << "\t|-> " ;
					for (int j = 0; j < n; ++j) {
						const reco::Candidate * d = p.daughter( j );
						cout << d->pdgId() << " " << d->status() << "; " ;
						}
					cout << "\n";
					}
				}
			}
		*/

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )


		iev++;
		totalEntries++;
		/*
		if (iev % treeStep == 0)
			{
			printf (".");
			if(!debug) fflush (stdout); // Otherwise debug messages are flushed
			}
		*/

		edm::EventBase const & myEvent = ev;

		// ---------------------------------- MC shaping to data
		// MC is weighted according to distributions of a bunch of data properties

		// NLO -1 corrections
		double weight_Gen = 1.;

		// Top pT reweighting
		double weight_TopPT = 1.;
		// top_pT_SF

		// there are also the (dissabled now, since NLO samples are used) HT-binned and pthat-binned stitching of LO and NLO

		// ---------------------------------- pileup weight
		double weight_PU         (1.0);
		double weight_PU_up      (1.0);
		double weight_PU_down    (1.0);

		// ---------------------------------- b-tagging SF weight
		map<systematic_shift, double> weight_bTaggingSFs = {{SYS_NOMINAL, 1.0}, {SYS_BTAG_UP, 1.0}, {SYS_BTAG_DOWN, 1.0}};
		//double weight_bTaggingSF_up   (1.0);
		//double weight_bTaggingSF_down (1.0);

		// --------------------------------- tau ID SF
		double weight_tauIDsf         (1.0);
		double weight_tauIDsf_up      (1.0);
		double weight_tauIDsf_down    (1.0);
		double weight_without_tauIDsf (1.0);

		// rawWeight is everything but Pile-Up
		double rawWeight        (1.0);
		// full NOMINAL weight
		//double weight        (1.0);

		double HLT_efficiency_sf = 1.0;

		// final weight of the event
		// is map of all possible systematic shifts of weight
		//double weight           (1.0);
		std::map<systematic_shift, double> weights_FULL;
		//weights_FULL[SYS_NOMINAL] = 1.0;
		//weights_FULL[SYS_PU_UP]   = 1.0;
		//weights_FULL[SYS_PU_DOWN] = 1.0;
		//weights_FULL[SYS_TOP_PT] = 1.0;
		for ( const auto s : weightSystematics )
			{
			weights_FULL[s] = 1.0;
			}


		// --------------------------------------------------- RHO variables
		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;

		double rhoCentral = 0;
		fwlite::Handle<double> rhoCentralHandle;
		rhoCentralHandle.getByLabel(ev, "fixedGridRhoFastjetCentral");
		if(rhoCentralHandle.isValid() ) rhoCentral = *rhoCentralHandle;

		double rhoCentralNeutral = 0;
		fwlite::Handle<double> rhoCentralNeutralHandle;
		rhoCentralNeutralHandle.getByLabel(ev, "fixedGridRhoFastjetCentralNeutral");
		if(rhoCentralNeutralHandle.isValid() ) rhoCentralNeutral = *rhoCentralNeutralHandle;

		double rhoCentralChargedPileUp = 0;
		fwlite::Handle<double> rhoCentralChargedPileUpHandle;
		rhoCentralChargedPileUpHandle.getByLabel(ev, "fixedGridRhoFastjetCentralChargedPileUp");
		if(rhoCentralChargedPileUpHandle.isValid() ) rhoCentralChargedPileUp = *rhoCentralChargedPileUpHandle;



		// -------------------------------------------------- FIRST SECTION OF MC WEIGHTS, [1, 10]

		// increment( string("weightflow_n_miniaod_events"), 1.0 );
		// iniweight 1
		//fill_1d(string("weightflow_mu"), 300, 0, 300,   1, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_el"), 300, 0, 300,   1, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elel"), 300, 0, 300, 1, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elmu"), 300, 0, 300, 1, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_mumu"), 300, 0, 300, 1, weights_FULL[SYS_NOMINAL]);
		// Saving weight-flow for all systematic shifts
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   1, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   1, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 1, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 1, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 1, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   1, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   1, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 1, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 1, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 1, 1);

		// ---------------------------------- Top pT reeventing
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting<Paste>
		// find the produced tops,
		// reweight according them
		if (isTTbarMC && genHandle.isValid())
			{
			if (debug) cout << "finding tops for Top pT reweighting\n";

			bool found_top=false, found_atop=false;
			if (debug) cout << "found_top = " << found_top << ", found_atop = " << found_atop << "\n";
			for(size_t i = 0; i < gen.size(); ++ i)
				{
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)

				if (id == 6 && (!found_top))
					{ // if it is first t quark
					if (debug) cout << "found top, pT = " << p.pt() << "\n";
					found_top = true;
					weight_TopPT *= TMath::Sqrt(top_pT_SF(p.pt()));
					if (debug) cout << "now weight_TopPT = " << weight_TopPT << "\n";
					}

				if (id == -6 && (!found_atop))
					{ // if it is first anti-t quark
					if (debug) cout << "found atop, pT = " << p.pt() << "\n";
					found_atop = true;
					weight_TopPT *= TMath::Sqrt(top_pT_SF(p.pt()));
					if (debug) cout << "now weight_TopPT = " << weight_TopPT << "\n";
					}

				if (found_top && found_atop) break;
				// TODO: should I apply the reweighting only when both top-atop are found?
				// it is always true in TTbar MC..
				}
			}

		fill_1d(string("weight_TopPT"), 200, 0., 2., weight_TopPT, 1);

		// how is the overall integral of MC?
		// TOP PT is a systematic
		//weights_FULL[SYS_NOMINAL] *= weight_TopPT;
		//weights_FULL[SYS_PU_UP]   *= weight_TopPT;
		//weights_FULL[SYS_PU_DOWN] *= weight_TopPT;
		weights_FULL[SYS_TOP_PT] *= weight_TopPT;
		// the MC is lumi-xsec scaled to weightflow_weighted_miniaod_events

		// wighttoppt 2
		//fill_1d(string("weightflow_mu"), 300, 0, 300,   2, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_el"), 300, 0, 300,   2, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elel"), 300, 0, 300, 2, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elmu"), 300, 0, 300, 2, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_mumu"), 300, 0, 300, 2, weights_FULL[SYS_NOMINAL]);
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   2, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   2, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 2, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 2, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 2, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   2, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   2, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 2, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 2, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 2, 1);

		// old bit on top pT reweighting:
		// FIXME: Top pT reweighting to be reactivated as soon as corrections are released
		// if(tPt>0 && tbarPt>0 && topPtWgt)
		//   {
		//   topPtWgt->computeWeight(tPt,tbarPt);
		//   topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown);
		//   wgtTopPtUp /= wgtTopPt;
		//   wgtTopPtDown /= wgtTopPt;
		//   }

		// ---------------------------------- these are weird NLO -1 events
		// TODO: figure out how exactly they correct for NLO
		// Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
		if(isNLOMC)
			{
			if (debug) cout << "seting the NLOMC Gen weight\n";

			fwlite::Handle<GenEventInfoProduct> evt;
			evt.getByLabel(ev, "generator");

			// get by token maybe?
			/*
			edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
			generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator")))
			const edm::Event& iEvent
			edm::Handle<GenEventInfoProduct> evt;                                                                                           
			iEvent.getByToken( generatorToken_,evt);
			*/

			//edm::Handle<GenEventInfoProduct> evt;                                                                                           
			//iEvent.getByToken( generatorToken_,evt);

			if(evt.isValid())
				{
				if (debug) cout << "evt is valid, evt->weight() = " << evt->weight() << "\n";
				weight_Gen = (evt->weight() > 0 ) ? 1. : -1. ;
				// weight_Gen = evt->weight();
				}

			// FIXME: this is for PDF uncertainties, must reactivate it at some point.
			//fwlite::Handle<LHEEventProduct> lheEvtProd;
			//lheEvtProd.getByLabel(ev, "externalLHEProducer");
			//if(lheEvtProd.isValid())
			//  {
			//    weightLhe=lheEvtProd->originalXWGTUP();
			//    
			//   //for(unsigned int i=0; i<evet->weights().size();i++){
			//   //  double asdde=evet->weights()[i].wgt;
			//   //  EventInfo.ttbar_w[EventInfo.ttbar_nw]=EventInfo.ttbar_w[0]*asdde/asdd;
			//   //  EventInfo.ttbar_nw++;
			//   //}
			//  }
			//cout << "Event " << iev << " has genweight: " << weight_Gen << " and LHE weight " << weightLhe << endl;

			}


		//std::vector < TString > tags (1, "all"); // Inclusive inclusiveness

		//
		// DERIVE WEIGHTS TO APPLY TO SAMPLE
		//


		fill_1d(string("weight_Gen"), 200, -2., 2., weight_Gen, 1);

		// it's a nominal weight, i.e. weight all weights with it
		for ( const auto s : weightSystematics )
			{
			weights_FULL[s] *= weight_Gen;
			}
		rawWeight *= weight_Gen;

		// weightgen
		//fill_1d(string("weightflow_el"), 300, 0, 300,   3, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_mu"), 300, 0, 300,   3, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elel"), 300, 0, 300, 3, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_elmu"), 300, 0, 300, 3, weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("weightflow_mumu"), 300, 0, 300, 3, weights_FULL[SYS_NOMINAL]);
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   3, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   3, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 3, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 3, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 3, weight);
			}

		fill_1d(string("eventflow_el"), 300, 0, 300,   3, 1);
		fill_1d(string("eventflow_mu"), 300, 0, 300,   3, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 3, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 3, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 3, 1);

		// ------------------------------- PRIMARY VERTEX
		// needed for particle selection/event classification later
		// and pile-up control-distribution for data
		reco::VertexCollection vtx;
		reco::Vertex goodPV;
		unsigned int nGoodPV(0);
		fwlite::Handle<reco::VertexCollection> vtxHandle;
		vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
		if(vtxHandle.isValid() ) vtx = *vtxHandle;
		// Clean up vertex collection
		// it seems utils::isGoodVertex is outdated
		//nGoodPV = vtx.size();
		//goodPV = vtx[0];
		for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
			{
			//if(utils::isGoodVertex(vtx[ivtx]))
			// directly from rumors
			// some use:
			// * at least 4 degrees of freedome (ndof) (>=4) (!)
			// * Rho < 2 (impact parameter to the beam spot
			// * z < 24
			bool its_good = (!vtx[ivtx].isFake()) && vtx[ivtx].ndof() > 4 && abs(vtx[ivtx].z()) < 24 && abs(vtx[ivtx].position().Rho()) < 2;
			// it should be equivalent to patUtils procedure
			// only they use reverse: ! > 24 etc -- but without the >=, thus there is 1 bit of discrepancy
			if(its_good)
				{
				if(nGoodPV==0) goodPV=vtx[ivtx];
				nGoodPV++;
				}
			}

		// ----------------------------------------- Apply pileup reweighting
		// why don't use nGoodPV for Pile-Up?
		unsigned int num_inters = 0, num_inters_raw = 0;
		double weight_pu_test = weights_FULL[SYS_NOMINAL];
		// tests v7-10+ pileup is back
		if(isMC)
			{
			int ngenITpu = 0;
			fwlite::Handle < std::vector < PileupSummaryInfo > >puInfoH;
			puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
			if (!puInfoH.isValid())
				{
				puInfoH.getByLabel( ev, "addPileupInfo" );
				if (!puInfoH.isValid()) {printf("collection PileupSummaryInfo with name slimmedAddPileupInfo or addPileupInfo does not exist\n"); exit(0);}
				}
			// so here we have valid puInfoH
			// otherwise exit was called
			for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
				{
				//if (it->getBunchCrossing () == 0) ngenITpu += it->getPU_NumInteractions ();
				// guys and Mara use getTrueNumInteractions :
				if (it->getBunchCrossing () == 0) ngenITpu += it->getTrueNumInteractions();
				}

			//ngenITpu = nGoodPV; // based on nvtx
			//weight_PU = LumiWeights->weight (ngenITpu) * PUNorm[0];
			// So, in Pietro's approach ngenITpu is number of vertices in the beam crossing?
			//weight_PU = direct_pileup_reweight[ngenITpu];
			// Mara does:
			//num_inters = puInfoH->at(0).getTrueNumInteractions(); // in 76 it seems to not work, returns 0 always
			// Using Pietro's PU number vertices:
			num_inters = ngenITpu;
			if (num_inters<100) {
				weight_PU = direct_pileup_reweight[num_inters];
				weight_PU_up = direct_pileup_reweight_up[num_inters];
				weight_PU_down = direct_pileup_reweight_down[num_inters];
			}
			else {//weight_PU = 0; weight_PU_up = 0; weight_PU_down = 0;
				continue; // just move on
			}

			// FIXME: testing raw vtx.size()
			num_inters_raw = vtx.size();
			if (num_inters_raw<100) {weight_pu_test *= direct_pileup_reweight[num_inters_raw];}
			else {weight_pu_test *= 1.5e-16;}
			// TODO: implement error margins of pile-up
			}
		else
			{
			// get data pile-up into num_inters
			// for now use number of good vertices for the data
			// should substitute it with something more appropriate
			// num_inters = nGoodPV;
			// let's try using the size of primary vertex collection (before selecting the good vertices)
			num_inters = vtx.size();
			}

		fill_1d(string("weight_PU"), 200, 0., 2., weight_PU, 1);

		// applying PU weight
		for ( const auto s : weightSystematics )
			{
			if (s == SYS_PU_UP)
				weights_FULL[s] *= weight_PU_up;
			else if (s == SYS_PU_DOWN)
				weights_FULL[s] *= weight_PU_down;
			else
				weights_FULL[s] *= weight_PU;
			}
		// not applying for rawWeight

		//weights_FULL[SYS_NOMINAL] *= weight_PU;
		//weights_FULL[SYS_PU_UP]   *= weight_PU_up;
		//weights_FULL[SYS_PU_DOWN] *= weight_PU_down;

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights
		// sum_weights += weight;
		// sum_weights_raw += rawWeight;

		// int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight);

		// puweight
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   4, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   4, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 4, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 4, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 4, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   4, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   4, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 4, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 4, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 4, 1);

		// pu distrs
		fill_1d( string("pileup_beforetrig_num_inters_rawWeight"), 100, 0, 100, num_inters, rawWeight);
		fill_1d( string("pileup_beforetrig_num_inters_weight"),    100, 0, 100, num_inters, weights_FULL[SYS_NOMINAL]);

		// vtx.size
		fill_1d( string("pileup_beforetrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
		fill_1d( string("pileup_beforetrig_nvtx_weight"),    100, 0, 100,      vtx.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d( string("pileup_beforetrig_nvtx_weight_up"),    100, 0, 100,   vtx.size(), weights_FULL[SYS_PU_UP]);
		fill_1d( string("pileup_beforetrig_nvtx_weight_down"),    100, 0, 100, vtx.size(), weights_FULL[SYS_PU_DOWN]);

		// pu distrs
		fill_1d( string("pileup_beforetrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_beforetrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weights_FULL[SYS_NOMINAL]);

		// RHO distrs
		fill_1d( string("rho_beforetrig_rawWeight"), 100, 0, 100, rho, rawWeight);
		fill_1d( string("rho_beforetrig_weight"),         100, 0, 100, rho, weights_FULL[SYS_NOMINAL]);
		fill_1d( string("rho_beforetrig_weight_up"),      100, 0, 100, rho, weights_FULL[SYS_PU_UP]);
		fill_1d( string("rho_beforetrig_weight_down"),    100, 0, 100, rho, weights_FULL[SYS_PU_DOWN]);

		fill_1d( string("rhoCentral_beforetrig_rawWeight"), 100, 0, 100, rhoCentral, rawWeight);
		fill_1d( string("rhoCentral_beforetrig_weight"),    100, 0, 100, rhoCentral, weights_FULL[SYS_NOMINAL]);

		fill_1d( string("rhoCentralNeutral_beforetrig_rawWeight"), 100, 0, 100, rhoCentralNeutral, rawWeight);
		fill_1d( string("rhoCentralNeutral_beforetrig_weight"),    100, 0, 100, rhoCentralNeutral, weights_FULL[SYS_NOMINAL]);

		fill_1d( string("rhoCentralChargedPileUp_beforetrig_rawWeight"), 100, 0, 100, rhoCentralChargedPileUp, rawWeight);
		fill_1d( string("rhoCentralChargedPileUp_beforetrig_weight"),    100, 0, 100, rhoCentralChargedPileUp, weights_FULL[SYS_NOMINAL]);

		// -------------------------------------   Basic event selection

		// -------------------------------------------------- FIRST SECTION OF MC WEIGHT is over
		// MC weights (weights that don't change the integral of N events)

		// FISRT SECTION SUM
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   10, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   10, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 10, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 10, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 10, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   10, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   10, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 10, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 10, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 10, 1);





		// -------------------------------------------------- SECOND SECTION -- cuts on MET filters, trigger, lumi etc [11, 20]

		// ------------------------------------------------- Apply MET FILTERS
		/*
		 * MET filters are data-only thing -- remove events before passing and counting lumi, since MC is then normalized to data lumi
		 * thus after passing lumi data and MC should only have the same cuts
		 *
		 * info on MET filters and their presence in MINIAOD:
		 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
		 *   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
		 *     filter	location	data	MC(fullSim)	MC(fastSim)	comment
		 *     primary vertex filter	available in miniAOD v2	DONE	suggested	suggested	 
		 *     beam halo filter	available in miniAOD v2	DONE	suggested	not suggested	Beam Halo Presentation
		 *     HBHE noise filter	available in miniAOD v2	DONE	suggested	suggested	HCAL DPG Presentation
		 *     HBHEiso noise filter	available in miniAOD v2	DONE	suggested	suggested	same as above
		 *     ECAL TP filter	available in miniAOD v2	DONE	suggested	suggested	ECAL DPG Presentation
		 *     Bad PF Muon Filter	to be run on the fly	DONE	suggested	suggested	PPD presentation
		 *     Bad Charged Hadron Filter	to be run on the fly	DONE	suggested	suggested	PPD presentation
		 *     ee badSC noise filter	available in miniAOD v2	DONE	not suggested	not suggested
		 *   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
		 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes
		 *
		 *   https://twiki.cern.ch/twiki/bin/view/CMS/MissingET
		 *   and their hypernews:
		 *   https://hypernews.cern.ch/HyperNews/CMS/get/met.html
		 */

		edm::TriggerResultsByName metFilters = ev.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO
		if(!metFilters.isValid()){metFilters = ev.triggerResultsByName("RECO");} //if not present, then it's part of RECO
		if(!metFilters.isValid()){       
			printf("TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event\n");    
		}

		if (! isMC && metFilters.isValid()) {
			if(debug){
				//cout << "Printing PAT/RECO trigger list for MET filters here" << endl;
				//for(edm::TriggerNames::Strings::const_iterator trnames = metFilters.triggerNames().begin(); trnames!=metFilters.triggerNames().end(); ++trnames)
					//cout << *trnames << endl;
				cout << "----------- End of MET filters trigger list ----------" << endl;
				//return 0;
			}

			// event is good if all filters ar true
			bool filters1 = utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*", "Flag_EcalDeadCellTriggerPrimitiveFilter*");
			bool good_vertices = utils::passTriggerPatterns(metFilters, "Flag_goodVertices");
			bool eebad = utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter"); // "not suggested" for 2017 Moriond, but also they say Flag_eeBadScFilter TO BE USED
			bool halo  = utils::passTriggerPatterns(metFilters, "Flag_globalTightHalo2016Filter");
			// 2016 thing: bad muons
			bool flag_noBadMuons = utils::passTriggerPatterns(metFilters, "Flag_noBadMuons");
			//bool flag_duplicateMuons = utils::passTriggerPatterns(metFilters, "Flag_duplicateMuons");
			/* from
			 * https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#Event_flags
			 * Three flags are saved in the event:
			      Flag_badMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as bad
			      Flag_duplicateMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as duplicate
			      Flag_noBadMuons: the event does not contain any PF muon of pT > 20 GeV flagged as bad or duplicate (i.e. the event is safe)

			 * --- thus the Flag_noBadMuons should be enough
			 */

			// TODO: de-hardcode
			// recording which filters are _not_ passed (filtered events by the filter) (less computing, more to the point)
			if (!filters1           ) fill_1d(string("control_MET_filters"), 10, 0, 10, 0, 1);
			if (!good_vertices      ) fill_1d(string("control_MET_filters"), 10, 0, 10, 1, 1);
			if (!eebad              ) fill_1d(string("control_MET_filters"), 10, 0, 10, 2, 1);
			if (!halo               ) fill_1d(string("control_MET_filters"), 10, 0, 10, 3, 1);
			if (!flag_noBadMuons    ) fill_1d(string("control_MET_filters"), 10, 0, 10, 4, 1);
			//if (!flag_duplicateMuons) fill_1d(string("control_MET_filters"), 10, 0, 10, 5, 1);

			if (! (filters1 & good_vertices & eebad & halo & flag_noBadMuons)) continue;
			// these Flag_noBadMuons/Flag_duplicateMuons are MET flags (the issue with bad muons in 2016),
			// they are true if the MET got corrected and event is fine

			/* 
			 * add: BadChHadron and BadPFMuon -- it seems their name should be Flag_BadChHadron etc
			 *
			 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
			 * Bad PF Muon Filter	to be run on the fly	
			 * -- "run on the fly", no this flag in the data itself
			 *
			 * but at the same time:
			 *
			 * Note that with the in the re-miniaod you will have (will rerun as pointed out below for) the following flags for the "bad muon" events:
			      Bad PF Muon Filter
			      Bad Charged Hadrons
			      Flag_badMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)
			      Flag_duplicateMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)

			 * aha https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
			 * "Note that many of the current recommended filters can be accessed directly from Miniaod using the flag stored in the TriggerResults,
			 *  with the exception of Bad Charged Hadron and Bad Muon Filters."
			 * --- so, 2 vs 1 that there should be no Flags for these two in MINIAOD
			 *  they should be run on the fly
			 */


			/*
			 * MET POG gives some names to their filters instead of givin the name in code
			 * apparently the actual name in the code can be found at:
			 * https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
			 *
			 * and there is no BadChHandron
			 * the closes to their names are:
			 * BadChargedCandidateFilter BadPFMuonFilter
			 *
			 * -- need to print out what actually is in 03Feb ReReco & ask on hypernews.
			 *
			 *  found these:
			 *  root [7] metFilters.triggerNames()
			 *  (const std::vector<std::string> &)
			 *  { "Flag_duplicateMuons", "Flag_badMuons", "Flag_noBadMuons",
			 *    "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
			 *    "Flag_CSCTightHaloFilter", "Flag_CSCTightHaloTrkMuUnvetoFilter", "Flag_CSCTightHalo2015Filter",
			 *    "Flag_globalTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter",
			 *    "Flag_HcalStripHaloFilter", "Flag_hcalLaserEventFilter",
			 *    "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellBoundaryEnergyFilter",
			 *    "Flag_goodVertices",
			 *    "Flag_eeBadScFilter",
			 *    "Flag_ecalLaserCorrFilter",
			 *    "Flag_trkPOGFilters",
			 *    "Flag_chargedHadronTrackResolutionFilter",
			 *    "Flag_muonBadTrackFilter",
			 *    "Flag_trkPOG_manystripclus53X", "Flag_trkPOG_toomanystripclus53X", "Flag_trkPOG_logErrorTooManyClusters",
			 *    "Flag_METFilters" }
			 */
			}

		if(debug)
			cout << "met filters applied here" << endl;

		// passmetfilters
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   11, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   11, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 11, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 11, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 11, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   11, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   11, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 11, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 11, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 11, 1);



		// -------------------------------------------------- Skip bad lumi
		// 80X, v2
		if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock())) continue; 
		// Notice: it is the first continue in the event loop
		// n_events_pass_lumi += 1;
		// there is no sum_weights_pass_lumi -- lumi is for data only..

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		// increment( string("weightflow_weight_passed_lumi"), weight ); // should not matter
		// increment( string("weightflow_weight_up_passed_lumi"), weight_up ); // should not matter
		// increment( string("weightflow_weight_down_passed_lumi"), weight_down ); // should not matter

		// passlumi
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   12, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   12, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 12, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 12, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 12, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   12, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   12, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 12, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 12, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 12, 1);

		// --------------------------------------------- HLT TRIGGER
		// ---------------- and require compatibilitiy of the event with the PD

		// HLT2 was a quirk of Spring16 MC campaigns (noHLT/reHLT/withHLT thing)
		edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT2");
		if (!tr.isValid ()){
			if(debug){
				cout << "HLT2 is NOT valid, switching to HLT!\n";
				}
			tr = ev.triggerResultsByName ("HLT");
			if (!tr.isValid ()){
				cout << "Trigger HLT is not valid, exiting" << endl;
				return false;
				}
			else if(debug){
				cout << "Trigger HLT is valid\n";
				}
			}
		else if(debug){
			cout << "Trigger HLT2 is valid\n";
			if(!tr.isValid()){
				cout << "And now trigger HLT2 is NOT valid\n";
				}
			}

		if(debug){
			cout << "Printing HLT trigger list" << endl;
			//cout << "-- Commented out --" << endl;
			//edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");
			//tr = ev.triggerResultsByName ("HLT");
			//if (!tr.isValid ()){
			//for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames)
				//cout << *trnames << endl;
				//cout << "1: " << *trnames << endl;
				//}
			/*
			else
				{
				cout << "HLT is not valid\n";
				}
			cout << "Printing HLT2 trigger list" << endl;
			tr = ev.triggerResultsByName ("HLT2");
			if (!tr.isValid ()){
				for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames)
					cout << "2: " << *trnames << endl;
				}
			else
				{
				cout << "HLT2 is not valid\n";
				tr = ev.triggerResultsByName ("HLT");
				}
			*/
			cout << "----------- End of trigger list ----------" << endl;
			//return 0;
			}

		// Need either to simulate the HLT (https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#How_to_easily_emulate_HLT_paths) to match triggers.
		// Mara's triggers: HLT_Ele23_WPLoose_Gsf for electrons
		//                  HLT_IsoMu20 or HLT_IsoTkMu20 for muons
		//HLT_Iso(Tk)Mu22_v3
		//HLT_Ele27_WPTight_Gsf_v2
		// wildcard * for data
		// specific version for MC
		// Setting no-reHLT triggers!
		bool muTrigger = false;

		if (tr.isValid())
			{
			muTrigger =   ( isMC ?
			utils::passTriggerPatterns (tr, muHLT_MC1, muHLT_MC2) :
			utils::passTriggerPatterns (tr, muHLT_Data1, muHLT_Data2)
			);
			}
		else return 233;

		// if data and SingleElectron dataset and both triggers -- skip event
		// i.e. SingleElectron + HLT mu are removed
		/* no need to orthogonalize -- using only SingleMuon trigger
		if (!debug) {
			if (!isMC && isSingleElectronDataset && eTrigger && muTrigger) continue;
		}
		*/
		if (! muTrigger) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

		// TODO: ----------------------------- HLT efficiency scale factors
		// one should run it on the fired trigger objects,
		// I run it on selection candidates now
		// which is done below, when the candidates are selected

		// double HLT_efficiency_sf = 1.0;

		//HLT_efficiency_sf *= eTrigger  ? eHLT_sf[] : 1 ;
		//HLT_efficiency_sf *= muTrigger ? muHLT_SF[] : 1 ;

		// passtrig 6
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   13, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   13, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 13, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 13, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 13, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   13, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   13, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 13, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 13, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 13, 1);
		// increment( string("eventflow_1_passed_trig"), weight ); // should not matter
		// increment( string("weightflow_weight_up_passed_trig"), weight_up ); // should not matter
		// increment( string("weightflow_weight_down_passed_trig"), weight_down ); // should not matter

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

		// pu distrs
		fill_1d( string("pileup_passtrig_num_inters_rawWeight"), 100, 0, 100, num_inters, rawWeight);
		fill_1d( string("pileup_passtrig_num_inters_weight"),    100, 0, 100, num_inters, weights_FULL[SYS_NOMINAL]);

		// vtx.size
		fill_1d( string("pileup_passtrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
		fill_1d( string("pileup_passtrig_nvtx_weight"),    100, 0, 100, vtx.size(), weights_FULL[SYS_NOMINAL]);

		// pu distrs
		fill_1d( string("pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weights_FULL[SYS_NOMINAL]);

		// RHO distributions
		fill_1d( string("rho_passtrig_rawWeight"), 100, 0, 100, rho, rawWeight);
		fill_1d( string("rho_passtrig_weight"),    100, 0, 100, rho, weights_FULL[SYS_NOMINAL]);

		fill_1d( string("rhoCentral_passtrig_rawWeight"), 100, 0, 100, rhoCentral, rawWeight);
		fill_1d( string("rhoCentral_passtrig_weight"),    100, 0, 100, rhoCentral, weights_FULL[SYS_NOMINAL]);

		fill_1d( string("rhoCentralNeutral_passtrig_rawWeight"), 100, 0, 100, rhoCentralNeutral, rawWeight);
		fill_1d( string("rhoCentralNeutral_passtrig_weight"),    100, 0, 100, rhoCentralNeutral, weights_FULL[SYS_NOMINAL]);

		fill_1d( string("rhoCentralChargedPileUp_passtrig_rawWeight"), 100, 0, 100, rhoCentralChargedPileUp, rawWeight);
		fill_1d( string("rhoCentralChargedPileUp_passtrig_weight"),    100, 0, 100, rhoCentralChargedPileUp, weights_FULL[SYS_NOMINAL]);

		// fill_pu( string("pileup_passtrig_rawweight_pernuminters"), num_inters, rawWeight);
		// fill_pu( string("pileup_passtrig_weight_pernuminters"), num_inters, weight);

		// fill_pu( string("pileup_passtrig_rawweight_pergoodpv"), nGoodPV, rawWeight);
		// fill_pu( string("pileup_passtrig_weight_pergoodpv"), nGoodPV, weight);
		// fill_pu( string("pileup_passtrig_weight_up_pergoodpv"), nGoodPV, weight_up);
		// fill_pu( string("pileup_passtrig_weight_down_pergoodpv"), nGoodPV, weight_down);

		// fill_pu( string("pileup_passtrig_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
		// fill_pu( string("pileup_passtrig_weight_perrawvtxsize"), vtx.size(), weight_pu_test);


		// -------------------------------------------------- SECOND SECTION OF event cuts is over
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   20, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   20, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 20, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 20, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 20, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   20, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   20, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 20, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 20, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 20, 1);

		// ------------------------- event physics and the corresponding selection

		//------------------------- PROCESS OBJECTS
		// ----------------------------------------- possible THIRD SECTION of MC WEIGHTS and corrections (with efficiency SFs used as event probability correction)


		// ------------------------------------ actual particles

		pat::MuonCollection muons;
		fwlite::Handle<pat::MuonCollection> muonsHandle;
		muonsHandle.getByLabel(ev, "slimmedMuons");
		if(muonsHandle.isValid() ) muons = *muonsHandle;

		pat::ElectronCollection electrons;
		fwlite::Handle<pat::ElectronCollection> electronsHandle;
		electronsHandle.getByLabel(ev, "slimmedElectrons");
		if(electronsHandle.isValid() ) electrons = *electronsHandle;

		pat::JetCollection jets;
		fwlite::Handle<pat::JetCollection>jetsHandle;
		jetsHandle.getByLabel(ev, "slimmedJets");
		if(jetsHandle.isValid() ) jets = *jetsHandle;


		std::vector<reco::GenJet> genJets;
		fwlite::Handle<std::vector<reco::GenJet>>genJetsHandle;
		genJetsHandle.getByLabel(ev, "slimmedGenJets"); // twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#GenJets
		if(genJetsHandle.isValid() ) genJets = *genJetsHandle;

		// testing WNJets
		if(debug){
			cout << "N Jets:" << jets.size() << "\n";
			}

		/* not using photons now
		pat::PhotonCollection photons;
		fwlite::Handle<pat::PhotonCollection> photonsHandle;
		photonsHandle.getByLabel(ev, "slimmedPhotons");
		if(photonsHandle.isValid() ) photons = *photonsHandle;
		*/

		pat::METCollection mets;
		fwlite::Handle<pat::METCollection> metsHandle;
		if (isMC)
			metsHandle.getByLabel(ev, "slimmedMETs"); // 2016: slimmedMETs are METs corrected by muons
		else // ReReco 03Feb data
			metsHandle.getByLabel(ev, "slimmedMETsMuEGClean");
		// 2016: slimmedMETsMuEGClean are corrected by muons and electrons, only in Data!
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes
		if(metsHandle.isValid() ) mets = *metsHandle;
		pat::MET MET = mets[0];
		// LorentzVector met = mets[0].p4 ();

		fill_1d(string("control_met_main_pt"),  200, 0., 200., MET.pt(),  weights_FULL[SYS_NOMINAL]);
		fill_1d(string("control_met_main_phi"), 200, 0., 200., MET.phi(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("control_met_corr1_pt"),  200, 0., 200., MET.corPt(pat::MET::Type1),  weights_FULL[SYS_NOMINAL]);
		fill_1d(string("control_met_corr1_phi"), 200, 0., 200., MET.corPhi(pat::MET::Type1), weights_FULL[SYS_NOMINAL]);

		// testing WNJets
		if(debug){
			cout << "got main MET" << endl;
			}

		// also for control let's get uncorrected met and compare the two:
		if (!isMC)
			{
			pat::METCollection mets_uncorrected;
			fwlite::Handle<pat::METCollection> mets_uncorrectedHandle;
			mets_uncorrectedHandle.getByLabel(ev, "slimmedMETsUncorrected");
			if(mets_uncorrectedHandle.isValid() ) mets_uncorrected = *mets_uncorrectedHandle;
			pat::MET met_uncorrected = mets_uncorrected[0];
			fill_1d(string("control_met_slimmedMETsUncorrected_pt"), 200, 0., 200., met_uncorrected.pt(), weights_FULL[SYS_NOMINAL]);
			fill_1d(string("control_met_slimmedMETsUncorrected_diff_slimmedMETsMuEGClean_pt"), 200, -20., 20., met_uncorrected.pt()  - MET.pt(), weights_FULL[SYS_NOMINAL]);
			fill_1d(string("control_met_slimmedMETsUncorrected_diff_slimmedMETsMuEGClean_phi"),128, -3.2, 3.2, met_uncorrected.phi() - MET.phi(), weights_FULL[SYS_NOMINAL]);
			}

		if(debug){
			// MET try:
			double mypt = mets[0].shiftedPt(pat::MET::METUncertainty::JetEnUp);
			cout << "MET = " << mets[0].pt() << ", JetEnUp: " << mypt << endl;
			LorentzVector myshiftedMet = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
			cout << "MET = " << mets[0].pt() << ", JetEnUp: " << myshiftedMet.pt() << endl;
			}

		pat::TauCollection taus;
		fwlite::Handle<pat::TauCollection> tausHandle;
		tausHandle.getByLabel(ev, "slimmedTaus");
		if(tausHandle.isValid() ) taus = *tausHandle;



		//
		//
		// BELOW FOLLOWS THE PARTICLE SELECTION
		//
		//
		


		if(debug){
			cout << "got objects from the event, starting the analysis" << endl;
			}


		/* Tried photons for QCD issue in e-tau (excces from 2 very-weighted events)
		 * it didn't change anything
		pat::PhotonCollection selPhotons;
		//pat::PhotonCollection photons;
		for(unsigned int n=0; n<photons.size (); ++n)
			{
			pat::Photon& photon = photons[n];

			fill_2d(string("control_ph_slimmedphotons_pt_eta"), 250, 0., 500., 200, -3., 3., photon.pt(), photon.eta(), weight);
			fill_1d(string("control_ph_slimmedphotons_phi"), 128, -3.2, 3.2, photon.phi(), weight);

			if(photon.pt()<55)continue;
			if(fabs(photon.superCluster()->eta())>1.4442 ) continue;
			if(!patUtils::passId(photon, rho, patUtils::llvvPhotonId::Tight)) continue;

			fill_2d(string("control_ph_photonsIDed_pt_eta"), 250, 0., 500., 200, -3., 3., photon.pt(), photon.eta(), weight);
			fill_1d(string("control_ph_photonsIDed_phi"), 128, -3.2, 3.2, photon.phi(), weight);

			selPhotons.push_back(photon); 
			}
		*/

		//
		// LEPTON SELECTION
		//

		// ---------------------------------- ELECTRONS SELECTION
		/* int processElectrons_ID_ISO_Kinematics(pat::ElectronCollection& electrons, reco::Vertex goodPV, double rho, double weight, // input
		 *         patUtils::llvvElecId::ElecId el_ID, patUtils::llvvElecId::ElecId veto_el_ID,                           // config/cuts
		 *         patUtils::llvvElecIso::ElecIso el_ISO, patUtils::llvvElecIso::ElecIso veto_el_ISO,
		 *         double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
		 *         pat::ElectronCollection& selElectrons, LorentzVector& elDiff, unsigned int& nVetoE,                    // output
		 *         bool record, bool debug) // more output
		 */

		LorentzVector elDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton>
		pat::ElectronCollection selElectrons;
		unsigned int nVetoE(0);

		processElectrons_ID_ISO_Kinematics(electrons, goodPV, rho, weights_FULL[SYS_NOMINAL], patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
			kino_cuts_electron_pt, kino_cuts_electron_eta, kino_cuts_electron_veto_pt, kino_cuts_electron_veto_eta, selElectrons, elDiff, nVetoE, false, debug);

		if(debug){
			cout << "processed electrons" << endl;
			}

		// applying electron ID SFs
		if (isMC && selElectrons.size()>0)
			{
			double electron_sfs_weight = 1;
			for (int i = 0; i<selElectrons.size(); i++)
				{
				pat::Electron& el = selElectrons[i];
				double eta = el.eta();
				double pt  = el.pt();

				// for control (TODO: remove this later)
				double weight_reco, weight_id;

				// here X axis is eta, Y axis is pt
				// X is from -2.5 to 2.5 -- our eta is up to 2.4, should be ok
				//double bin_x = (pt < electron_effs_tracking_all_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_id_BCDEF_histo->GetXaxis()->GetXmax() - 1);
				double bin_x = eta;
				double bin_y = (pt < electron_effs_tracking_all_histo->GetYaxis()->GetXmax() ? pt : electron_effs_tracking_all_histo->GetYaxis()->GetXmax() - 1);
				weight_reco = electron_effs_tracking_all_histo->GetBinContent (electron_effs_tracking_all_histo->FindBin(bin_x, bin_y));

				//bin_x = eta;
				bin_y = (pt < electron_effs_id_all_histo->GetYaxis()->GetXmax() ? pt : electron_effs_id_all_histo->GetYaxis()->GetXmax() - 1);
				weight_id = electron_effs_id_all_histo->GetBinContent (electron_effs_id_all_histo->FindBin(bin_x, bin_y));

				fill_1d(string("weight_electron_effs_all_reco"),  200, 0., 1.1,   weight_reco, 1);
				fill_1d(string("weight_electron_effs_all_id"),    200, 0., 1.1,   weight_id,  1);
				electron_sfs_weight *= weight_reco * weight_id;
				}
			fill_1d(string("weight_electron_effs_all"),  200, 0., 1.1,   electron_sfs_weight, 1);

			// and apply to all weights:
			for ( const auto s : weightSystematics )
				{
				weights_FULL[s] *= electron_sfs_weight;
				}
			rawWeight *= electron_sfs_weight;
			//weights_FULL[SYS_NOMINAL] *= electron_sfs_weight * el_trig_weight;
			//weights_FULL[SYS_PU_UP]   *= electron_sfs_weight * el_trig_weight;
			//weights_FULL[SYS_PU_DOWN] *= electron_sfs_weight * el_trig_weight;
			}

		if(debug)
			cout << "applied electron SFs to event weight" << endl;



		// ---------------------------------- MUONS SELECTION
		LorentzVector muDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton> selLeptons;
		pat::MuonCollection selMuons;
		unsigned int nVetoMu(0);
		// unsigned int count_idiso_muons = 0;
		/*
		 * int processMuons_ID_ISO_Kinematics(pat::MuonCollection& muons, reco::Vertex goodPV,            // input
		 *         patUtils::llvvMuonId::MuonId mu_ID, patUtils::llvvMuonId::MuonId veto_mu_ID,       // config/cuts
		 *         patUtils::llvvMuonIso::MuonIso mu_ISO, patUtils::llvvMuonIso::MuonIso veto_mu_ISO,
		 *         double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
		 *         pat::MuonCOllection& selMuons, LorentzVector& muDiff, unsigned int& nVetoMu,       //output
		 *         bool record, bool debug) // more output
		 */
		processMuons_ID_ISO_Kinematics(muons, goodPV, weights_FULL[SYS_NOMINAL], patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,
			kino_cuts_muon_pt, kino_cuts_muon_eta, kino_cuts_muon_veto_pt, kino_cuts_muon_veto_eta, selMuons, muDiff, nVetoMu, false, debug);

		if(debug){
			cout << "processed muons" << endl;
			}

		// applying muon ID SFs
		if (isMC && selMuons.size()>0)
			{
			double muon_sfs_weight = 1;
			double bcdef_weight = 1;
			double gh_weight = 1;
			for (int i = 0; i<selMuons.size(); i++)
				{
				pat::Muon& mu = selMuons[i];
				double abs_eta = abs(mu.eta());
				double pt = mu.pt();

				// for control (TODO: remove this later)
				double bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso,
					gh_weight_trk, gh_weight_id, gh_weight_iso;

				// hopefully tracking won't overflow in eta range:
				bcdef_weight_trk = muon_effs_tracking_BCDEF_graph->Eval(abs_eta);
				// the id-s totally can overflow:
				double bin_x = (pt < muon_effs_id_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_id_BCDEF_histo->GetXaxis()->GetXmax() - 1);
				double bin_y = (abs_eta < muon_effs_id_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_id_BCDEF_histo->GetYaxis()->GetXmax() - 1);
				bcdef_weight_id = muon_effs_id_BCDEF_histo->GetBinContent (muon_effs_id_BCDEF_histo->FindBin(bin_x, bin_y));
				// these too:
				bin_x = (pt < muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax() - 1);
				bin_y = (abs_eta < muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() - 1);
				bcdef_weight_iso = muon_effs_iso_BCDEF_histo->GetBinContent (muon_effs_iso_BCDEF_histo->FindBin(bin_x, bin_y));

				fill_1d(string("weight_muon_effs_BCDEF_trk"),  200, 0., 1.1,   bcdef_weight_trk, 1);
				fill_1d(string("weight_muon_effs_BCDEF_id"),   200, 0., 1.1,   bcdef_weight_id,  1);
				fill_1d(string("weight_muon_effs_BCDEF_iso"),  200, 0., 1.1,   bcdef_weight_iso, 1);
				bcdef_weight *= bcdef_weight_trk * bcdef_weight_id * bcdef_weight_iso;

				// same ... for GH era:
				gh_weight_trk = muon_effs_tracking_GH_graph->Eval(abs_eta);
				// the id-s totally can overflow:
				bin_x = (pt < muon_effs_id_GH_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_id_GH_histo->GetXaxis()->GetXmax() - 1);
				bin_y = (abs_eta < muon_effs_id_GH_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_id_GH_histo->GetYaxis()->GetXmax() - 1);
				gh_weight_id = muon_effs_id_GH_histo->GetBinContent (muon_effs_id_GH_histo->FindBin(bin_x, bin_y));
				// these too:
				bin_x = (pt < muon_effs_iso_GH_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_iso_GH_histo->GetXaxis()->GetXmax() - 1);
				bin_y = (abs_eta < muon_effs_iso_GH_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_iso_GH_histo->GetYaxis()->GetXmax() - 1);
				gh_weight_iso = muon_effs_iso_GH_histo->GetBinContent (muon_effs_iso_GH_histo->FindBin(bin_x, bin_y));

				fill_1d(string("weight_muon_effs_GH_trk"),  200, 0., 1.1,   gh_weight_trk, 1);
				fill_1d(string("weight_muon_effs_GH_id"),   200, 0., 1.1,   gh_weight_id,  1);
				fill_1d(string("weight_muon_effs_GH_iso"),  200, 0., 1.1,   gh_weight_iso, 1);
				gh_weight *= gh_weight_trk * gh_weight_id * gh_weight_iso;
				}
			// and merge them:
			muon_sfs_weight = SingleMuon_data_bcdef_fraction * bcdef_weight + SingleMuon_data_gh_fraction * gh_weight;
			fill_1d(string("weight_muon_effs_BCDEF"), 200, 0., 1.1,   bcdef_weight, 1);
			fill_1d(string("weight_muon_effs_GH"),    200, 0., 1.1,   gh_weight, 1);
			fill_1d(string("weight_muon_effs_FULL"),  200, 0., 1.1,   muon_sfs_weight, 1);

			// and apply to all weights:
			for ( const auto s : weightSystematics )
				{
				weights_FULL[s] *= muon_sfs_weight;
				}
			rawWeight *= muon_sfs_weight;
			}

		if(debug)
			cout << "applied muon SFs to event weight" << endl;

		// Finally, merge leptons for cross-cleaning with taus and jets, and other conveniences:

		std::vector<patUtils::GenericLepton> selLeptons;
		for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
		for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
		std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);


		// -------------------- weights for trigger SFs
		if (isMC)
			{
			// calculate no-trig probability, then multiply them and get 1 - multiplication
			// should work for all cases
			// single trig
			// two leptons of the same type or different
			double lep_trig_weight = 1;

			/* no el trigger anymore
			double no_ele_trig = 1;
			if (eTrigger)
				{
				// calculate it the inverse-probbility way
				for (int i = 0; i<selElectrons.size(); i++)
                                	{
                                	pat::Electron& el = selElectrons[i];
					double eta = el.superCluster()->position().eta();
					double pt = el.pt();
					// here X axis is pt, Y axis is eta (from -2.5 to 2.5)
					double bin_x = (pt  < electron_effs_trg_all_histo->GetXaxis()->GetXmax() ? pt  : electron_effs_trg_all_histo->GetXaxis()->GetXmax() - 1);
					double bin_y = (eta < electron_effs_trg_all_histo->GetYaxis()->GetXmax() ? eta : electron_effs_trg_all_histo->GetYaxis()->GetXmax() - 1);
					no_ele_trig *= 1 - electron_effs_trg_all_histo->GetBinContent( electron_effs_trg_all_histo->FindBin(bin_x, bin_y) );
					}
				//el_trig_weight = 1 - no_trig; // so for 1 lepton it will = to the SF, for 2 there will be a mix
				fill_1d(string("weight_trigger_no_electron"),  200, 0., 1.1,   no_ele_trig, 1);
				}
			*/

			// also, if there was muon trigger apply the trigger weight:
			double no_mu_trig = 1;
			//double mu_trig_weight = 1;
			if (muTrigger)
				{
				// calculate it the inverse-probbility way
				for (int i = 0; i<selMuons.size(); i++)
                                	{
                                	pat::Muon& mu = selMuons[i];
					double abs_eta = abs(mu.eta());
					double pt = mu.pt();
					double bin_x = (pt < muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax() - 1);
					double bin_y = (abs_eta < muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() - 1);
					no_mu_trig *= SingleMuon_data_bcdef_fraction * (1 - muon_effs_trg_BCDEF_histo->GetBinContent( muon_effs_trg_BCDEF_histo->FindBin(bin_x, bin_y) )) +
							SingleMuon_data_gh_fraction * (1 - muon_effs_trg_GH_histo->GetBinContent( muon_effs_trg_GH_histo->FindBin(bin_x, bin_y) ));
					}
				//mu_trig_weight = 1 - no_trig; // so for 1 muon it will = to the SF, for 2 there will be a mix
				fill_1d(string("weight_trigger_no_muon"),  200, 0., 1.1,   no_mu_trig, 1);
				}

			lep_trig_weight = 1 - no_mu_trig; //*no_ele_trig;

			fill_1d(string("weight_trigger_weight"),  200, 0., 1.1,   lep_trig_weight, 1);
			/*
			if (eTrigger && selElectrons.size() == 1 && ! muTrigger)
				fill_1d(string("weight_trigger_weight_el"),  200, 0., 1.1,   lep_trig_weight, 1);
			if (eTrigger && selElectrons.size() == 2 && ! muTrigger)
				fill_1d(string("weight_trigger_weight_elel"),200, 0., 1.1,   lep_trig_weight, 1);
			if (!eTrigger && muTrigger && selMuons.size() == 1)
				fill_1d(string("weight_trigger_weight_mu"),  200, 0., 1.1,   lep_trig_weight, 1);
			if (!eTrigger && muTrigger && selMuons.size() == 2)
				fill_1d(string("weight_trigger_weight_mumu"),200, 0., 1.1,   lep_trig_weight, 1);
			if (eTrigger && muTrigger)
				fill_1d(string("weight_trigger_weight_elmu"),  200, 0., 1.1,   lep_trig_weight, 1);
			*/

			for ( const auto s : weightSystematics )
				{
				weights_FULL[s] *= lep_trig_weight;
				}
			rawWeight *= lep_trig_weight;
			}

		// ------------------------------------- Propagate lepton energy scale to MET
		/* no lepton corrections propagation yet
		met.setP4(met.p4() - muDiff - elDiff); // TODO: note this also propagates to all MET uncertainties -- does it??
		met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
		*/


		// ------------------------------------------ TAUS SELECTION

		//int processTaus_ID_ISO(pat::TauCollection& taus, double weight, // input
		//	string& tauID_decayMode, string& tauID,               // config/cuts
		//	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
		//	pat::TauCollection& selTaus,                          // output
		//	bool record, bool debug) // more output

		pat::TauCollection IDtaus, selTaus, selTaus_JetTauFakeRate,
			NoIsoIDtaus,  selNoIsoIDtaus,
			IDVLoosetaus, selVLooseTaus,
			IDLoosetaus,  selLooseTaus,
			IDMediumtaus, selMediumTaus,
			IDTighttaus,  selTightTaus,
			IDVTighttaus, selVTightTaus;
		//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
		/* all tau IDs:
		 * byVLooseIsolationMVArun2v1DBoldDMwLT
		 * byLooseIsolationMVArun2v1DBoldDMwLT
		 * byMediumIsolationMVArun2v1DBoldDMwLT
		 * byTightIsolationMVArun2v1DBoldDMwLT
		 * byVTightIsolationMVArun2v1DBoldDMwLT
		 */

		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, debug);

		// taus with no isolation (pre-iso control)
		for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
			{
			pat::Tau& tau = taus[n];

			//if (record)
			//	{
			//	fill_2d(string("control_tau_slimmedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
			//	fill_1d(string("control_tau_slimmedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
			//	}

			if (tau.tauID(tau_decayMode) < 0.5) continue;
			//if (tau.tauID(tauID) < 0.5)           continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
			if (tau.tauID(tau_againstMuon) < 0.5)  continue; // Medium working point not usable. Available values: Loose, Tight
			if (tau.tauID(tau_againstElectron) < 0.5) continue;

			NoIsoIDtaus.push_back(tau);

			//if (record)
				{
				fill_1d(string("control_tau_NoIsoIDtaus_pt"),  250,  0., 500. , tau.pt(), weights_FULL[SYS_NOMINAL]);
				fill_1d(string("control_tau_NoIsoIDtaus_eta"), 200, -4.,   4. , tau.eta(), weights_FULL[SYS_NOMINAL]);
				fill_1d(string("control_tau_NoIidedtaus_phi"), 128, -3.2,  3.2, tau.phi(), weights_FULL[SYS_NOMINAL]);
				}
			}

		std::sort (NoIsoIDtaus.begin(), NoIsoIDtaus.end(), utils::sort_CandidatesByPt);

		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_VLoose_ID, tau_againstMuon, tau_againstElectron, IDVLoosetaus, false, debug);
		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_Loose_ID,  tau_againstMuon, tau_againstElectron, IDLoosetaus,  false, debug);
		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_Medium_ID, tau_againstMuon, tau_againstElectron, IDMediumtaus, false, debug);
		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_Tight_ID,  tau_againstMuon, tau_againstElectron, IDTighttaus,  false, debug);
		processTaus_ID_ISO(taus, weights_FULL[SYS_NOMINAL], tau_decayMode, tau_VTight_ID, tau_againstMuon, tau_againstElectron, IDVTighttaus, false, debug);

		if(debug){
			cout << "selected taus [individual]" << endl;
			}

		//int processTaus_Kinematics(pat::TauCollection& taus,          // input
		//	double weight,
		//	double pt_cut, double eta_cut,
		//	pat::TauCollection& selTaus,                          // output
		//	bool record, bool debug) // more output

		processTaus_Kinematics(IDtaus, weights_FULL[SYS_NOMINAL],      tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, debug);
		processTaus_Kinematics(IDtaus, weights_FULL[SYS_NOMINAL],      jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selTaus_JetTauFakeRate, false, debug);
		// all-ID taus for measurement of el-mu fake rates
		processTaus_Kinematics(NoIsoIDtaus,  weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selNoIsoIDtaus, false, debug);
		processTaus_Kinematics(IDVLoosetaus, weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selVLooseTaus, false, debug);
		processTaus_Kinematics(IDLoosetaus,  weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selLooseTaus, false, debug);
		processTaus_Kinematics(IDMediumtaus, weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selMediumTaus, false, debug);
		processTaus_Kinematics(IDTighttaus,  weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selTightTaus, false, debug);
		processTaus_Kinematics(IDVTighttaus, weights_FULL[SYS_NOMINAL], jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selVTightTaus, false, debug);

		// ------------------------------------------ select the taus cleaned from leptons


		//int crossClean_in_dR(pat::TauCollection& selTaus, std::vector<patUtils::GenericLepton>& leptons,
		//	float min_dR,
		//	pat::TauCollection& selTausNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		pat::TauCollection selTausNoLep, selTaus_JetTauFakeRate_NoLep,
			selNoIsoIDtausNoLep,
			selVLooseTausNoLep, selLooseTausNoLep, selMediumTausNoLep, selTightTausNoLep, selVTightTausNoLep;

		crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weights_FULL[SYS_NOMINAL], string("selTausNoLep"),        false, debug);
		crossClean_in_dR(selTaus_JetTauFakeRate, selLeptons, 0.4, selTaus_JetTauFakeRate_NoLep, weights_FULL[SYS_NOMINAL], string("selTaus_JetTauFakeRate_NoLep"), false, debug);

		crossClean_in_dR(selNoIsoIDtaus, selLeptons, 0.4, selNoIsoIDtausNoLep, weights_FULL[SYS_NOMINAL], string("selNoIsoIDtausNoLep"),   false, debug);
		crossClean_in_dR(selVLooseTaus,  selLeptons, 0.4, selVLooseTausNoLep,  weights_FULL[SYS_NOMINAL], string("selVLooseTausNoLep"),   false, debug);
		crossClean_in_dR(selLooseTaus,   selLeptons, 0.4, selLooseTausNoLep,   weights_FULL[SYS_NOMINAL], string("selLooseTausNoLep"),   false, debug);
		crossClean_in_dR(selMediumTaus,  selLeptons, 0.4, selMediumTausNoLep,  weights_FULL[SYS_NOMINAL], string("selMediumTausNoLep"),   false, debug);
		crossClean_in_dR(selTightTaus,   selLeptons, 0.4, selTightTausNoLep,   weights_FULL[SYS_NOMINAL], string("selTightTausNoLep"),   false, debug);
		crossClean_in_dR(selVTightTaus,  selLeptons, 0.4, selVTightTausNoLep,  weights_FULL[SYS_NOMINAL], string("selVTightTausNoLep"),   false, debug);

		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Measurement_in_Z_tautau_events
		// Medium MVA (no dR03) is 0.97 +- 0.05

		/*
		 * tau ID SF needs:
		 * the SF
		 * the up/down variation
		 */
		if (isMC && selTausNoLep.size() > 0) // 2016 data/MC tau ID efficiency (all discriminators, all pt and eta ranges) = 0.83 +- 0.06
			{
			/*
			 * this calculates probability to get at least 1 tau from the given event, if the tdu ID efficiency was shifted by SF
			 * usually there is only 1 tau in the event
			 * rarely two are found (SingleElectron Run2016G):
			 *             mean    dev
			 * slimmedTaus 1.14722 0.758783
			 *              ratio 2-tau events to 1-tau
			 * selTaus      0.00079
			 * selTausNoLep 0.00077
			 *
			 * -- so, let's assume there are only 1-tau events:
			 */

			/*
			double pre_weight_tauIDsf      = 1;
			double pre_weight_tauIDsf_up   = 1;
			double pre_weight_tauIDsf_down = 1;
			for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
				{
				//pre_weight_tauIDsf *= (1 - 0.83 + r3->Gaus(0, 0.06)); // gaussian +- 0.06
				// recommendations update (noticed 8-11-2016, last page update 2016-11-07)
				// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_ID_efficiency
				// in ReReco it's different
				// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Measurement_in_Z_tautau_events
				// previous, randomized SF
				//pre_weight_tauIDsf *= (1 - 0.97 + r3->Gaus(0, 0.05)); // gaussian +- 0.05
				// now the SF shifts are accounted as systematics
				pre_weight_tauIDsf      *= (1 -  0.97);
				pre_weight_tauIDsf_up   *= (1 - (0.97 + 0.05)); // 1.02
				pre_weight_tauIDsf_down *= (1 - (0.97 - 0.05));
				}
			*/
			weight_tauIDsf      *= tauIDsf;
			weight_tauIDsf_up   *= tauIDsf + tauIDsf_shift;
			weight_tauIDsf_down *= tauIDsf - tauIDsf_shift;
			}
			// TODO: should here be a normalization to all MC events?
		fill_1d(string("weight_tauIDsf"), 100, 0., 2.,   weight_tauIDsf, 1); // const value now

		weight_without_tauIDsf = weights_FULL[SYS_NOMINAL];
		// weight *= weight_tauIDsf;
		// ------- apply tau weight at single-lepton selection

		if(debug){
			cout << "processed taus" << " N selTausNoLep = " << selTausNoLep.size() << endl;
			cout << "ID SF = " << weight_tauIDsf << endl;
			}

		//
		// JET/MET SELECTION
		//

		if(debug) cout << "Now update Jet Energy Corrections" << endl;




		// ----------------------- JETS CORRECTIONS, JEC, JER
		// ----------------------- UPDATE JEC

		//int processJets_CorrectJES_SmearJERnJES_ID_ISO(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
		//	bool isMC, double weight,
		//	double rho, unsigned int nGoodPV,
		//	FactorizedJetCorrector *jesCor,
		//	JetCorrectionUncertainty *totalJESUnc,
		//	double dR_max, // for jet matching in jet corrections smearing for MC
		//	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
		//	jet_id   & jetID,
		//	pu_jet_id& jetPUID,
		//	bool with_PUID,
		//	//double pt_cut, double eta_cut,
		//	TRandom3 *r3,   // the randomizer for the smearing
		//	LorentzVector& full_jet_corr, pat::JetCollection& selJets,                          // output
		//	bool record, bool debug) // more output

		//int processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
		//	bool isMC, double weight,
		//	double rho, unsigned int nGoodPV,
		//	FactorizedJetCorrector *jesCor,
		//	JetCorrectionUncertainty *totalJESUnc,
		//	double dR_max, // for jet matching in jet corrections smearing for MC
		//	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
		//	jet_id   & jetID,
		//	pu_jet_id& jetPUID,
		//	bool with_PUID,
		//	//double pt_cut, double eta_cut,
		//	TRandom3 *r3,   // the randomizer for the smearing
		//	LorentzVector& full_jet_corr,
		//	map<systematic_shift, pat::JetCollection>& selJets,                          // output
		//	bool record, bool debug) // more output

		//LorentzVector full_jet_corr(0., 0., 0., 0.);
		map<systematic_shift, LorentzVector> full_jet_corr;
		//pat::JetCollection IDjets;
		map<systematic_shift, pat::JetCollection> IDjets;
		// it's filled with jetSystematics by processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics
		//string jetID("Loose");
		//string jetPUID("MediumPU");

		processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics(jets, genJets, isMC, weights_FULL[SYS_NOMINAL], rho, nGoodPV, jesCor, totalJESUnc, 0.4/2,
			jet_resolution_in_pt, jet_resolution_sf_per_eta, /*jet_m_systematic_variation,*/ jetID, jetPUID, /*with_PU*/ false, r3, full_jet_corr, IDjets, true, debug);


		fill_3d(string("control_jet_full_jet_corr_pX_pY_pZ"), 10, -100., 100., 10, -100., 100., 10, -100., 100.,  full_jet_corr[SYS_NOMINAL].X(), full_jet_corr[SYS_NOMINAL].Y(), full_jet_corr[SYS_NOMINAL].Z(), weights_FULL[SYS_NOMINAL]);
		// 1000 bins

		//fill_2d(string("control_jet_full_jet_corr_pX_pY"), 100, -50., 50., 100, -50., 50.,  full_jet_corr[SYS_NOMINAL].X(), full_jet_corr[SYS_NOMINAL].Y(), weights_FULL[SYS_NOMINAL]);
		//fill_1d(string("control_jet_full_jet_corr_pZ"),    100, -50., 50., full_jet_corr[SYS_NOMINAL].Z(), weights_FULL[SYS_NOMINAL]);
		// 10 000 and 100 bins

		//met.setP4(met.p4() - full_jet_corr[SYS_NOMINAL]); // just return the full correction and propagate in place

		// -------------------------------------------------- different mets with jet corrections and some nominal met corrections
		/*
		 * the original pat::MET is stored in object MET
		 * it corresponds to jet corrections at the production of data
		 * (and it hass all the corrections: jet/met etc; default operations return Type1 values MET.pt() == MET.corPt(pat::MET::Type1) -- and same for .p4())
		 *
		 * now the jets were recorrected and these changes should be propagated to METs
		 */
		map<systematic_shift, LorentzVector> systematic_mets;

		// the NOMINAL and jet systematic shifts
		for ( const auto s : jetSystematics )
			{
			const char* name = systematic_shift_names[s];
			systematic_mets[s] = MET.p4() - full_jet_corr[s];
			fill_2d(string("control_jet_full_jet_corr_pX_pY") + name, 100, -50., 50., 100, -50., 50.,  full_jet_corr[s].X(), full_jet_corr[s].Y(), weights_FULL[SYS_NOMINAL]);
			}

		// also the original MET contains some Unclustered Energy UP/DOWN variations
		// UnclusteredEnUp=10, UnclusteredEnDown=11,
		// need to up-port them into new NOMINAL MET
		// with this
		// LorentzVector shiftedP4(METUncertainty shift, METCorrectionLevel level=Type1)  const ;
		// so I want: shifted + (new_nominal - old_nominal) = shifted - nominal + new_nominal
		systematic_mets[SYS_MET_UNCLUSTERED_EN_UP]   = MET.shiftedP4(pat::MET::UnclusteredEnUp)   - MET.p4() + systematic_mets[SYS_NOMINAL];
		systematic_mets[SYS_MET_UNCLUSTERED_EN_DOWN] = MET.shiftedP4(pat::MET::UnclusteredEnDown) - MET.p4() + systematic_mets[SYS_NOMINAL];


		// ok, let's save them all
		for ( const auto s : jetSystematics )
			{
			const char* name = systematic_shift_names[s];
			fill_1d(string("control_met_slimmedMETsMuEGClean_fulljetcorrs_pt") + name, 100, 0., 200., systematic_mets[s].pt(), weights_FULL[SYS_NOMINAL]);
			//fill_2d(string("control_met_slimmedMETs_fulljetcorrs_pt"), 200, 0., 200., 200, -4., 4., met.pt(), met.eta(), weight);
			}
		fill_1d(string("control_met_slimmedMETsMuEGClean_fulljetcorrs_pt") + systematic_shift_names[SYS_MET_UNCLUSTERED_EN_UP],   100, 0., 200., systematic_mets[SYS_MET_UNCLUSTERED_EN_UP].pt(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("control_met_slimmedMETsMuEGClean_fulljetcorrs_pt") + systematic_shift_names[SYS_MET_UNCLUSTERED_EN_DOWN], 100, 0., 200., systematic_mets[SYS_MET_UNCLUSTERED_EN_DOWN].pt(), weights_FULL[SYS_NOMINAL]);

		//int processJets_Kinematics(pat::JetCollection& jets, // input
		//	//bool isMC,
		//	double weight,
		//	double pt_cut, double eta_cut,
		//	pat::JetCollection& selJets,                 // output
		//	bool record, bool debug) // more output

		//pat::JetCollection selJets;
		map<systematic_shift, pat::JetCollection> selJets;
		for ( const auto s : jetSystematics )
			{
			processJets_Kinematics(IDjets[s], /*bool isMC,*/ weights_FULL[SYS_NOMINAL], jettaufr_jet_kino_cuts_pt, jettaufr_jet_kino_cuts_eta, selJets[s], true, debug);
			}

		pat::JetCollection selJets_JetTauFakeRate; // for fake rates in dileptons
		processJets_Kinematics(IDjets[SYS_NOMINAL], /*bool isMC,*/ weights_FULL[SYS_NOMINAL], jettaufr_jet_kino_cuts_pt, jettaufr_jet_kino_cuts_eta, selJets_JetTauFakeRate, true, debug);

		// ---------------------------- Clean jet collections from selected leptons
		// TODO: add gamma-cleaning as well?

		//int crossClean_in_dR(pat::JetCollection& selJets, std::vector<patUtils::GenericLepton>& selLeptons,
		//	float min_dR,
		//	pat::JetCollection& selJetsNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		//pat::JetCollection selJetsNoLep;
		map<systematic_shift, pat::JetCollection> selJetsNoLep;
		for ( const auto s : jetSystematics )
			{
			crossClean_in_dR(selJets[s], selLeptons, 0.4, selJetsNoLep[s], weights_FULL[SYS_NOMINAL],
					string("selJetsNoLep_") + systematic_shift_names[s], false, debug);
			}

		pat::JetCollection selJets_JetTauFakeRate_NoLep;
		crossClean_in_dR(selJets_JetTauFakeRate, selLeptons, 0.4, selJets_JetTauFakeRate_NoLep, weights_FULL[SYS_NOMINAL], string("selJets_JetTauFakeRate_NoLep"), false, debug);




		// ---------------------------- Clean jet collection from selected taus

		//int crossClean_in_dR(pat::JetCollection& selJets, pat::TauCollection& selTaus,
		//	float min_dR,
		//	pat::JetCollection& selJetsOut, // output
		//	string control_name,
		//	bool record, bool debug) // more output
		//pat::JetCollection selJetsNoLepNoTau;
		map<systematic_shift, pat::JetCollection> selJetsNoLepNoTau;
		for ( const auto s : jetSystematics )
			{
			crossClean_in_dR(selJetsNoLep[s], selTausNoLep, 0.4, selJetsNoLepNoTau[s], weights_FULL[SYS_NOMINAL],
				string("selJetsNoLepNoTau_") + systematic_shift_names[s], false, debug);
			}

		if(debug){
			cout << "processed jets" << endl;
			}


		// --------------------------- B-TAGGED JETS BJETS B-JETS

		//int processBJets_BTag(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
		//	BTagCalibrationReader& btagCal,
		//	struct bTaggingEfficiencyHistograms& bEffs,
		//	string& b_tagger_label, float b_tag_WP,
		//	pat::JetCollection& selBJets,                          // output
		//	bool record, bool debug) // more output

		//int processBJets_BTag_for_sys(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
		//	BTagCalibrationReader& btagCal, // BTagSFUtil& btsfutil, old b-tag SF weighting, done with bEffs now
		//	struct bTaggingEfficiencyHistograms& bEffs,
		//	string& b_tagger_label, float b_tag_WP,
		//	pat::JetCollection& selBJets,                          // output
		//	string& btag_sys_point,
		//	bool record, bool record_control, bool debug) // more output

		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
		string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		float btag_WP = 0.8484; // medium

		//pat::JetCollection selBJets;
		// JET_SHIFT, BTAG_SHIFT
		// in multiselect I'll use NOMINAL BTAG_SHIFT for jet corrections
		// and record BTAG_SHIFT-s for NOMINAL jets
		map<systematic_shift, map<systematic_shift, pat::JetCollection>> selBJets;
		for ( const auto jet_sys : jetSystematics )
			{
			for ( const auto btag_sys : btagSystematics )
				{
				bool record_btag_processing = jet_sys == SYS_NOMINAL && btag_sys == SYS_NOMINAL;
				double temp_jet_sys_b_tagging_sf = 1;
				// TODO: maybe selJetsNoLep are more appropriate here?
				processBJets_BTag_for_sys(selJetsNoLep[jet_sys], isMC, weights_FULL[SYS_NOMINAL],
					// for nominal jets find the btag SF weights
					// for other jets -- temporary weight
					(jet_sys == SYS_NOMINAL ? weight_bTaggingSFs[btag_sys] : temp_jet_sys_b_tagging_sf),
					btagCal, bEffs, btagger_label, btag_WP,
					selBJets[jet_sys][btag_sys],
					btag_sys_points[btag_sys],
					record_btag_processing, record_btag_processing, debug);
				}
			}

		//weight *= weight_bTaggingSF;
		for ( const auto s : weightSystematics )
			{
			if (s == SYS_BTAG_UP)
				weights_FULL[s] *= weight_bTaggingSFs[SYS_BTAG_UP];
			else if (s == SYS_BTAG_DOWN)
				weights_FULL[s] *= weight_bTaggingSFs[SYS_BTAG_DOWN];
			else
				weights_FULL[s] *= weight_bTaggingSFs[SYS_NOMINAL];
			}
		rawWeight *= weight_bTaggingSFs[SYS_NOMINAL];
		fill_1d(string("weight_bTaggingSF"), 200, 0., 2.,   weight_bTaggingSFs[SYS_NOMINAL], 1.);
		fill_1d(string("weight_bTaggingSF_UP"), 200, 0., 2.,   weight_bTaggingSFs[SYS_BTAG_UP], 1.);
		fill_1d(string("weight_bTaggingSF_DOWN"), 200, 0., 2.,   weight_bTaggingSFs[SYS_BTAG_DOWN], 1.);

		if(debug){
			cout << "processed b-tagged jets" << endl;
			}


		// also for referense:
		pat::TauCollection selTausNoLepNoJet;
		crossClean_in_dR(selTausNoLep, selJetsNoLep[SYS_NOMINAL], 0.4, selTausNoLepNoJet, weights_FULL[SYS_NOMINAL], string("selTausNoLepNoJet"), false, debug);


		// -------------------------------------------------- all particles are selected

		// -------------------------------------------------- cCONTROL VALUES

		// fill_1d(string("weightflow_mu_passmetfilters"), 300, 0, 300,   7, weight);

		fill_1d(string("n_slimmedjets"),  10, 0, 10,   jets.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selJets"),      10, 0, 10,   selJets[SYS_NOMINAL].size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selJetsNoLep"), 10, 0, 10,   selJetsNoLep[SYS_NOMINAL].size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selJetsNoLepNoTau"), 10, 0, 10,   selJetsNoLepNoTau[SYS_NOMINAL].size(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("n_selBJets"), 10, 0, 10,   selBJets[SYS_NOMINAL][SYS_NOMINAL].size(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("n_slimmedtaus"),  10, 0, 10,   taus.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTaus"),      10, 0, 10,   selTaus.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTausNoLep"), 10, 0, 10,   selTausNoLep.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTausNoLepNoJet"), 10, 0, 10,   selTausNoLepNoJet.size(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("n_slimmedtaus"),  10, 0, 10,   taus.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTaus"),      10, 0, 10,   selTaus.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTausNoLep"), 10, 0, 10,   selTausNoLep.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selTausNoLepNoJet"), 10, 0, 10,   selTausNoLepNoJet.size(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("n_slimmedmuons"),  10, 0, 10,   muons.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selMuons"),      10, 0, 10,   selMuons.size(), weights_FULL[SYS_NOMINAL]);

		fill_1d(string("n_slimmedelectrons"),  10, 0, 10,   electrons.size(), weights_FULL[SYS_NOMINAL]);
		fill_1d(string("n_selElectrons"),      10, 0, 10,   selElectrons.size(), weights_FULL[SYS_NOMINAL]);





		// -------------------------------------------------- THIRD SECTION OF mc weights for event probability SF-s
		for ( const auto s : weightSystematics )
			{
			double weight = (weights_FULL.find(s) != weights_FULL.end() ? weights_FULL[s] : weights_FULL[SYS_NOMINAL]);
			const char* name = systematic_shift_names[s];
			fill_1d(string("weightflow_mu_") + name, 300, 0, 300,   30, weight);
			fill_1d(string("weightflow_el_") + name, 300, 0, 300,   30, weight);
			fill_1d(string("weightflow_elel_") + name, 300, 0, 300, 30, weight);
			fill_1d(string("weightflow_elmu_") + name, 300, 0, 300, 30, weight);
			fill_1d(string("weightflow_mumu_") + name, 300, 0, 300, 30, weight);
			}

		fill_1d(string("eventflow_mu"), 300, 0, 300,   30, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   30, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 30, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 30, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 30, 1);


		if(debug){
			cout << "all particle-objects are processed, checking channel selection" << endl;
			}

		//
		// -------------------------------------------------- ASSIGN CHANNEL
		//
		// -------------------------------------------------- channel weightflow
		bool 
			isSingleMu(false),
			isSingleE(false),
			isDoubleMu(false),
			isDoubleE(false),
			isEMu(false);
		// int multiChannel(0);


		// int slepId(0);
		// if(selLeptons.size()>0)
			// slepId=selLeptons[0].pdgId();

		// bool iso_lep = nVetoE==0 && nVetoMu==0 && selLeptons.size() == 1 && nGoodPV != 0; // 2^5
		//if(selLeptons.size()!=1 || nGoodPV==0) continue; // Veto requirement alredy applied during the event categoriziation
		bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
		// maybe it is worth trying not considering the nGoodPV?

		// TODO: test if trigger needed here at all
		// FIXME: TESTING new, trigger-less channel assignment
		isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && clean_lep_conditions;
		isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && clean_lep_conditions;

		// TODO: this kind of thing:
		// if (isSingleMu) fill_pt_e( string("met0_all_leptoncorr_jetcorr_singlemu_pt"), met.pt(), weight);
		// if (isSingleE)  fill_pt_e( string("met0_all_leptoncorr_jetcorr_singleel_pt"), met.pt(), weight);

		/* HLT efficiencies, redo them sometime
		if (isSingleE)
			{
			fill_pt_e( string("singleel_electrons_pt"), selElectrons[0].pt(), weight);
			fill_pt_e( string("met0_singleel_slimmed_pt"), met.pt(), weight);
			// TODO: no HLT efficiency SF for electron HLT yet:
			float el_eta = fabs(selElectrons[0].eta());
			float el_pt  = fabs(selElectrons[0].pt());

			// Double_t muon_HLTeff_SF1 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l1_eta, l1_pt) );
			// Double_t muon_HLTeff_SF2 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l2_eta, l2_pt) );
			Double_t electron_HLTeff_SF = 1;

			weight *= electron_HLTeff_SF;
			}
		if (isSingleMu)
			{
			fill_pt_e( string("singlemu_muons_pt"),     selMuons[0].pt(), weight);
			fill_pt_e( string("met0_singlemu_slimmed_pt"), met.pt(), weight);
			// muon_HLTeff_TH2F->FindBin
			float mu_eta = fabs(selMuons[0].eta());
			float mu_pt  = fabs(selMuons[0].pt());

			// Double_t muon_HLTeff_SF1 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l1_eta, l1_pt) );
			// Double_t muon_HLTeff_SF2 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l2_eta, l2_pt) );
			//Double_t muon_HLTeff_SF = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(mu_eta, mu_pt) );
			Double_t muon_HLTeff_SF = 1;
			weight *= muon_HLTeff_SF;
			}
		*/

		if(debug){
			cout << "assigned lepton channel" << endl;
			}

		// ------------------------------------- Propagate lepton energy scale to MET
		// Propagate now (v13)
		/* no lepton corrections propagation v13.1
		met.setP4(met.p4() - muDiff - elDiff); // TODO: note this also propagates to all MET uncertainties -- does it??
		met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
		*/

		//fill_pt_e( string("met0_all_leptoncorr_pt"), met.pt(), weight);

		//if (isSingleE)  fill_pt_e( string("met0_singleel_leptoncorr_pt"), met.pt(), weight);
		//if (isSingleMu) fill_pt_e( string("met0_singlemu_leptoncorr_pt"), met.pt(), weight);

		if(debug){
			cout << "propagated lepton corrections to met" << endl;
			}


		if(debug){
			cout << "merged selected leptons" << endl;
			}

		// TODO: more conditions for double-lepton channel? No veto leptons etc?
		// in progress...

		// TODO: add selTaus.size() == 0?
		if (selLeptons.size()==2 && clean_lep_conditions)
			{
			// so this should be a double-lepton channel
			// TODO: remove tau-selection below?
			int dilep_ids = selLeptons[0].pdgId() * selLeptons[1].pdgId();

			if (fabs(dilep_ids) == 121 )
				{
				isDoubleE = true;
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[0].pt(), weights_FULL[SYS_NOMINAL]);
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[1].pt(), weights_FULL[SYS_NOMINAL]);

				// float l1_eta = fabs(selElectrons[0].eta());
				// float l1_pt  = fabs(selElectrons[0].pt());
				// float l2_eta = fabs(selElectrons[1].eta());
				// float l2_pt  = fabs(selElectrons[1].pt());

				// Double_t muon_HLTeff_SF1 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l1_eta, l1_pt) );
				// Double_t muon_HLTeff_SF2 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l2_eta, l2_pt) );
				Double_t electron_HLTeff_SF1 = 1;
				Double_t electron_HLTeff_SF2 = 1;

				//weight *= 1 - (1 - electron_HLTeff_SF1)*(1 - electron_HLTeff_SF2);
				}
			else if (fabs(dilep_ids) == 169 )
				{
				isDoubleMu = true;
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[0].pt(), weights_FULL[SYS_NOMINAL]);
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[1].pt(), weights_FULL[SYS_NOMINAL]);

				// Double_t electron_HLTeff_SF = 1;
				float l1_eta = fabs(selMuons[0].eta());
				float l1_pt  = fabs(selMuons[0].pt());
				float l2_eta = fabs(selMuons[1].eta());
				float l2_pt  = fabs(selMuons[1].pt());

				//Double_t muon_HLTeff_SF1 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l1_eta, l1_pt) );
				//Double_t muon_HLTeff_SF2 = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(l2_eta, l2_pt) );
				Double_t muon_HLTeff_SF1 = 1;
				Double_t muon_HLTeff_SF2 = 1;
				//weight *= 1 - (1 - muon_HLTeff_SF1)*(1 - muon_HLTeff_SF2);
				}
			else
				{
				isEMu = true;
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[0].pt(), weights_FULL[SYS_NOMINAL]);
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[1].pt(), weights_FULL[SYS_NOMINAL]);

				float mu_eta = fabs(selMuons[0].eta());
				float mu_pt  = fabs(selMuons[0].pt());
				//float el_eta = fabs(selElectrons[0].eta());
				//float el_pt  = fabs(selElectrons[0].pt());

				Double_t electron_HLTeff_SF = 1;
				//Double_t muon_HLTeff_SF = muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(mu_eta, mu_pt) );
				Double_t muon_HLTeff_SF = 1;
				//weight *= 1 - (1 - electron_HLTeff_SF)*(1 - muon_HLTeff_SF);
				}
			}




		unsigned int n_leptons = selLeptons.size();
		// unsigned int n_taus = selTaus.size();
		unsigned int n_taus = selTausNoLep.size();
		// unsigned int n_taus = selTausNoLepNoJet.size(); // Try "reverse" tau-jet cleanning logic
		//unsigned int n_bjets = selBJets[SYS_NOMINAL].size();
		//unsigned int n_jets = selJets.size();
		//unsigned int n_jets = selJetsNoLepNoTau[SYS_NOMINAL].size();
		unsigned int n_jets = selJetsNoLep[SYS_NOMINAL].size(); // noLep jet as in jet fake-rate study
		// unsigned int n_bjets = selSingleLepBJets.size();
		unsigned int n_bjets = selBJets[SYS_NOMINAL][SYS_NOMINAL].size();


		// --------------------------------------------- DILEPTON CHANNELS

		if (isEMu) // reducing the computation isDoubleE || isDoubleMu || isEMu)
			{
			// nominal weight for all the control distr-s
			// without PU:
			/*
			double rawWeight = weights_FULL[SYS_NOMINAL];
			// apply PU and get the weight with PU
			for ( const auto s : weightSystematics )
				{
				if (s == SYS_PU_UP)
					weights_FULL[s] *= weight_PU_up;
				else if (s == SYS_PU_DOWN)
					weights_FULL[s] *= weight_PU_down;
				else
					weights_FULL[s] *= weight_PU;
				}
			*/
			// in principle it spoils the control of the MC integral...
			// let's do it well...

			// rawWeight should be calculated together with weights_FULL
			// nominal weight:
			double weight = weights_FULL[SYS_NOMINAL];

			// vtx.size
			fill_1d( string("pileup_dilepton_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
			fill_1d( string("pileup_dilepton_nvtx_weight"),    100, 0, 100, vtx.size(), weight);

			// pu distrs
			fill_1d( string("pileup_dilepton_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
			fill_1d( string("pileup_dilepton_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);

			// RHO distributions
			fill_1d( string("rho_dilepton_rawWeight"), 100, 0, 100, rho, rawWeight);
			fill_1d( string("rho_dilepton_weight"),    100, 0, 100, rho, weight);

			fill_1d( string("rhoCentral_dilepton_rawWeight"), 100, 0, 100, rhoCentral, rawWeight);
			fill_1d( string("rhoCentral_dilepton_weight"),    100, 0, 100, rhoCentral, weight);

			fill_1d( string("rhoCentralNeutral_dilepton_rawWeight"), 100, 0, 100, rhoCentralNeutral, rawWeight);
			fill_1d( string("rhoCentralNeutral_dilepton_weight"),    100, 0, 100, rhoCentralNeutral, weight);

			fill_1d( string("rhoCentralChargedPileUp_dilepton_rawWeight"), 100, 0, 100, rhoCentralChargedPileUp, rawWeight);
			fill_1d( string("rhoCentralChargedPileUp_dilepton_weight"),    100, 0, 100, rhoCentralChargedPileUp, weight);

			int dilId (1);
			// slepId(0);
			LorentzVector dileptonSystem (0, 0, 0, 0);
			if(selLeptons.size()>=2)
				{
				for (size_t ilep = 0; ilep < 2; ilep++)
					{
					dilId *= selLeptons[ilep].pdgId();
					int id(abs (selLeptons[ilep].pdgId()));
					dileptonSystem += selLeptons[ilep].p4();
					}
				}

			// Event selection booleans
			// or dilepton mass > 20GeV?
			// and only for ee/mumu?
			bool passMllVeto(isEMu ? dileptonSystem.mass()>20. : (fabs(dileptonSystem.mass()-91.)>15 && dileptonSystem.mass()>20. ) );
			// bool passJetSelection(selJets.size()>1);
			bool passJetSelection(n_jets>1);
			bool passMetSelection(systematic_mets[SYS_NOMINAL].pt()>40.);
			bool passOS(selLeptons[0].pdgId() * selLeptons[1].pdgId() < 0 );
			// bool passBtagsSelection(selBJets[SYS_NOMINAL].size()>1);
			bool passBtagsSelection(n_bjets>0);

			record_jets_fakerate_distrs(string("dilep_"), string("pass2leps"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

			if (passJetSelection)
				{
				record_jets_fakerate_distrs(string("dilep_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
				}

			if (passJetSelection && passBtagsSelection)
				{
				record_jets_fakerate_distrs(string("dilep_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
				}

			if (isDoubleE)
				{

				if (passJetSelection)
					{
					record_jets_fakerate_distrs(string("elel_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					// loose taus for fake-factor method: selLooseTausNoLep
					record_jets_fakerate_distrs(string("elel_"), string("passjets_looseTaus"), selJets_JetTauFakeRate_NoLep, selLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs(string("elel_"), string("passjets_tightTaus"), selJets_JetTauFakeRate_NoLep, selTightTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_large_bins(string("elel_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					fill_1d(string("elel_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  1, selLooseTausNoLep.size() * weight);
					fill_1d(string("elel_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  2, selTausNoLep.size() * weight);
					fill_1d(string("elel_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  3, selTightTausNoLep.size() * weight);
					}

				if (passJetSelection && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("elel_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					// pu distrs
					fill_1d( string("pileup_elel_passos_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
					fill_1d( string("pileup_elel_passos_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
					// nGoodPV = vtx.size() now

					record_jets_fakerate_distrs(string("elel_"), string("passbtagfinal"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					//fill_1d( string("elel_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("elel_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					//fill_1d( string("elel_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("elel_selection_njets"), 10, 0, 10, selJetsNoLepNoTau[SYS_NOMINAL].size(), weight);
					fill_1d( string("elel_selection_nbjets"), 10, 0, 10, selBJets[SYS_NOMINAL][SYS_NOMINAL].size(), weight);

					// pt
					fill_1d( string("elel_selection_el1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("elel_selection_el2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("elel_selection_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("elel_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][0].pt(), weight);
					fill_1d( string("elel_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][1].pt(), weight);
					fill_1d( string("elel_selection_bjet_pt"), 200, 0, 400, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].pt(), weight);
					// fill_1d( string("elel_selection_met_pt"), 200, 0, 400, n_systematic_mets[SYS_NOMINAL].pt(), weight);
					fill_1d( string("elel_selection_met_pt"), 200, 0, 400, systematic_mets[SYS_NOMINAL].pt(), weight);

					// energies
					fill_1d( string("elel_selection_el1_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("elel_selection_el2_energy"), 200, 0, 400, selLeptons[1].energy(), weight);
					// fill_1d( string("elel_selection_tau_energy"), 200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("elel_selection_jet1_energy"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][0].energy(), weight);
					fill_1d( string("elel_selection_jet2_energy"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][1].energy(), weight);
					fill_1d( string("elel_selection_bjet_energy"), 200, 0, 400, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].energy(), weight);
					// fill_1d( string("elel_selection_met_energy"), 200, 0, 400, n_systematic_mets[SYS_NOMINAL].energy(), weight);
					fill_1d( string("elel_selection_met_energy"), 200, 0, 400, systematic_mets[SYS_NOMINAL].energy(), weight);

					//etas:
					fill_1d( string("elel_selection_el1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("elel_selection_el2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("elel_selection_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("elel_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[SYS_NOMINAL][0].eta(), weight);
					fill_1d( string("elel_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[SYS_NOMINAL][1].eta(), weight);
					fill_1d( string("elel_selection_bjet_eta"), 300, -3, 3, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].eta(), weight);
					// fill_1d( string("elel_selection_met_eta"), 300, -3, 3, n_systematic_mets[SYS_NOMINAL].eta(), weight);
					fill_1d( string("elel_selection_met_eta"), 300, -3, 3, systematic_mets[SYS_NOMINAL].eta(), weight);

					weightflow_control_elel_selection += weight;
					fill_1d( string("weightflow_control_elel_selection"),  300, -3, 3, weight, weight);
					}
				}

			if (isDoubleMu)
				{
				if (passJetSelection)
					{
					record_jets_fakerate_distrs(string("mumu_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					// loose taus for fake-factor method: selLooseTausNoLep
					record_jets_fakerate_distrs(string("mumu_"), string("passjets_looseTaus"), selJets_JetTauFakeRate_NoLep, selLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs(string("mumu_"), string("passjets_tightTaus"), selJets_JetTauFakeRate_NoLep, selTightTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_large_bins(string("mumu_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					fill_1d(string("mumu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  1, selLooseTausNoLep.size() * weight);
					fill_1d(string("mumu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  2, selTausNoLep.size() * weight);
					fill_1d(string("mumu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  3, selTightTausNoLep.size() * weight);
					}

				if (passJetSelection && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("mumu_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					// pu distrs
					fill_1d( string("pileup_mumu_passos_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
					fill_1d( string("pileup_mumu_passos_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
					// nGoodPV = vtx.size() now

					record_jets_fakerate_distrs(string("mumu_"), string("passbtagfinal"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					//fill_1d( string("mumu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("mumu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					//fill_1d( string("mumu_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("mumu_selection_njets"), 10, 0, 10, selJetsNoLepNoTau[SYS_NOMINAL].size(), weight);
					fill_1d( string("mumu_selection_nbjets"), 10, 0, 10, selBJets[SYS_NOMINAL][SYS_NOMINAL].size(), weight);

					// pt
					fill_1d( string("mumu_selection_mu1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("mumu_selection_mu2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("mumu_selection_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("mumu_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][0].pt(), weight);
					fill_1d( string("mumu_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[SYS_NOMINAL][1].pt(), weight);
					fill_1d( string("mumu_selection_bjet_pt"), 200, 0, 400, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].pt(), weight);
					// fill_1d( string("mumu_selection_met_pt"), 200, 0, 400, n_systematic_mets[SYS_NOMINAL].pt(), weight);
					fill_1d( string("mumu_selection_met_pt"), 200, 0, 400, systematic_mets[SYS_NOMINAL].pt(), weight);

					// energies
					fill_1d( string("mumu_selection_mu1_energy"), 200, 0, 200, selLeptons[0].energy(), weight);
					fill_1d( string("mumu_selection_mu2_energy"), 200, 0, 200, selLeptons[1].energy(), weight);
					// fill_1d( string("mumu_selection_tau_energy"), 200, 0, 200, selTausNoLep[0].energy(), weight);
					fill_1d( string("mumu_selection_jet1_energy"), 200, 0, 200, selJetsNoLepNoTau[SYS_NOMINAL][0].energy(), weight);
					fill_1d( string("mumu_selection_jet2_energy"), 200, 0, 200, selJetsNoLepNoTau[SYS_NOMINAL][1].energy(), weight);
					fill_1d( string("mumu_selection_bjet_energy"), 200, 0, 200, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].energy(), weight);
					// fill_1d( string("mumu_selection_met_energy"), 200, 0, 200, n_systematic_mets[SYS_NOMINAL].energy(), weight);
					fill_1d( string("mumu_selection_met_energy"), 200, 0, 200, systematic_mets[SYS_NOMINAL].energy(), weight);

					//etas:
					fill_1d( string("mumu_selection_mu1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("mumu_selection_mu2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("mumu_selection_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("mumu_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[SYS_NOMINAL][0].eta(), weight);
					fill_1d( string("mumu_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[SYS_NOMINAL][1].eta(), weight);
					fill_1d( string("mumu_selection_bjet_eta"), 300, -3, 3, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].eta(), weight);
					// fill_1d( string("mumu_selection_met_eta"), 300, -3, 3, n_systematic_mets[SYS_NOMINAL].eta(), weight);
					fill_1d( string("mumu_selection_met_eta"), 300, -3, 3, systematic_mets[SYS_NOMINAL].eta(), weight);

					weightflow_control_mumu_selection += weight;
					fill_1d( string("weightflow_control_mumu_selection"),  300, -3, 3, weight, weight);
					}
				}

			if (isEMu)
				{

				if (passJetSelection)
					{ // THE ELMU FAKE RATES
					// loose taus for fake-factor method: selLooseTausNoLep
					/* record fake rates at JER, JES and PU systematics
					 */
					{ // NOMINAL
					fill_1d( string("elmu_passjets_lep1_pt_NOMINAL"), 100, 0, 300, selLeptons[0].pt(), weight);
					fill_1d( string("elmu_passjets_lep2_pt_NOMINAL"), 100, 0, 300, selLeptons[1].pt(), weight);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_NoIsoIDtaus"), selJetsNoLep[SYS_NOMINAL], selNoIsoIDtausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_vlooseTaus"), selJetsNoLep[SYS_NOMINAL], selVLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_looseTaus"),  selJetsNoLep[SYS_NOMINAL], selLooseTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_mediumTaus"), selJetsNoLep[SYS_NOMINAL], selMediumTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_tightTaus"),  selJetsNoLep[SYS_NOMINAL], selTightTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_vtightTaus"), selJetsNoLep[SYS_NOMINAL], selVTightTausNoLep, visible_gen_taus, weight, isMC);
					}
					{ // with Top pT
					double weight_with_top = weights_FULL[SYS_TOP_PT];
					fill_1d( string("elmu_passjets_lep1_pt_TOP_PT"), 100, 0, 300, selLeptons[0].pt(), weight_with_top);
					fill_1d( string("elmu_passjets_lep2_pt_TOP_PT"), 100, 0, 300, selLeptons[1].pt(), weight_with_top);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_NoIsoIDtaus"), selJetsNoLep[SYS_NOMINAL], selNoIsoIDtausNoLep, visible_gen_taus, weight_with_top, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_vlooseTaus"),  selJetsNoLep[SYS_NOMINAL], selVLooseTausNoLep, visible_gen_taus, weight_with_top, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_looseTaus"),   selJetsNoLep[SYS_NOMINAL], selLooseTausNoLep,  visible_gen_taus, weight_with_top, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_mediumTaus"),  selJetsNoLep[SYS_NOMINAL], selMediumTausNoLep, visible_gen_taus, weight_with_top, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_tightTaus"),   selJetsNoLep[SYS_NOMINAL], selTightTausNoLep,  visible_gen_taus, weight_with_top, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_TOP_PT_vtightTaus"),  selJetsNoLep[SYS_NOMINAL], selVTightTausNoLep, visible_gen_taus, weight_with_top, isMC);
					}
					{ // PU UP
					double weight_up = weights_FULL[SYS_PU_UP];
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_UP_vlooseTaus"), selJetsNoLep[SYS_NOMINAL], selVLooseTausNoLep, visible_gen_taus, weight_up, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_UP_looseTaus"),  selJetsNoLep[SYS_NOMINAL], selLooseTausNoLep,  visible_gen_taus, weight_up, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_UP_mediumTaus"), selJetsNoLep[SYS_NOMINAL], selMediumTausNoLep, visible_gen_taus, weight_up, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_UP_tightTaus"),  selJetsNoLep[SYS_NOMINAL], selTightTausNoLep,  visible_gen_taus, weight_up, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_UP_vtightTaus"), selJetsNoLep[SYS_NOMINAL], selVTightTausNoLep, visible_gen_taus, weight_up, isMC);

					fill_1d( string("elmu_passjets_lep1_pt_PU_UP"), 200, 0, 400, selLeptons[0].pt(), weight_up);
					fill_1d( string("elmu_passjets_lep2_pt_PU_UP"), 200, 0, 400, selLeptons[1].pt(), weight_up);
					fill_1d( string("elmu_passjets_jet1_pt_PU_UP"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][0].pt(), weight_up);
					fill_1d( string("elmu_passjets_jet2_pt_PU_UP"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][1].pt(), weight_up);
					}
					{ // PU DOWN
					double weight_down = weights_FULL[SYS_PU_DOWN];
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_DOWN_vlooseTaus"), selJetsNoLep[SYS_NOMINAL], selVLooseTausNoLep, visible_gen_taus, weight_down, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_DOWN_looseTaus"),  selJetsNoLep[SYS_NOMINAL], selLooseTausNoLep,  visible_gen_taus, weight_down, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_DOWN_mediumTaus"), selJetsNoLep[SYS_NOMINAL], selMediumTausNoLep, visible_gen_taus, weight_down, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_DOWN_tightTaus"),  selJetsNoLep[SYS_NOMINAL], selTightTausNoLep,  visible_gen_taus, weight_down, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_PU_DOWN_vtightTaus"), selJetsNoLep[SYS_NOMINAL], selVTightTausNoLep, visible_gen_taus, weight_down, isMC);

					fill_1d( string("elmu_passjets_lep1_pt_PU_DOWN"), 200, 0, 400, selLeptons[0].pt(), weight_down);
					fill_1d( string("elmu_passjets_lep2_pt_PU_DOWN"), 200, 0, 400, selLeptons[1].pt(), weight_down);
					fill_1d( string("elmu_passjets_jet1_pt_PU_DOWN"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][0].pt(), weight_down);
					fill_1d( string("elmu_passjets_jet2_pt_PU_DOWN"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][1].pt(), weight_down);
					}

					{ // JER UP and DOWN
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_UP_vlooseTaus"), selJetsNoLep[SYS_JER_UP], selVLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_UP_looseTaus"),  selJetsNoLep[SYS_JER_UP], selLooseTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_UP_mediumTaus"), selJetsNoLep[SYS_JER_UP], selMediumTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_UP_tightTaus"),  selJetsNoLep[SYS_JER_UP], selTightTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_UP_vtightTaus"), selJetsNoLep[SYS_JER_UP], selVTightTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_DOWN_vlooseTaus"), selJetsNoLep[SYS_JER_DOWN], selVLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_DOWN_looseTaus"),  selJetsNoLep[SYS_JER_DOWN], selLooseTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_DOWN_mediumTaus"), selJetsNoLep[SYS_JER_DOWN], selMediumTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_DOWN_tightTaus"),  selJetsNoLep[SYS_JER_DOWN], selTightTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JER_DOWN_vtightTaus"), selJetsNoLep[SYS_JER_DOWN], selVTightTausNoLep, visible_gen_taus, weight, isMC);
					}
					{ // JES UP and DOWN
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_UP_vlooseTaus"),   selJetsNoLep[SYS_JES_UP], selVLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_UP_looseTaus"),    selJetsNoLep[SYS_JES_UP], selLooseTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_UP_mediumTaus"),   selJetsNoLep[SYS_JES_UP], selMediumTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_UP_tightTaus"),    selJetsNoLep[SYS_JES_UP], selTightTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_UP_vtightTaus"),   selJetsNoLep[SYS_JES_UP], selVTightTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_DOWN_vlooseTaus"), selJetsNoLep[SYS_JES_DOWN], selVLooseTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_DOWN_looseTaus"),  selJetsNoLep[SYS_JES_DOWN], selLooseTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_DOWN_mediumTaus"), selJetsNoLep[SYS_JES_DOWN], selMediumTausNoLep, visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_DOWN_tightTaus"),  selJetsNoLep[SYS_JES_DOWN], selTightTausNoLep,  visible_gen_taus, weight, isMC);
					record_jets_fakerate_distrs_1D_2D(string("elmu_"), string("passjets_JES_DOWN_vtightTaus"), selJetsNoLep[SYS_JES_DOWN], selVTightTausNoLep, visible_gen_taus, weight, isMC);

					if (selJetsNoLep[SYS_JES_UP].size() > 0) fill_1d( string("elmu_passjets_jet1_pt_JES_UP"), 200, 0, 400, selJetsNoLep[SYS_JES_UP][0].pt(), weight);
					if (selJetsNoLep[SYS_JES_UP].size() > 1) fill_1d( string("elmu_passjets_jet2_pt_JES_UP"), 200, 0, 400, selJetsNoLep[SYS_JES_UP][1].pt(), weight);
					if (selJetsNoLep[SYS_JES_DOWN].size() > 0) fill_1d( string("elmu_passjets_jet1_pt_JES_DOWN"), 200, 0, 400, selJetsNoLep[SYS_JES_DOWN][0].pt(), weight);
					if (selJetsNoLep[SYS_JES_DOWN].size() > 1) fill_1d( string("elmu_passjets_jet2_pt_JES_DOWN"), 200, 0, 400, selJetsNoLep[SYS_JES_DOWN][1].pt(), weight);
					}

					// also large bins for (medium) jet fake rate
					record_jets_fakerate_distrs_large_bins(string("elmu_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					// met
					fill_1d( string("met_elmu_passjets"),    100, 0, 500, systematic_mets[SYS_NOMINAL].pt(), weight);
					//systematic_mets[SYS_NOMINAL];

					// pu distrs
					fill_1d( string("pileup_elmu_passjets_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
					fill_1d( string("pileup_elmu_passjets_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);

					// RHO distributions
					fill_1d( string("rho_elmu_passjets_rawWeight"), 100, 0, 100, rho, rawWeight);
					fill_1d( string("rho_elmu_passjets_weight"),    100, 0, 100, rho, weight);
					fill_1d( string("rhoCentral_elmu_passjets_rawWeight"), 100, 0, 100, rhoCentral, rawWeight);
					fill_1d( string("rhoCentral_elmu_passjets_weight"),    100, 0, 100, rhoCentral, weight);
					fill_1d( string("rhoCentralNeutral_elmu_passjets_rawWeight"), 100, 0, 100, rhoCentralNeutral, rawWeight);
					fill_1d( string("rhoCentralNeutral_elmu_passjets_weight"),    100, 0, 100, rhoCentralNeutral, weight);
					fill_1d( string("rhoCentralChargedPileUp_elmu_passjets_rawWeight"), 100, 0, 100, rhoCentralChargedPileUp, rawWeight);
					fill_1d( string("rhoCentralChargedPileUp_elmu_passjets_weight"),    100, 0, 100, rhoCentralChargedPileUp, weight);

					fill_1d(string("elmu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  1, selLooseTausNoLep.size() * weight);
					fill_1d(string("elmu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  2, selTausNoLep.size() * weight);
					fill_1d(string("elmu_passjets_fakefactor_n_loose_tight_taus"), 5, 0, 5,  3, selTightTausNoLep.size() * weight);

					// pt
					fill_1d( string("elmu_passjets_lep1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("elmu_passjets_lep2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("elmu_passjets_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("elmu_passjets_jet1_pt"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][0].pt(), weight);
					fill_1d( string("elmu_passjets_jet2_pt"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][1].pt(), weight);
					// fill_1d( string("elmu_passjets_met_pt"), 200, 0, 400, n_systematic_mets[SYS_NOMINAL].pt(), weight);
					fill_1d( string("elmu_passjets_met_pt"), 200, 0, 400, systematic_mets[SYS_NOMINAL].pt(), weight);

					// energies
					fill_1d( string("elmu_passjets_lep1_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("elmu_passjets_lep2_energy"), 200, 0, 400, selLeptons[1].energy(), weight);
					// fill_1d( string("elmu_passjets_tau_energy"), 200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("elmu_passjets_jet1_energy"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][0].energy(), weight);
					fill_1d( string("elmu_passjets_jet2_energy"), 200, 0, 400, selJetsNoLep[SYS_NOMINAL][1].energy(), weight);
					// fill_1d( string("elmu_passjets_met_energy"), 200, 0, 400, n_systematic_mets[SYS_NOMINAL].energy(), weight);
					fill_1d( string("elmu_passjets_met_energy"), 200, 0, 400, systematic_mets[SYS_NOMINAL].energy(), weight);

					//etas:
					fill_1d( string("elmu_passjets_lep1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("elmu_passjets_lep2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("elmu_passjets_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("elmu_passjets_jet1_eta"), 300, -3, 3, selJetsNoLep[SYS_NOMINAL][0].eta(), weight);
					fill_1d( string("elmu_passjets_jet2_eta"), 300, -3, 3, selJetsNoLep[SYS_NOMINAL][1].eta(), weight);
					// fill_1d( string("elmu_passjets_met_eta"), 300, -3, 3, n_systematic_mets[SYS_NOMINAL].eta(), weight);
					fill_1d( string("elmu_passjets_met_eta"), 300, -3, 3, systematic_mets[SYS_NOMINAL].eta(), weight);

					// the optional b-jets:
					if (selBJets[SYS_NOMINAL][SYS_NOMINAL].size()>0)
						{
						fill_1d( string("elmu_passjets_bjet_pt"), 200, 0, 400,     selBJets[SYS_NOMINAL][SYS_NOMINAL][0].pt(), weight);
						fill_1d( string("elmu_passjets_bjet_energy"), 200, 0, 400, selBJets[SYS_NOMINAL][SYS_NOMINAL][0].energy(), weight);
						fill_1d( string("elmu_passjets_bjet_eta"), 300, -3, 3,     selBJets[SYS_NOMINAL][SYS_NOMINAL][0].eta(), weight);
						}
					}

				}
			}

		// TODO: properly count multichannel?
		if (isSingleE && isSingleMu) nMultiChannel++;

		if(debug){
			cout << "end of event" << endl;
			}


		} // End single file event loop

	delete file;
	} // End loop on files

printf("Done processing the job of files\n");

printf("End of (file loop) the job.\n");

// Controls distributions of processed particles



if(nMultiChannel>0) cout << "Warning! There were " << nMultiChannel << " multi-channel events out of " << totalEntries << " events!" << endl;
printf ("\n");

//##############################################
//########    SAVING HISTOs TO FILE     ########
//##############################################
//save control plots to file
printf ("Results save in %s\n", outUrl.Data());


// CONTROL DISTRS, ROOT OUTPUT

for(std::map<string, std::map<string, TH1D>>::iterator it = th1d_distr_maps_control.begin(); it != th1d_distr_maps_control.end(); ++it)
	{
	// const std::pair <string,string> *key = &it->first;
	string channel = it->first;

	//outUrl.Data() is dtag_jobnum
	// use them separately, take from: dtag_s, job_num
	// TFile* out_f = TFile::Open (TString(outUrl.Data() + string("_") + channel + string(".root")), "CREATE");
	TFile* out_f = TFile::Open (outdir + TString(string("/") + dtag_s + string("_") + job_num + '.' + channel + string(".root")), "CREATE");
	// string mc_decay_suffix = key->first;
	std::map<string, TH1D> * th1d_controlpoints = & it->second;

	for(std::map<string, TH1D>::iterator it = th1d_controlpoints->begin(); it != th1d_controlpoints->end(); ++it)
		{
		string controlpoint_name = it->first;
		TH1D * distr = & it->second;
		distr->SetName(controlpoint_name.c_str());
		distr->Write();
		out_f->Write(controlpoint_name.c_str());
		//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
		}

	std::map<string, TH1I> * th1i_controlpoints = & th1i_distr_maps_control[channel];

	for(std::map<string, TH1I>::iterator it = th1i_controlpoints->begin(); it != th1i_controlpoints->end(); ++it)
		{
		string controlpoint_name = it->first;
		TH1I * distr = & it->second;
		distr->SetName(controlpoint_name.c_str());
		distr->Write();
		out_f->Write(controlpoint_name.c_str());
		//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
		}

	std::map<string, TH2D> * th2d_controlpoints = & th2d_distr_maps_control[channel];

	for(std::map<string, TH2D>::iterator it = th2d_controlpoints->begin(); it != th2d_controlpoints->end(); ++it)
		{
		string controlpoint_name = it->first;
		TH2D * distr = & it->second;
		distr->SetName(controlpoint_name.c_str());
		distr->Write();
		out_f->Write(controlpoint_name.c_str());
		//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
		}

	std::map<string, TH3D> * th3d_controlpoints = & th3d_distr_maps_control[channel];

	for(std::map<string, TH3D>::iterator it = th3d_controlpoints->begin(); it != th3d_controlpoints->end(); ++it)
		{
		string controlpoint_name = it->first;
		TH3D * distr = & it->second;
		distr->SetName(controlpoint_name.c_str());
		distr->Write();
		out_f->Write(controlpoint_name.c_str());
		//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
		}

	out_f->Close();
	}

printf ("New output results saved in %s\n", (outdir.Data() + string("/") + dtag_s + string("_") + job_num + string(".<channel>") + string(".root")).c_str());


// JOB_DONE file
// needed for precise accounting of done jobs, since a job can output many root files

FILE *csv_out;
string FileName = ((outUrl.ReplaceAll(".root",""))+".job_done").Data();
csv_out = fopen(FileName.c_str(), "w");

fprintf(csv_out, "The job is done!\n\n");
fprintf(csv_out, "# wheightflow control with plain doubles\n");
fprintf(csv_out, "job_num,channel,weight\n");
fprintf(csv_out, "%s,weightflow_control_el_selection,%g\n",   job_num.c_str(), weightflow_control_el_selection);
fprintf(csv_out, "%s,weightflow_control_mu_selection,%g\n",   job_num.c_str(), weightflow_control_mu_selection);
fprintf(csv_out, "%s,weightflow_control_elel_selection,%g\n", job_num.c_str(), weightflow_control_elel_selection);
fprintf(csv_out, "%s,weightflow_control_mumu_selection,%g\n", job_num.c_str(), weightflow_control_mumu_selection);
fprintf(csv_out, "%s,weightflow_control_elmu_selection,%g\n", job_num.c_str(), weightflow_control_elmu_selection);

// printout_counters(csv_out, job_def);
// printout_distrs(csv_out, job_def);

fclose(csv_out);


// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}


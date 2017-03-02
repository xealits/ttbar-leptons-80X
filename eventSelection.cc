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
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

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
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

#include <map>
#include <string>

//#include "jetDistrs.h"
#include "recordFuncs.h"

#include "ProcessingMuons.cc"
#include "ProcessingElectrons.cc"
#include "ProcessingTaus.cc"
#include "ProcessingJets.cc"
#include "ProcessingBJets.cc"
#include "ProcessingDRCleaning.cc"

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
 * string find_W_decay(const reco::Candidate * W) {
 *
 * returns a substring {el, mu, tau, q} identifing which decay happend
 * or returns "" string if something weird happend
 *
 * Usage:
 * const reco::Candidate * W = p.daughter( W_num );
 * mc_decay += find_W_decay(W);
*/
string find_W_decay(const reco::Candidate * W) {
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
			if (d0_id == 15 || d1_id == 15 ) return string("tau");
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


//jetToTauFakeRate(TH3F * tau_fake_rate_jets_histo, TH3F * tau_fake_rate_taus_histo, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]))
//jetToTauFakeRate(tau_fake_rate_jets_histo1, tau_fake_rate_taus_histo1, tau_fake_rate_jets_histo2, tau_fake_rate_taus_histo2, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]));
//jetToTauFakeRate(TH3F * tau_fake_rate_jets_histo1, TH3F * tau_fake_rate_taus_histo1, TH3F * tau_fake_rate_jets_histo2, TH3F * tau_fake_rate_taus_histo2, Double_t tau_fake_rate_histo1_fraction, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius)
// later tau_fake_rate_histo1_fraction can be a TH3F histogram with fractions per pt, eta, radius
double jetToTauFakeRate_vanila(TH3F * tau_fake_rate_jets_histo1, TH3F * tau_fake_rate_taus_histo1, TH3F * tau_fake_rate_jets_histo2, TH3F * tau_fake_rate_taus_histo2, Double_t tau_fake_rate_histo1_fraction, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	// the tau_fake_rate_jets_histo and tau_fake_rate_taus_histo
	// are identical TH3F histograms
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





/* New fake rates
 * via extension of existing distributions
 * assuming pt, eta and radius distributions are independent.
 */
typedef struct {
	TH1D* x;
	TH1D* y;
	TH1D* z;
} FakeRateProjections;

double jetToTauFakeRate_Projections(
		FakeRateProjections & tau_fake_rate_jets_histo,
		FakeRateProjections & tau_fake_rate_taus_histo,
		Double_t tau_fake_rate_histo1_fraction,
		Double_t jet_pt, Double_t jet_eta, Double_t jet_radius,
		bool debug)

	{
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
	double fakerate = (taus_ratex/jets_ratex) * (taus_ratey/jets_ratey) * (taus_ratez/jets_ratez);

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
int debug_len           = runProcess.getParameter<int>  ("debug_len");
bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
bool isMC            = runProcess.getParameter<bool>  ("isMC");
double xsec          = runProcess.getParameter<double>("xsec");
int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
TString dtag         = runProcess.getParameter<std::string>("dtag");
string dtag_s        = runProcess.getParameter<std::string>("dtag");
string job_num       = runProcess.getParameter<std::string>("job_num");

JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

TString outUrl = runProcess.getParameter<std::string>("outfile");
TString outdir = runProcess.getParameter<std::string>("outdir");

string  muHLT_MC1   = runProcess.getParameter<std::string>("muHLT_MC1"),   muHLT_MC2   = runProcess.getParameter<std::string>("muHLT_MC2"),
	muHLT_Data1 = runProcess.getParameter<std::string>("muHLT_Data1"), muHLT_Data2 = runProcess.getParameter<std::string>("muHLT_Data2"),
	elHLT_Data  = runProcess.getParameter<std::string>("elHLT_Data"),  elHLT_MC    = runProcess.getParameter<std::string>("elHLT_MC");

cout << "Triggers:" << endl;
cout << muHLT_MC1 << '\t' << muHLT_MC2 << '\t' << muHLT_Data1 << '\t' << muHLT_Data2 << endl;
cout << elHLT_Data << '\t' << elHLT_MC << endl;

// Kino cuts
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
/*
string tau_decayMode("decayModeFinding"),
	//tau_ID("byMediumCombinedIsolationDeltaBetaCorr3Hits"),
	tau_ID("byMediumIsolationMVArun2v1DBoldDMwLT"),
	tau_againstMuon("againstMuonTight3"),
	tau_againstElectron("againstElectronTightMVA6");
*/
string tau_decayMode = runProcess.getParameter<std::string>("tau_decayMode"),
	tau_ID       = runProcess.getParameter<std::string>("tau_ID"),
	tau_againstMuon     = runProcess.getParameter<std::string>("tau_againstMuon"),
	tau_againstElectron = runProcess.getParameter<std::string>("tau_againstElectron");

cout << "Tau IDs:" << tau_decayMode << '\t' << tau_ID << '\t' << tau_againstMuon << '\t' << tau_againstElectron << endl;


cout << "Output directory: " << outdir << endl;
	
const edm::ParameterSet& myVidElectronIdConf = runProcess.getParameterSet("electronidparas");
const edm::ParameterSet& myVidElectronMainIdWPConf = myVidElectronIdConf.getParameterSet("tight");
const edm::ParameterSet& myVidElectronVetoIdWPConf = myVidElectronIdConf.getParameterSet("loose");
	
VersionedPatElectronSelector electronVidMainId(myVidElectronMainIdWPConf);
VersionedPatElectronSelector electronVidVetoId(myVidElectronVetoIdWPConf);
	
TString suffix = runProcess.getParameter < std::string > ("suffix");
std::vector < std::string > urls = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
//TString baseDir = runProcess.getParameter < std::string > ("dirName");
//  if (mctruthmode != 0) //FIXME
//    {
//      outFileUrl += "_filt";
//      outFileUrl += mctruthmode;
//    }

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

/* not in use now
 * I need some method to check which jobs run now
TString outTxtUrl = outUrl + ".ran";
FILE *outTxtFile = NULL;
if (!isMC) outTxtFile = fopen (outTxtUrl.Data (), "w");
printf ("StartFile URL = %s\n", outTxtUrl.Data ());
if (outTxtFile) fclose (outTxtFile);
*/

//tree info
TString dirname = runProcess.getParameter < std::string > ("dirName");

//systematics
std::vector<TString> systVars(1,"");
if(runSystematics && isMC)
	{
	systVars.push_back("jerup" );     systVars.push_back("jerdown"   );
	systVars.push_back("jesup" );     systVars.push_back("jesdown"   );
	//systVars.push_back("lesup" );   systVars.push_back("lesdown"   );
	systVars.push_back("leffup");     systVars.push_back("leffdown"  );
	systVars.push_back("puup"  );     systVars.push_back("pudown"   );
	systVars.push_back("umetup");     systVars.push_back("umetdown" );
	systVars.push_back("btagup");     systVars.push_back("btagdown" );
	systVars.push_back("unbtagup");   systVars.push_back("unbtagdown" );
	if(isTTbarMC) {systVars.push_back("topptuncup"); systVars.push_back("topptuncdown"); }
	//systVars.push_back(); systVars.push_back();

	if(isTTbarMC) { systVars.push_back("pdfup"); systVars.push_back("pdfdown"); }
	cout << "Systematics will be computed for this analysis - this will take a bit" << endl;
	}

size_t nSystVars(systVars.size());
	

// TODO: what is this: allWeightsURL ... "weightsFile"??
std::vector < std::string > allWeightsURL = runProcess.getParameter < std::vector < std::string > >("weightsFile");
std::string weightsDir (allWeightsURL.size ()? allWeightsURL[0] : "");
// weightsDir is not used
//  //shape uncertainties for dibosons
//  std::vector<TGraph *> vvShapeUnc;
//  if(isMC_ZZ || isMC_WZ)
//    {
//      TString weightsFile=weightsDir+"/zzQ2unc.root";
//      TString dist("zzpt");
//      if(isMC_WZ) { weightsFile.ReplaceAll("zzQ2","wzQ2"); dist.ReplaceAll("zzpt","wzpt"); }
//      gSystem->ExpandPathName(weightsFile);
//      TFile *q2UncF=TFile::Open(weightsFile);
//      if(q2UncF){
//    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_up") ) );
//    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_down") ) );
//    q2UncF->Close();
//      }
//    }




	
//##############################################
//######## GET READY FOR THE EVENT LOOP ########
//##############################################
size_t totalEntries(0);

TFile* summaryFile = NULL;
TTree* summaryTree = NULL; //ev->;

//  
//  if(saveSummaryTree)
//    {
//      TDirectory* cwd = gDirectory;
//      std::string summaryFileName(outUrl); 
//      summaryFileName.replace(summaryFileName.find(".root", 0), 5, "_summary.root");
//      
//      summaryFile = new TFile(summaryFileName.c_str() "recreate");
//      
//      summaryTree = new TTree("Events", "Events");
//      KEY: TTreeMetaData;1
//      KEY: TTreeParameterSets;1
//      KEY: TTreeParentage;1
//      KEY: TTreeEvents;1
//      KEY: TTreeLuminosityBlocks;1
//      KEY: TTreeRuns;
//      summaryTree->SetDirectory(summaryFile);  // This line is probably not needed
//      
//      summmaryTree->Branch(
//
//      cwd->cd();
//    }
//


// Data-driven tau fakerate background FAKERATES

TFile * tau_fake_rate1_file = TFile::Open(runProcess.getParameter < std::string > ("dataDriven_tauFakeRates1") .c_str() );

// TODO: so the two fakerates are done 2ce -- two files and each file has qcd and wjets jistos
// rate1 = file1 = JetHT data file
TH3F * tau_fake_rate1_jets_histo_q = (TH3F *) tau_fake_rate1_file->Get("HLTjet_qcd_jets_distr");
TH3F * tau_fake_rate1_taus_histo_q = (TH3F *) tau_fake_rate1_file->Get("HLTjet_qcd_tau_jets_distr");
TH3F * tau_fake_rate1_jets_histo_w = (TH3F *) tau_fake_rate1_file->Get("HLTmu_qcd_jets_distr");
TH3F * tau_fake_rate1_taus_histo_w = (TH3F *) tau_fake_rate1_file->Get("HLTmu_qcd_tau_jets_distr");


TFile * tau_fake_rate2_file = TFile::Open(runProcess.getParameter < std::string > ("dataDriven_tauFakeRates2") .c_str() );
// rate2 = file2 = SingleMuon data file
TH3F * tau_fake_rate2_jets_histo_q = (TH3F *) tau_fake_rate2_file->Get("HLTmu_qcd_jets_distr");
TH3F * tau_fake_rate2_taus_histo_q = (TH3F *) tau_fake_rate2_file->Get("HLTmu_qcd_tau_jets_distr");
TH3F * tau_fake_rate2_jets_histo_w = (TH3F *) tau_fake_rate2_file->Get("HLTmu_wjets_jets_distr");
TH3F * tau_fake_rate2_taus_histo_w = (TH3F *) tau_fake_rate2_file->Get("HLTmu_wjets_tau_jets_distr");


Double_t tau_fake_rate_histo1_fraction = runProcess.getParameter < Double_t > ("tau_fake_rate_histo1_fraction");


// dilepton fake rates
TFile * tau_fake_rate_file_dileptons = TFile::Open(runProcess.getParameter < std::string > ("dataDriven_tauFakeRates_dileptons") .c_str() );
TH3F * tau_fake_rate_jets_histo_elmu = (TH3F *) tau_fake_rate_file_dileptons->Get("elmu_passjets_jets_distr");
TH3F * tau_fake_rate_taus_histo_elmu = (TH3F *) tau_fake_rate_file_dileptons->Get("elmu_passjets_tau_jets_distr");
TH3F * tau_fake_rate_jets_histo_mumu = (TH3F *) tau_fake_rate_file_dileptons->Get("mumu_passjets_jets_distr");
TH3F * tau_fake_rate_taus_histo_mumu = (TH3F *) tau_fake_rate_file_dileptons->Get("mumu_passjets_tau_jets_distr");

// Make fake rate projections for smooth procedure
//TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);

// QCD
TH1D* tau_fake_rate1_jets_histo_q_x = (TH1D*) tau_fake_rate1_jets_histo_q->Project3D("x");
TH1D* tau_fake_rate1_jets_histo_q_y = (TH1D*) tau_fake_rate1_jets_histo_q->Project3D("y");
TH1D* tau_fake_rate1_jets_histo_q_z = (TH1D*) tau_fake_rate1_jets_histo_q->Project3D("z");
TH1D* tau_fake_rate1_taus_histo_q_x = (TH1D*) tau_fake_rate1_taus_histo_q->Project3D("x");
TH1D* tau_fake_rate1_taus_histo_q_y = (TH1D*) tau_fake_rate1_taus_histo_q->Project3D("y");
TH1D* tau_fake_rate1_taus_histo_q_z = (TH1D*) tau_fake_rate1_taus_histo_q->Project3D("z");
FakeRateProjections frates_qcd_jets_proj;
frates_qcd_jets_proj.x = tau_fake_rate1_jets_histo_q_x;
frates_qcd_jets_proj.y = tau_fake_rate1_jets_histo_q_y;
frates_qcd_jets_proj.z = tau_fake_rate1_jets_histo_q_z;
FakeRateProjections frates_qcd_taus_proj;
frates_qcd_taus_proj.x = tau_fake_rate1_taus_histo_q_x;
frates_qcd_taus_proj.y = tau_fake_rate1_taus_histo_q_y;
frates_qcd_taus_proj.z = tau_fake_rate1_taus_histo_q_z;

// WJets
TH1D* tau_fake_rate2_jets_histo_w_x = (TH1D*) tau_fake_rate2_jets_histo_w->Project3D("x");
TH1D* tau_fake_rate2_jets_histo_w_y = (TH1D*) tau_fake_rate2_jets_histo_w->Project3D("y");
TH1D* tau_fake_rate2_jets_histo_w_z = (TH1D*) tau_fake_rate2_jets_histo_w->Project3D("z");
TH1D* tau_fake_rate2_taus_histo_w_x = (TH1D*) tau_fake_rate2_taus_histo_w->Project3D("x");
TH1D* tau_fake_rate2_taus_histo_w_y = (TH1D*) tau_fake_rate2_taus_histo_w->Project3D("y");
TH1D* tau_fake_rate2_taus_histo_w_z = (TH1D*) tau_fake_rate2_taus_histo_w->Project3D("z");
FakeRateProjections frates_wjets_jets_proj;
frates_wjets_jets_proj.x = tau_fake_rate2_jets_histo_w_x;
frates_wjets_jets_proj.y = tau_fake_rate2_jets_histo_w_y;
frates_wjets_jets_proj.z = tau_fake_rate2_jets_histo_w_z;
FakeRateProjections frates_wjets_taus_proj;
frates_wjets_taus_proj.x = tau_fake_rate2_taus_histo_w_x;
frates_wjets_taus_proj.y = tau_fake_rate2_taus_histo_w_y;
frates_wjets_taus_proj.z = tau_fake_rate2_taus_histo_w_z;

// ELMU
TH1D* tau_fake_rate_jets_histo_elmu_x = (TH1D*) tau_fake_rate_jets_histo_elmu->Project3D("x");
TH1D* tau_fake_rate_jets_histo_elmu_y = (TH1D*) tau_fake_rate_jets_histo_elmu->Project3D("y");
TH1D* tau_fake_rate_jets_histo_elmu_z = (TH1D*) tau_fake_rate_jets_histo_elmu->Project3D("z");
TH1D* tau_fake_rate_taus_histo_elmu_x = (TH1D*) tau_fake_rate_taus_histo_elmu->Project3D("x");
TH1D* tau_fake_rate_taus_histo_elmu_y = (TH1D*) tau_fake_rate_taus_histo_elmu->Project3D("y");
TH1D* tau_fake_rate_taus_histo_elmu_z = (TH1D*) tau_fake_rate_taus_histo_elmu->Project3D("z");
FakeRateProjections frates_elmu_jets_proj;
frates_elmu_jets_proj.x = tau_fake_rate_jets_histo_elmu_x;
frates_elmu_jets_proj.y = tau_fake_rate_jets_histo_elmu_y;
frates_elmu_jets_proj.z = tau_fake_rate_jets_histo_elmu_z;
FakeRateProjections frates_elmu_taus_proj;
frates_elmu_taus_proj.x = tau_fake_rate_taus_histo_elmu_x;
frates_elmu_taus_proj.y = tau_fake_rate_taus_histo_elmu_y;
frates_elmu_taus_proj.z = tau_fake_rate_taus_histo_elmu_z;

// MUMU
TH1D* tau_fake_rate_jets_histo_mumu_x = (TH1D*) tau_fake_rate_jets_histo_mumu->Project3D("x");
TH1D* tau_fake_rate_jets_histo_mumu_y = (TH1D*) tau_fake_rate_jets_histo_mumu->Project3D("y");
TH1D* tau_fake_rate_jets_histo_mumu_z = (TH1D*) tau_fake_rate_jets_histo_mumu->Project3D("z");
TH1D* tau_fake_rate_taus_histo_mumu_x = (TH1D*) tau_fake_rate_taus_histo_mumu->Project3D("x");
TH1D* tau_fake_rate_taus_histo_mumu_y = (TH1D*) tau_fake_rate_taus_histo_mumu->Project3D("y");
TH1D* tau_fake_rate_taus_histo_mumu_z = (TH1D*) tau_fake_rate_taus_histo_mumu->Project3D("z");
FakeRateProjections frates_mumu_jets_proj;
frates_mumu_jets_proj.x = tau_fake_rate_jets_histo_mumu_x;
frates_mumu_jets_proj.y = tau_fake_rate_jets_histo_mumu_y;
frates_mumu_jets_proj.z = tau_fake_rate_jets_histo_mumu_z;
FakeRateProjections frates_mumu_taus_proj;
frates_mumu_taus_proj.x = tau_fake_rate_taus_histo_mumu_x;
frates_mumu_taus_proj.y = tau_fake_rate_taus_histo_mumu_y;
frates_mumu_taus_proj.z = tau_fake_rate_taus_histo_mumu_z;

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

Variation jet_m_systematic_variation = Variation::NOMINAL; // FIXME: it should be in some headers, included before... but remake it somehow


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

// Electron energy scale, based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer and adapted to this framework
// v1
//string EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015";
//EpCombinationTool theEpCombinationTool;
//theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC 
//ElectronEnergyCalibratorRun2 ElectronEnCorrector(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
//ElectronEnCorrector.initPrivateRng(new TRandom(1234));


// --------------------------------------- lepton efficiencies
LeptonEfficiencySF lepEff;

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
TH2F* bTaggingEfficiencies_b_alljet    = (TH2F*) bTaggingEfficiencies_file->Get("btag_b_hadronFlavour_candidates");
TH2F* bTaggingEfficiencies_b_tagged    = (TH2F*) bTaggingEfficiencies_file->Get("btag_b_hadronFlavour_candidates_tagged");
TH2F* bTaggingEfficiencies_c_alljet    = (TH2F*) bTaggingEfficiencies_file->Get("btag_c_hadronFlavour_candidates");
TH2F* bTaggingEfficiencies_c_tagged    = (TH2F*) bTaggingEfficiencies_file->Get("btag_c_hadronFlavour_candidates_tagged");
TH2F* bTaggingEfficiencies_udsg_alljet = (TH2F*) bTaggingEfficiencies_file->Get("btag_udsg_hadronFlavour_candidates");
TH2F* bTaggingEfficiencies_udsg_tagged = (TH2F*) bTaggingEfficiencies_file->Get("btag_udsg_hadronFlavour_candidates_tagged");


//struct bTaggingEfficiencyHistograms {
//	TH2F* b_alljet   ;
//	TH2F* b_tagged   ;
//	TH2F* c_alljet   ;
//	TH2F* c_tagged   ;
//	TH2F* udsg_alljet;
//	TH2F* udsg_tagged;
//	};
struct bTaggingEfficiencyHistograms bEffs;

bEffs.b_alljet    = bTaggingEfficiencies_b_alljet   ;
bEffs.b_tagged    = bTaggingEfficiencies_b_tagged   ;
bEffs.c_alljet    = bTaggingEfficiencies_c_alljet   ;
bEffs.c_tagged    = bTaggingEfficiencies_c_tagged   ;
bEffs.udsg_alljet = bTaggingEfficiencies_udsg_alljet;
bEffs.udsg_tagged = bTaggingEfficiencies_udsg_tagged;

// Prescriptions taken from: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
// TODO: update to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X

// b-tagging working points for 50ns 
//   (pfC|c)ombinedInclusiveSecondaryVertexV2BJetTags
//      v2CSVv2L 0.605
//      v2CSVv2M 0.890
//      v2CSVv2T 0.970
double
	btagLoose(0.605), // not used anywhere in the code
	//btagMedium(0.890), // used twice in the code
	btagMedium(0.8), // new medium working point
	btagTight(0.970); // not used anywhere in the code

//b-tagging: scale factors
//beff and leff must be derived from the MC sample using the discriminator vs flavor
//the scale factors are taken as average numbers from the pT dependent curves see:
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
BTagSFUtil btsfutil;
double beff(0.68), sfb(0.99), sfbunc(0.015);
double leff(0.13), sfl(1.05), sflunc(0.12);

// Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
sfb = 0.861; // SF is not used --- BTagCalibrationReader btagCal instead
// sbbunc =;
beff = 0.559;
beff = 0.747;

// new btag calibration
// TODO: check callibration readers in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
// and latest standalone callibrator:
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h

// Setup calibration readers
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
// BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_76X.csv");
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_ichep_80X.csv");
BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_Moriond17_B_H.csv");

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




// ------------------------------ electron IDs
// does not appear anywhere in the code at all
// TString
	// electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
	// electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

// --------------------------------------- pileup weighting
// pile-up is done directly with direct_pileup_reweight
// v7-10 pile-up is back
std::vector<double> direct_pileup_reweight = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct");
std::vector<double> direct_pileup_reweight_up = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_up");
std::vector<double> direct_pileup_reweight_down = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_down");
	
gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

// --------------------------------------- hardcoded MET filter
// trying met-filters in HLT paths (76X MINIAODs should have it)
/*
patUtils::MetFilter metFiler;
if(!isMC)
	{
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_ecalscn1043093.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_ecalscn1043093.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_ecalscn1043093.txt");
	}
*/

// ---------------------------------------- HLT trigger efficiencies
// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Muon_reconstruction_identificati
// -- SingleMu Triggers
// ??? https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
// ??? -- Triggering electrons MVA-based WPs, scale factors for 76X
//
// in principle, one needs the electron/muon which triggered the HLT
// and according to its' pt, eta the scale factor is extracted
// TODO: procedure of getting the fired lepton
//
// only single-lepton HLTs are used, and their scale factors are taken accordingly
// though events with both HLTs on are also considered
// how to deal with them having only single-lepton SF?
//
// the formula can be such:
// e = 1 - (1 - e_mu)(1 - e_el)
// --- e_mu and e_el are trigger eff-s for single lepton
// the formula is taken from TOP-16-017
// https://indico.cern.ch/event/546389/contributions/2218845/attachments/1299941/1940277/khvastunov_28Jun_2016_TopPAG.pdf
//
// the study measures the efficiencies themselves
// I use the values provided by POGs
// muon POG only has efficiency SF for IsoMu20_OR_IsoTkMu20 yet
// (runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins)
// EGamma doesn't have any, thus there will be 1 on its' place

// v1
//TString muon_HLTeff_filename(string(std::getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/analysis/hlt-triggers/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root");
// .Data() returns char *
// c_str as well?
//TFile* muon_HLTeff_file = TFile::Open(muon_HLTeff_filename.Data());
//TH2F* muon_HLTeff_TH2F = (TH2F*) muon_HLTeff_file->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio");

// I do access the muon_HLTeff_TH2F histo with muon_HLTeff_TH2F->GetBin(eta, pt)

// To USE:
/*
muon_HLTeff_TH2F->FindBin
Double_t weight *= muon_HLTeff_TH2F->GetBinContent( muon_HLTeff_TH2F->FindBin(fabs(leta), electron.pt()) );

-- and it should return the scale factor for the HLT
*/


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

// Removed the SmartSelectionMonitor
// SmartSelectionMonitor mon;

TH1D* singlelep_ttbar_initialevents  = (TH1D*) new TH1D("singlelep_ttbar_init",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_preselectedevents = (TH1D*) new TH1D("singlelep_ttbar_presele",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
/* These counter histograms are dissabled -- use int for that
TH1D* singlelep_ttbar_selected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 

TH1D* singlelep_ttbar_maraselected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_maraselected_el_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
*/

// Kinematic parameters of the decay
TLorentzVector pl, plb, pb, pbb, prest;

// -------------------------------
// Here the output histograms and other object should be initialized




//##############################################
//########           EVENT LOOP         ########
//##############################################
//loop on all the events
printf ("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");

int nMultiChannel(0);

// FIXME: it's initialization of a rare control point, make it automatic somehow? initialize all?
// to check if processing goes well now
// increment( string("weight_passed_oursel"), 0. );

unsigned int iev = 0;
double weightflow_control_el_selection = 0;
double weightflow_control_mu_selection = 0;
double weightflow_control_elel_selection = 0;
double weightflow_control_mumu_selection = 0;
double weightflow_control_elmu_selection = 0;

for(size_t f=0; f<urls.size();++f)
	{
	cout << "Processing file: " << urls[f].c_str() << "\n";
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());

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
						mc_decay += find_W_decay(W);
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
						mc_decay += find_W_decay(W);
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
		mc_decay = string("_") + mc_decay; // so we'll have "_" or "_mcdecay"

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


		//singlelep_ttbar_initialevents->Fill(1);
		iev++;
		totalEntries++;
		if (iev % treeStep == 0)
			{
			printf (".");
			if(!debug) fflush (stdout); // Otherwise debug messages are flushed
			}

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
		double weight_bTaggingSF (1.0);

		// --------------------------------- tau ID SF
		double weight_tauIDsf         (1.0);
		double weight_without_tauIDsf (1.0);

		// rawWeight is everything but Pile-Up
		double rawWeight        (1.0);

		double HLT_efficiency_sf = 1.0;

		// final weight of the event
		double weight           (1.0);
		double weight_up        (1.0);
		double weight_down      (1.0);
		// and systematic corrections? TODO: check how TotalWeight_plus is used?




		// -------------------------------------------------- FIRST SECTION OF MC WEIGHTS, [1, 10]

		// increment( string("weightflow_n_miniaod_events"), 1.0 );
		// iniweight 1
		fill_1d(string("weightflow_mu"), 300, 0, 300,   1, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   1, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 1, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 1, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 1, weight);

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

		weight *= weight_TopPT; // how is the overall integral of MC?
		// the MC is lumi-xsec scaled to weightflow_weighted_miniaod_events

		// wighttoppt 2
		fill_1d(string("weightflow_mu"), 300, 0, 300,   2, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   2, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 2, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 2, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 2, weight);

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


		// -------------- it should be the creeppish merging of LO and NLO sets
		// not used now at all?
		// This must remain deactivated if you use HT-binned samples (it was for pthat-binned samples)
		// if (isV0JetsMC)
		//   {
		//   fwlite::Handle < LHEEventProduct > lheEPHandle;
		//   lheEPHandle.getByLabel (ev, "externalLHEProducer");
		//   mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
		//   if (lheEPHandle->hepeup ().NUP > 5)  continue;
		//   mon.fillHisto ("nupfilt", "", lheEPHandle->hepeup ().NUP, 1);
		//   }

		// HT-binned samples stitching: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples
		/*
		if(isV0JetsMC)
			{
			// access generator level HT               
			fwlite::Handle<LHEEventProduct> lheEventProduct;
			lheEventProduct.getByLabel(ev, "externalLHEProducer");
			//edm::Handle<LHEEventProduct> lheEventProduct;
			//ev.getByLabel( 'externalLHEProducer', lheEventProduct);
			const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup(); 
			std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
			double lheHt = 0.;
			size_t numParticles = lheParticles.size();
			for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle )
				{
				int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
				int status = lheEvent.ISTUP[idxParticle];
				if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) )
					{
					// quarks and gluons
					lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.));
					// first entry is px, second py
					}                                        
				}
			if(debug) cout << "Sample: " << dtag << ", lheHt: " << lheHt << ", scale factor from spreadsheet: " << patUtils::getHTScaleFactor(dtag, lheHt) << endl;
			// getHTScaleFactor works on combining several LO datasets with NLO
			// now one 1 NLO dataset is used for both WJets and DYJets
			// thus it is commented out here
			//weight_Gen *=   patUtils::getHTScaleFactor(dtag, lheHt);
			}
		*/

		fill_1d(string("weight_Gen"), 200, -2., 2., weight_Gen, 1);

		weight *= weight_Gen;
		weight_up *= weight_Gen;
		weight_down *= weight_Gen;
		rawWeight *=weight_Gen;

		// weightgen
		fill_1d(string("weightflow_el"), 300, 0, 300,   3, weight);
		fill_1d(string("weightflow_mu"), 300, 0, 300,   3, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 3, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 3, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 3, weight);

		fill_1d(string("eventflow_el"), 300, 0, 300,   3, 1);
		fill_1d(string("eventflow_mu"), 300, 0, 300,   3, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 3, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 3, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 3, 1);
		// ------------------------------- count N good verteces
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
		nGoodPV = vtx.size();
		goodPV = vtx[0];
		/*
		for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
			{
			if(utils::isGoodVertex(vtx[ivtx]))
				{
				if(nGoodPV==0) goodPV=vtx[ivtx];
				nGoodPV++;
				}
			}
		*/

		// ----------------------------------------- Apply pileup reweighting
		// why don't use nGoodPV for Pile-Up?
		unsigned int num_inters = 0, num_inters_raw = 0;
		double weight_pu_test = weight;
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

		weight *= weight_PU;
		weight_up *= weight_PU_up;
		weight_down *= weight_PU_down;

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights
		// sum_weights += weight;
		// sum_weights_raw += rawWeight;

		// int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight);

		// puweight
		fill_1d(string("weightflow_mu"), 300, 0, 300,   4, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   4, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 4, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 4, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 4, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   4, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   4, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 4, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 4, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 4, 1);

		// pu distrs
		fill_1d( string("pileup_beforetrig_num_inters_rawWeight"), 100, 0, 100, num_inters, rawWeight);
		fill_1d( string("pileup_beforetrig_num_inters_weight"),    100, 0, 100, num_inters, weight);
		// nGoodPV = vtx.size() now

		// vtx.size
		// now nGoodPV = vtx.size()
		//fill_1d( string("pileup_beforetrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
		//fill_1d( string("pileup_beforetrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);

		// pu distrs
		fill_1d( string("pileup_beforetrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_beforetrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
		// nGoodPV = vtx.size() now

		//fill_1d( string("pileup_passtrig_rawweight_pernuminters"), 100, 0, 100, nGoodPV, rawWeight);
		//fill_1d( string("pileup_passtrig_weight_pernuminters"),    100, 0, 100, nGoodPV, weight);

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		// increment( string("weightflow_weighted_raw_miniaod_events"), rawWeight );
		// increment( string("weightflow_weighted_miniaod_events"), weight );
		// increment( string("weightflow_weighted_up_miniaod_events"), weight_up );
		// increment( string("weightflow_weighted_down_miniaod_events"), weight_down );

		// // fill_pu( string("pileup_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
		// fill_pu( string("pileup_weight_perrawvtxsize"), vtx.size(), weight_pu_test);

		// fill_pu( string("pileup_ini_rawweight_pergoodpv"), nGoodPV, rawWeight);
		// fill_pu( string("pileup_ini_weight_pergoodpv"), nGoodPV, weight);
		// fill_pu( string("pileup_ini_weight_up_pergoodpv"), nGoodPV, weight_up);
		// fill_pu( string("pileup_ini_weight_down_pergoodpv"), nGoodPV, weight_down);

		// fill_pu( string("pileup_rawweight_pernuminters"), num_inters, rawWeight);
		// fill_pu( string("pileup_weight_pernuminters"), num_inters, weight);
		// fill_pu( string("pileup_weight_up_pernuminters"), num_inters, weight_up);
		// fill_pu( string("pileup_weight_down_pernuminters"), num_inters, weight_down);

		//num_inters = 1;
		// TODO: turn it back on for sakes of additional check
		// -- the weights are rarely negative now

		// if (num_inters>99) num_inters = 99;
		// if (nGoodPV>100) nGoodPV = 99;
		// event_pergoodpv_weight[nGoodPV] += weight;
		// //if (num_inters<0)  num_inters = 0;
		// if (weight_Gen<0)
		// 	{
		// 	increment( string("negative_events"), 1 );
		// 	fill_pu( string("pileup_negative_weight_pernuminters"), num_inters, weight);
		// 	fill_pu( string("pileup_negative_rawweight_pernuminters"), num_inters, rawWeight);

		// 	negative_event_nvtx[num_inters] += 1;
		// 	negative_event_pernvtx_weight[num_inters] += weight;
		// 	negative_event_pergoodpv_weight[nGoodPV] += weight;
		// 	}
		// else
		// 	{
		// 	increment( string("positive_events"), 1 );
		// 	fill_pu( string("pileup_positive_weight_pernuminters"), num_inters, weight);
		// 	fill_pu( string("pileup_positive_rawweight_pernuminters"), num_inters, rawWeight);
		// 	positive_event_nvtx[num_inters] += 1;
		// 	positive_event_pernvtx_weight[num_inters] += weight;
		// 	positive_event_pergoodpv_weight[nGoodPV] += weight;
		// 	}

		// -------------------------------------   Basic event selection

		// ---------------------- Orthogonalize Run2015B PromptReco+17Jul15 mix
		// let's remove Run2015B
		// if(isRun2015B)
		// {
		// if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
		// }

		// it's not needed with the latest versions of RunB rereconstruction

		// -------------------------------------------------- FIRST SECTION OF MC WEIGHT is over

		// FISRT SECTION SUM
		fill_1d(string("weightflow_mu"), 300, 0, 300,   10, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   10, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 10, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 10, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 10, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   10, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   10, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 10, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 10, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 10, 1);





		// -------------------------------------------------- SECOND SECTION -- cuts on trigger, lumi etc [11, 20]

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

		// passlumi 5
		fill_1d(string("weightflow_mu"), 300, 0, 300,   11, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   11, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 11, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 11, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 11, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   11, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   11, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 11, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 11, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 11, 1);

		// --------------------------------------------- HLT TRIGGER
		// ---------------- and require compatibilitiy of the event with the PD

		/* passed from job cfg.py
		string //jetHLT("HLT_PFJet40_v*"), // jetHLT("HLT_AK4PFJet30_v*"),
			muHLT_MC1("HLT_IsoMu24_v4"), muHLT_MC2("HLT_IsoTkMu24_v4"),
			muHLT_Data1("HLT_IsoMu24_v*"), muHLT_Data2("HLT_IsoTkMu24_v*"),
			elHLT_MC("HLT_Ele32_eta2p1_WPTight_Gsf_v8"), elHLT_Data("HLT_Ele32_eta2p1_WPTight_Gsf_v*");
		*/

			/* more triggers:
			for Run2016B,C, D, E, F and G 25 ns data with RunIISpring16reHLT80 MC (7th October Update)
			https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2016B_C_D_E_F_and_G_25_ns_dat
			          Data                              MC
			Muon      HLT_Iso(Tk)Mu24_v*                HLT_Iso(Tk)Mu24_v2
			Electron  HLT_Ele32_eta2p1_WPTight_Gsf_v*   HLT_Ele32_eta2p1_WPTight_Gsf_v3 
			*/

			/*
			// 2015, 76X MC
			// utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v*") :
			// 2016, 80X MC
			//true : // for noHLT MC
			// utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v2") : // recommended inn ttbar trig for reHLT
			//utils::passTriggerPatterns(tr, "HLT_Ele32_eta2p1_WPTight_Gsf_v3") : // recommended inn ttbar trig for reHLT Spring 16 MC
			utils::passTriggerPatterns(tr, "") : // recommended inn ttbar trig for Summer16 MC
			//utils::passTriggerPatterns(tr, "HLT_Ele27_eta2p1_WPTight_Gsf_v*") : // Using no-reHLT MC for now
			//utils::passTriggerPatterns(tr, "HLT_Ele*") : // Using no-reHLT MC for now
			// other trigger HLT_Ele27_eta2p1_WPTight_Gsf_v2
			// utils::passTriggerPatterns(tr, "HLT_Ele27_eta2p1_WPTight_Gsf_v*") :
			utils::passTriggerPatterns(tr, "*")

			// 2015, 76X MC
			// utils::passTriggerPatterns (tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*") // the efficiency scale factor are available for these only
			// utils::passTriggerPatterns (tr, "HLT_IsoMu18_v*", "HLT_IsoTkMu18_v*")
			// utils::passTriggerPatterns (tr, "HLT_IsoMu18_v*")
			// 2016, 80X MC
			//true : // for noHLT MC
			//utils::passTriggerPatterns (tr, "HLT_IsoMu24_v2", "HLT_IsoTkMu24_v2") : // recommended in ttbar trig for reHLT Spring16 MC
			utils::passTriggerPatterns (tr, "HLT_IsoMu24_v4", "HLT_IsoTkMu24_v4") : // recommended in ttbar trig Summer16 MC
			//utils::passTriggerPatterns (tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*") :
			utils::passTriggerPatterns (tr, "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*")
			*/

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
		bool eTrigger = false;
		bool muTrigger = false;

		if (isNoHLT)
			{
			eTrigger = true;
			muTrigger = true;
			}
		else if (tr.isValid())
			{
			eTrigger =    ( isMC ?
			utils::passTriggerPatterns(tr, elHLT_MC) :
			utils::passTriggerPatterns(tr, elHLT_Data)
			);
			muTrigger =   ( isMC ?
			utils::passTriggerPatterns (tr, muHLT_MC1, muHLT_MC2) :
			utils::passTriggerPatterns (tr, muHLT_Data1, muHLT_Data2)
			);
			}
		else return 233;

		//if(filterOnlySINGLEMU) {                    eTrigger = false; }
		//if(filterOnlySINGLEE)  { muTrigger = false;                   }
		// SingleMuon-dataset jobs process double-HLT events
		// SingleElectron-dataset jobs skip them

		// if data and SingleElectron dataset and both triggers -- skip event
		if (!debug) {
			if (!isMC && isSingleElectronDataset && eTrigger && muTrigger) continue;
		
			if (!(eTrigger || muTrigger)) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
		}

		// TODO: ----------------------------- HLT efficiency scale factors
		// one should run it on the fired trigger objects,
		// I run it on selection candidates now
		// which is done below, when the candidates are selected

		// double HLT_efficiency_sf = 1.0;

		//HLT_efficiency_sf *= eTrigger  ? eHLT_sf[] : 1 ;
		//HLT_efficiency_sf *= muTrigger ? muHLT_SF[] : 1 ;

		// passtrig 6
		fill_1d(string("weightflow_mu"), 300, 0, 300,   12, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   12, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 12, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 12, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 12, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   12, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   12, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 12, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 12, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 12, 1);
		// increment( string("eventflow_1_passed_trig"), weight ); // should not matter
		// increment( string("weightflow_weight_up_passed_trig"), weight_up ); // should not matter
		// increment( string("weightflow_weight_down_passed_trig"), weight_down ); // should not matter

		// sum_weights_passtrig_raw += rawWeight;
		// sum_weights_passtrig += weight;


		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )


		// if (eTrigger)
		// 	{
		// 	increment( string("weightflow_weight_passed_electron_trigger"), weight );
		// 	}

		// if (muTrigger)
		// 	{
		// 	increment( string("weightflow_weight_passed_muon_trigger"), weight );
		// 	}

		// if (eTrigger && muTrigger)
		// 	{
		// 	increment( string("weightflow_weight_passed_electron_and_muon_triggers"), weight );
		// 	}

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

		// pu distrs
		fill_1d( string("pileup_passtrig_num_inters_rawWeight"), 100, 0, 100, num_inters, rawWeight);
		fill_1d( string("pileup_passtrig_num_inters_weight"),    100, 0, 100, num_inters, weight);
		// nGoodPV = vtx.size() now

		// pu distrs
		fill_1d( string("pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
		// nGoodPV = vtx.size() now


		// fill_pu( string("pileup_passtrig_rawweight_pernuminters"), num_inters, rawWeight);
		// fill_pu( string("pileup_passtrig_weight_pernuminters"), num_inters, weight);

		// fill_pu( string("pileup_passtrig_rawweight_pergoodpv"), nGoodPV, rawWeight);
		// fill_pu( string("pileup_passtrig_weight_pergoodpv"), nGoodPV, weight);
		// fill_pu( string("pileup_passtrig_weight_up_pergoodpv"), nGoodPV, weight_up);
		// fill_pu( string("pileup_passtrig_weight_down_pergoodpv"), nGoodPV, weight_down);

		// fill_pu( string("pileup_passtrig_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
		// fill_pu( string("pileup_passtrig_weight_perrawvtxsize"), vtx.size(), weight_pu_test);


		// ------------------------------------------------- Apply MET filters
		//if( !isMC && !metFiler.passMetFilter( ev, isPromptReco)) continue;
		// New passMetFilter procedure from patUtils:
		//if( !isMC && !metFiler.passMetFilter( ev )) continue;
		// trying the HLTs path for metfilters:
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_76X_v2_produced_with_the
		// recommendations (6-2016):
		// Flag_HBHENoiseFilter TO BE USED
		// Flag_HBHENoiseIsoFilter TO BE USED
		// Flag_CSCTightHalo2015Filter TO BE USED
		// Flag_EcalDeadCellTriggerPrimitiveFilter TO BE USED
		// Flag_goodVertices TO BE USED
		// Flag_eeBadScFilter TO BE USED 
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
				cout << "----------- End of trigger list ----------" << endl;
				//return 0;
			}
			//std::vector<std::string>& patterns("Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_CSCTightHalo2015Filter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter");
			if (!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*", "Flag_EcalDeadCellTriggerPrimitiveFilter*"))
				continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_goodVertices")) continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter")) continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_globalTightHalo2016Filter")) continue;
		}
		

		if(debug)
			{
			cout << "met filters applied here" << endl;
			}

		// passmetfilters 7
		fill_1d(string("weightflow_mu"), 300, 0, 300,   13, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   13, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 13, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 13, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 13, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   13, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   13, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 13, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 13, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 13, 1);

		// -------------------------------------------------- SECOND SECTION OF event cuts is over
		fill_1d(string("weightflow_mu"), 300, 0, 300,   20, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   20, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 20, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 20, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 20, weight);

		fill_1d(string("eventflow_mu"), 300, 0, 300,   20, 1);
		fill_1d(string("eventflow_el"), 300, 0, 300,   20, 1);
		fill_1d(string("eventflow_elel"), 300, 0, 300, 20, 1);
		fill_1d(string("eventflow_elmu"), 300, 0, 300, 20, 1);
		fill_1d(string("eventflow_mumu"), 300, 0, 300, 20, 1);

		// ------------------------- event physics and the corresponding selection

		//------------------------- PROCESS OBJECTS
		// -------------------------------------------------- possible THIRD SECTION of MC WEIGHTS and corrections (with efficiency SFs used as event probability correction)

		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;




		// TODO: figure out quark status/tau status in MC
		// FIXME: Save time and don't load the rest of the objects when selecting by mctruthmode :)
		//bool hasTop(false);
		//int
		//ngenLeptonsStatus3(0),
		//ngenLeptonsNonTauSonsStatus3(0),
		//ngenTausStatus3(0),
		//ngenQuarksStatus3(0);

		//double tPt(0.), tbarPt(0.); // top pt reweighting - dummy value results in weight equal to 1 if not set in loop
		//float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
		// TODO: what is this??
		// there was some wague answer from Pietro.....
/*
		if(isMC)
			{
			// FIXME: Considering add support for different generators (based on PYTHIA6) for comparison.
			for(size_t igen=0; igen<gen.size(); igen++)
				{
				// FIXME: Should pass to the new status scheme from: https://github.com/cms-sw/cmssw/pull/7791
				// ////// if(!gen[igen].isHardProcess() && !gen[igen].isPromptFinalState()) continue;

				if(gen[igen].status() != 1 &&  gen[igen].status() !=2 && gen[igen].status() !=62 ) continue;
				int absid=abs(gen[igen].pdgId());
				// OK, so taus should be checked as status 2, and quarks as 71 or 23. More testing needed
				//if( absid==15 && hasWasMother(gen[igen]) ) cout << "Event " << iev << ", Particle " << igen << " has " << gen[igen].numberOfDaughters() << " daughters, pdgId " << gen[igen].pdgId() << " and status " << gen[igen].status() << ", mothers " << gen[igen].numberOfMothers() << ", pt " << gen[igen].pt() << ", eta " << gen[igen].eta() << ", phi " << gen[igen].phi() << ". isHardProcess is " << gen[igen].isHardProcess() << ", and isPromptFinalState is " << gen[igen].isPromptFinalState() << endl;


				////// if(absid==6 && gen[igen].isHardProcess()){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass)
				if(absid==6 && gen[igen].status()==62)
					{
					// particles of the hardest subprocess 22 : intermediate (intended to have preserved mass). Josh says 62 (last in chain)
					hasTop=true;
					// FIXME: Top pT reweighting. 13 TeV values not propagated yet, so not using.
					//if(isTTbarMC){
					//  if(gen[igen].get("id") > 0) tPt=gen[igen].pt();
					//  else                        tbarPt=gen[igen].pt();
					//}
					} 


				//if(!gen[igen].isPromptFinalState() ) continue;
				if( (gen[igen].status() != 1 && gen[igen].status()!= 2 ) || !hasWasMother(gen[igen])) continue;

				if((absid==11 || absid==13) && hasLeptonAsDaughter(gen[igen]))
					cout << "Electron or muon " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;

				if((absid==11 || absid==13) && gen[igen].status()==1)
					{
					ngenLeptonsStatus3++;

					if(!hasTauAsMother(gen[igen]))
						ngenLeptonsNonTauSonsStatus3++;
					}

				if(absid==15 && gen[igen].status()==2 )
					{
					ngenTausStatus3++; // This should be summed to ngenLeptonsStatus3 for the dilepton final states, not summed for the single lepton final states.
					//if(hasLeptonAsDaughter(gen[igen]))
					//	cout << "Tau " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
					}

				if(debug && (ngenTausStatus3==1 && ngenLeptonsStatus3==1 )  ) cout << "Event: " << iev << ". Leptons: " << ngenLeptonsStatus3 << ". Leptons notaus: " << ngenLeptonsNonTauSonsStatus3 << ". Taus: " << ngenTausStatus3 << ". Quarks: " << ngenQuarksStatus3 << endl;
						
				// Dileptons:
				//    ttbar dileptons --> 1
				//    ttbar other     --> 2
				if(mctruthmode==1 && (ngenLeptonsStatus3+ngenTausStatus3!=2 || !hasTop )) continue;
				if(mctruthmode==2 && (ngenLeptonsStatus3+ngenTausStatus3==2 || !hasTop )) continue;
				// FIXME: port tt+bb splitting from 8 TeV (check the reference to the matched genjet)
				//if(mcTruthMode==1 && (ngenLeptonsStatus3!=2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
				//if(mcTruthMode==2 && (ngenLeptonsStatus3==2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
				//if(mcTruthMode==3 && (ngenBQuarksStatus23<4 || !hasTop))                           continue;
						
				// lepton-tau:
				//    ttbar ltau      --> 3
				//    ttbar dileptons --> 4
				//    ttbar ljets     --> 5
				//    ttbar hadrons   --> 6
				if(mctruthmode==3 && (ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1  || !hasTop )) continue; // This is bugged, as it is obvious
				if(mctruthmode==4 && (ngenLeptonsNonTauSonsStatus3!=2                        || !hasTop )) continue;
				if(mctruthmode==5 && (ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1        || !hasTop )) continue;
						
				bool isHad(false);
				if (
					(ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1 ) &&
					(ngenLeptonsNonTauSonsStatus3!=2                      ) &&
					(ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1      ) 
					)
				isHad=true;
					
				//if(mctruthmode==6 && (ngenLeptonsNonTauSonsStatus3!=0 || ngenTausStatus3!=0  || !hasTop )) continue;
				if(mctruthmode==6 && (!isHad || !hasTop )) continue;
				}
			}

		if(debug) cout << "DEBUG: Event was not stopped by the ttbar sample categorization (either success, or it was not ttbar)" << endl;
*/


		// ------------------------------------ actual particles?

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
		metsHandle.getByLabel(ev, "slimmedMETs");
		if(metsHandle.isValid() ) mets = *metsHandle;
		pat::MET met = mets[0];
		// LorentzVector met = mets[0].p4 ();

		fill_1d(string("control_met_slimmedMETs_pt"), 200, 0., 200., met.pt(), weight);


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

		processElectrons_ID_ISO_Kinematics(electrons, goodPV, rho, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
			35., 2.4, 15., 2.5, selElectrons, elDiff, nVetoE, false, debug);

		if(debug){
			cout << "processed electrons" << endl;
			}


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
		processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,
			30., 2.4, 10., 2.5, selMuons, muDiff, nVetoMu, false, debug);

		if(debug){
			cout << "processed muons" << endl;
			}

		// Finally, merge leptons for cross-cleaning with taus and jets, and other conveniences:

		std::vector<patUtils::GenericLepton> selLeptons;
		for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
		for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
		std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);


		// ------------------------------------- Propagate lepton energy scale to MET
		// Propagate now (v13)
		/* no lepton corrections propagation v13.1
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

		pat::TauCollection IDtaus, selTaus, selTaus_JetTauFakeRate;

		processTaus_ID_ISO(taus, weight, tau_decayMode, tau_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, debug);

		if(debug){
			cout << "selected taus [individual]" << endl;
			}

		//int processTaus_Kinematics(pat::TauCollection& taus,          // input
		//	double weight,
		//	double pt_cut, double eta_cut,
		//	pat::TauCollection& selTaus,                          // output
		//	bool record, bool debug) // more output

		processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,       false, debug);
		processTaus_Kinematics(IDtaus, weight, jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selTaus_JetTauFakeRate, false, debug);

		// ------------------------------------------ select the taus cleaned from leptons


		//int crossClean_in_dR(pat::TauCollection& selTaus, std::vector<patUtils::GenericLepton>& leptons,
		//	float min_dR,
		//	pat::TauCollection& selTausNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		pat::TauCollection selTausNoLep, selTaus_JetTauFakeRate_NoLep;
		crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weight, string("selTausNoLep"),        false, debug);
		crossClean_in_dR(selTaus_JetTauFakeRate, selLeptons, 0.4, selTaus_JetTauFakeRate_NoLep, weight, string("selTaus_JetTauFakeRate_NoLep"), false, debug);

		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Measurement_in_Z_tautau_events
		// Medium MVA (no dR03) is 0.97 +- 0.05

		if (isMC && selTausNoLep.size() > 0) // 2016 data/MC tau ID efficiency (all discriminators, all pt and eta ranges) = 0.83 +- 0.06
			{
			double pre_weight_tauIDsf = 1;
			for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
				//pre_weight_tauIDsf *= (1 - 0.83 + r3->Gaus(0, 0.06)); // gaussian +- 0.06
				// recommendations update (noticed 8-11-2016, last page update 2016-11-07)
				// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_ID_efficiency
				// in ReReco it's different
				// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Measurement_in_Z_tautau_events
				pre_weight_tauIDsf *= (1 - 0.97 + r3->Gaus(0, 0.05)); // gaussian +- 0.05
			weight_tauIDsf = 1 - pre_weight_tauIDsf;
			}
			// TODO: should here be a normalization to all MC events?
		fill_1d(string("weight_tauIDsf"), 200, 0., 2.,   weight_tauIDsf, 1);

		weight_without_tauIDsf = weight;
		// weight *= weight_tauIDsf;// apply tau weight at selection

		if(debug){
			cout << "processed taus" << " N selTausNoLep = " << selTausNoLep.size() << endl;
			cout << "ID SF = " << weight_tauIDsf << endl;
			}

		//
		// JET/MET SELECTION
		//

		if(debug) cout << "Now update Jet Energy Corrections" << endl;
		//add scale/resolution uncertainties and propagate to the MET
		//utils::cmssw::updateJEC(jets, jesCor, totalJESUnc, rho, nGoodPV, isMC);

		// up to here jets were not processed in any way
		// now goes the procedure of corrections to jets and then METs




		// ----------------------- JETS CORRECTIONS, JEC, JER
		// ----------------------- UPDATE JEC

		//int processJets_CorrectJES_SmearJERnJES_ID_ISO(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
		//	bool isMC, double weight,
		//	double rho, unsigned int nGoodPV,
		//	FactorizedJetCorrector *jesCor,
		//	JetCorrectionUncertainty *totalJESUnc,
		//	double dR_max, // for jet matching in jet corrections smearing for MC
		//	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
		//	string& jetID,
		//	//double pt_cut, double eta_cut,
		//	TRandom3 *r3,   // the randomizer for the smearing
		//	LorentzVector& full_jet_corr, pat::JetCollection& selJets,                          // output
		//	bool record, bool debug) // more output

		LorentzVector full_jet_corr(0., 0., 0., 0.);
		pat::JetCollection IDjets;
		string jetID("Loose");

		processJets_CorrectJES_SmearJERnJES_ID_ISO(jets, genJets, isMC, weight, rho, nGoodPV, jesCor, totalJESUnc, 0.4/2,
			jet_resolution_in_pt, jet_resolution_sf_per_eta, jet_m_systematic_variation, jetID, r3, full_jet_corr, IDjets, true, debug);


		fill_3d(string("control_jet_full_jet_corr_pX_pY_pZ"), 10, -100., 100., 10, -100., 100., 10, -100., 100.,  full_jet_corr.X(), full_jet_corr.Y(), full_jet_corr.Z(), weight);
		// 1000 bins

		fill_2d(string("control_jet_full_jet_corr_pX_pY"), 100, -50., 50., 100, -50., 50.,  full_jet_corr.X(), full_jet_corr.Y(), weight);
		fill_1d(string("control_jet_full_jet_corr_pZ"),    100, -50., 50., full_jet_corr.Z(), weight);
		// 10 000 and 100 bins

		met.setP4(met.p4() - full_jet_corr); // just return the full correction and propagate in place

		fill_1d(string("control_met_slimmedMETs_fulljetcorrs_pt"), 200, 0., 200., met.pt(), weight);
		//fill_2d(string("control_met_slimmedMETs_fulljetcorrs_pt"), 200, 0., 200., 200, -4., 4., met.pt(), met.eta(), weight);

		//int processJets_Kinematics(pat::JetCollection& jets, // input
		//	//bool isMC,
		//	double weight,
		//	double pt_cut, double eta_cut,
		//	pat::JetCollection& selJets,                 // output
		//	bool record, bool debug) // more output

		pat::JetCollection selJets;
		processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jet_kino_cuts_pt, jet_kino_cuts_eta, selJets, true, debug);

		pat::JetCollection selJets_JetTauFakeRate; // for fake rates in dileptons
		processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jettaufr_jet_kino_cuts_pt, jettaufr_jet_kino_cuts_eta, selJets_JetTauFakeRate, true, debug);

		// ---------------------------- Clean jet collections from selected leptons
		// TODO: add gamma-cleaning as well?

		//int crossClean_in_dR(pat::JetCollection& selJets, std::vector<patUtils::GenericLepton>& selLeptons,
		//	float min_dR,
		//	pat::JetCollection& selJetsNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		pat::JetCollection selJetsNoLep;
		crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, debug);

		pat::JetCollection selJets_JetTauFakeRate_NoLep;
		crossClean_in_dR(selJets_JetTauFakeRate, selLeptons, 0.4, selJets_JetTauFakeRate_NoLep, weight, string("selJets_JetTauFakeRate_NoLep"), false, debug);

		/*
		// so, just for the fake rates:

		pat::JetCollection selJets30GeVNoLep;
		for (size_t ijet = 0; ijet < selJets30GeV.size(); ++ijet)
			{
			pat::Jet& jet = selJets30GeV[ijet];

			double minDRlj (9999.);

			for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
				minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));

			if (minDRlj < 0.4) continue;

			selJets30GeVNoLep.push_back(jet);
			}

		pat::JetCollection selJets20GeVNoLep;
		for (size_t ijet = 0; ijet < selJets20GeV.size(); ++ijet)
			{
			pat::Jet& jet = selJets20GeV[ijet];

			double minDRlj (9999.);

			for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
				minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));

			if (minDRlj < 0.4) continue;

			selJets20GeVNoLep.push_back(jet);
			}
		*/





		// ---------------------------- Clean jet collection from selected taus

		//int crossClean_in_dR(pat::JetCollection& selJets, pat::TauCollection& selTaus,
		//	float min_dR,
		//	pat::JetCollection& selJetsOut, // output
		//	string control_name,
		//	bool record, bool debug) // more output
		pat::JetCollection selJetsNoLepNoTau;
		crossClean_in_dR(selJetsNoLep, selTausNoLep, 0.4, selJetsNoLepNoTau, weight, string("selJetsNoLepNoTau"), false, debug);

		if(debug){
			cout << "processed jets" << endl;
			}


		// --------------------------- B-TAGGED JETS
		pat::JetCollection selBJets;

		//int processBJets_BTag(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
		//	BTagCalibrationReader& btagCal, BTagSFUtil& btsfutil,
		//	struct bTaggingEfficiencyHistograms& bEffs,
		//	string& b_tagger_label, float b_tag_WP,
		//	pat::JetCollection& selBJets,                          // output
		//	bool record, bool debug) // more output

		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
		string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		float btag_WP = 0.8484; // medium
		processBJets_BTag(selJetsNoLepNoTau, isMC, weight, weight_bTaggingSF, btagCal, btsfutil, bEffs, btagger_label, btag_WP, selBJets, true, debug);

		weight *= weight_bTaggingSF;
		fill_1d(string("weight_bTaggingSF"), 200, 0., 2.,   weight_bTaggingSF, 1.);

		if(debug){
			cout << "processed b-tagged jets" << endl;
			}


		// also for referense:
		pat::TauCollection selTausNoLepNoJet;
		crossClean_in_dR(selTausNoLep, selJetsNoLep, 0.4, selTausNoLepNoJet, weight, string("selTausNoLepNoJet"), false, debug);


		// -------------------------------------------------- all particles are selected

		// -------------------------------------------------- cCONTROL VALUES

		// fill_1d(string("weightflow_mu_passmetfilters"), 300, 0, 300,   7, weight);

		fill_1d(string("n_slimmedjets"),  10, 0, 10,   jets.size(), weight);
		fill_1d(string("n_selJets"),      10, 0, 10,   selJets.size(), weight);
		fill_1d(string("n_selJetsNoLep"), 10, 0, 10,   selJetsNoLep.size(), weight);
		fill_1d(string("n_selJetsNoLepNoTau"), 10, 0, 10,   selJetsNoLepNoTau.size(), weight);

		fill_1d(string("n_selBJets"), 10, 0, 10,   selBJets.size(), weight);

		fill_1d(string("n_slimmedtaus"),  10, 0, 10,   taus.size(), weight);
		fill_1d(string("n_selTaus"),      10, 0, 10,   selTaus.size(), weight);
		fill_1d(string("n_selTausNoLep"), 10, 0, 10,   selTausNoLep.size(), weight);
		fill_1d(string("n_selTausNoLepNoJet"), 10, 0, 10,   selTausNoLepNoJet.size(), weight);

		fill_1d(string("n_slimmedtaus"),  10, 0, 10,   taus.size(), weight);
		fill_1d(string("n_selTaus"),      10, 0, 10,   selTaus.size(), weight);
		fill_1d(string("n_selTausNoLep"), 10, 0, 10,   selTausNoLep.size(), weight);
		fill_1d(string("n_selTausNoLepNoJet"), 10, 0, 10,   selTausNoLepNoJet.size(), weight);

		fill_1d(string("n_slimmedmuons"),  10, 0, 10,   muons.size(), weight);
		fill_1d(string("n_selMuons"),      10, 0, 10,   selMuons.size(), weight);

		fill_1d(string("n_slimmedelectrons"),  10, 0, 10,   electrons.size(), weight);
		fill_1d(string("n_selElectrons"),      10, 0, 10,   selElectrons.size(), weight);





		// last record was:
		// fill_1i(string("weightflow_mu_passmetfilters"), 300, 0, 300,   7, weight);
		// fill_1i(string("weightflow_el_passmetfilters"), 300, 0, 300,   7, weight);
		// fill_1i(string("weightflow_elel_passmetfilters"), 300, 0, 300, 7, weight);
		// fill_1i(string("weightflow_elmu_passmetfilters"), 300, 0, 300, 7, weight);
		// fill_1i(string("weightflow_mumu_passmetfilters"), 300, 0, 300, 7, weight);

		// -------------------------------------------------- THIRD SECTION OF mc weights for event probability SF-s
		fill_1d(string("weightflow_mu"), 300, 0, 300,   30, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   30, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 30, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 30, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 30, weight);

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
		// isSingleMu = (abs(slepId)==13) && muTrigger && iso_lep;
		// isSingleE  = (abs(slepId)==11) && eTrigger  && iso_lep;
		bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
		// maybe it is worth trying not considering the nGoodPV?

		// TODO: test if trigger needed here at all
		//isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && muTrigger && clean_lep_conditions;
		//isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && eTrigger  && clean_lep_conditions;
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
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[1].pt(), weight);

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
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[1].pt(), weight);

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
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[1].pt(), weight);

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
		//unsigned int n_jets = selJets.size();
		//unsigned int n_bjets = selBJets.size();
		// unsigned int n_jets = selJetsNoLepNoTau.size();
		unsigned int n_jets = selJetsNoLep.size(); // noLep jet as in jet fake-rate study
		// unsigned int n_bjets = selSingleLepBJets.size();
		unsigned int n_bjets = selBJets.size();


		// ------------------------------------------ SINGLE LEPTON CHANNELS
		if (isSingleMu || isSingleE)
			{
			// in-channel selection/multiselect for leptons

			//singlelep_ttbar_preselectedevents->Fill(1);


			// lepton+tau mass > 12 GeV, as in dilepton case
			LorentzVector dileptonSystem (0, 0, 0, 0);
			if(n_taus>0)
				{
				dileptonSystem = selLeptons[0].p4() + selTausNoLep[0].p4();
				/*
				if (isMC && withTauIDSFs)
					{
					// weight_tauIDsf = 1 - (1 -  0.83 + r3->Gaus(0, 0.06))^n_taus;
					weight_tauIDsf = 1 - pow(1 -  0.83 + r3->Gaus(0, 0.06), n_taus);
					//weight_tauIDsf = 0.83;
					// weight *= weight_tauIDsf;
					}
				*/
				}
			weight *= weight_tauIDsf;
			//fill_1d(string("weight_tauIDsf_2"), 200, 0., 2.,   weight_tauIDsf, 1);

			// bool passJetRawSelection(selSingleLepJets.size()>1); // 2 jets
			//bool passJetSelection(selSingleLepJets.size()>1); // 2 jets // 2^4
			//bool passJetSelection(selJets.size()>1); // 2 jets // 2^4
			//bool passJetSelection(n_jets>1); // 2 jets // 2^4
			bool passJetSelection(n_jets>2); // >= 3 jets
			bool passMetSelection(met.pt()>40.); // MET > 40 // 2^3
			// bool passMetSelection(n_met.pt()>40.); // MET > 40 // 2^3
			//bool passBtagsSelection(selSingleLepBJets.size()>0); // 1 b jet // 2^2
			//bool passBtagsSelection(selBJets.size()>0); // 1 b jet // 2^2
			bool passBtagsSelection(n_bjets>0); // 1 b jet // 2^2
			// bool passTauSelection(n_taus==1); // only 1 tau // 2^1
			//bool passTauSelection(n_taus>0); // >= 1 tau in v8.8

			//bool passTauSelection(n_taus>0 && dileptonSystem.mass()>12.); // >= 1 tau in v8.8
			bool passTauSelection(n_taus>0); // >= 1 tau in v8.8
			bool passOS( n_taus>0 && n_leptons>0 ? selLeptons[0].pdgId() * selTausNoLep[0].pdgId() < 0 : 0); // Oposite sign // 2^0
			// bool passOS( n_taus>0 && n_leptons>0 ? selLeptons[0].pdgId() * selTausNoLepNoJet[0].pdgId() < 0 : 0); // Oposite sign // 2^0

			// check no w-mass di-jet system for fake rates
			// thus, the tau should not couple with another jet into W-mass system
			bool noWmassDiJet = true;
			//if (passTauSelection && passJetSelection && passBtagsSelection) // to compute less
			if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS)
				{
				LorentzVector dijetSystem (0, 0, 0, 0);
				for (int i=0; i<selJetsNoLep.size(); i++)
					{
					// record jets around tau in dijet mass
					// scope all jets
					// record their dijet_mass VS inverse dR to tau
					// inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
					//   return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2));
					// }
					pat::Tau& tau = selTausNoLep[0];
					pat::Jet& jet = selJetsNoLep[i];
					dijetSystem = jet.p4() + tau.p4();
					double dijet_mass = dijetSystem.mass();
					double dijet_momentum = dijetSystem.P(); // square of spacial part, use M() for magnitude of spacial part
					double inverse_dR = reco::deltaR(tau.eta(), tau.phi(), -jet.eta(), -jet.phi());
					//if (dijet_mass > 60 && dijet_mass < 100) // or make the window narrower?
					//fill_2d(string("slep_vanila_selection_dijet_mass_VS_dR"), 100, 0, 200, 50, 0, 4,  dijet_mass, inverse_dR, weight);
					fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
					fill_1d(string("slep_vanila_selection_njets"), 10, 0, 10,   n_jets, weight);
					// and bin in N jets
					if (n_jets == 3)
						fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum_nj3"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
						//fill_1d(string("slep_vanila_selection_dijet_mass_nj3"), 100, 0, 100,   dijet_mass, weight);
					else if (n_jets == 4)
						fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum_nj4"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
						//fill_1d(string("slep_vanila_selection_dijet_mass_nj4"), 100, 0, 100,   dijet_mass, weight);
					else if (n_jets == 5)
						fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum_nj5"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
						//fill_1d(string("slep_vanila_selection_dijet_mass_nj5"), 100, 0, 100,   dijet_mass, weight);
					else if (n_jets == 6)
						fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum_nj6"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
						//fill_1d(string("slep_vanila_selection_dijet_mass_nj6"), 100, 0, 100,   dijet_mass, weight);
					else //if (n_jets >  6)
						fill_2d(string("slep_vanila_selection_dijet_mass_VS_momentum_nj6p"), 100, 0, 200, 100, 0, 300,  dijet_mass, dijet_momentum, weight);
						//fill_1d(string("slep_vanila_selection_dijet_mass_nj6p"), 100, 0, 100,   dijet_mass, weight);
					//noWmassDiJet = false;
					//break;
					}
				}
			// add it into tau selection
			//passTauSelection &= noWmassDiJet;

			// MULTISELECT
			unsigned int multisel = 0;
			// multisel += (isSingleMu ? 1 : 0); //! should be 1
			// multisel += (isSingleE ? 2 : 0);
			multisel += (passJetSelection ? 1 : 0);
			multisel += (passMetSelection ? 2 : 0);
			multisel += (passBtagsSelection ? 4 : 0);
			multisel += (passTauSelection ? 8 : 0);
			multisel += (passOS ? 16 : 0);


			{
				fill_2d(string("control_met_singlelep_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
				for (int i=0; i<selJetsNoLep.size(); i++)
					fill_2d(string("control_jet_singlelep_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
				for (int i=0; i<selTausNoLep.size(); i++)
					fill_2d(string("control_tau_singlelep_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
				for (int i=0; i<selLeptons.size(); i++)
					fill_2d(string("control_lep_singlelep_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
			}

			if (passJetSelection) {
				fill_2d(string("control_met_singlelep_passjet_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
				for (int i=0; i<selJetsNoLep.size(); i++)
					fill_2d(string("control_jet_singlelep_passjet_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
				for (int i=0; i<selTausNoLep.size(); i++)
					fill_2d(string("control_tau_singlelep_passjet_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
				for (int i=0; i<selLeptons.size(); i++)
					fill_2d(string("control_lep_singlelep_passjet_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
				}

			if (isSingleMu)
				{
				fill_2d(string("control_met_singlemu_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
				for (int i=0; i<selJetsNoLep.size(); i++)
					fill_2d(string("control_jet_singlemu_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
				for (int i=0; i<selTausNoLep.size(); i++)
					fill_2d(string("control_tau_singlemu_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
				for (int i=0; i<selLeptons.size(); i++)
					fill_2d(string("control_lep_singlemu_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);

				if (passJetSelection) {
					record_jets_fakerate_distrs(string("singlemu_"), string("control_passjet"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					fill_2d(string("control_met_singlemu_passjet_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singlemu_passjet_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singlemu_passjet_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singlemu_passjet_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
					}
				if (passJetSelection && passMetSelection && passBtagsSelection)
					{
					// pre-tau selection
					// control fake rates (might be useful)
					record_jets_fakerate_distrs(string("singlemu_"), string("control_pretau"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					// and the region of our selection fake rates
					// (used for jet origins mainly, plus for control)
					record_jets_fakerate_distrs(string("singlemu_"), string("pretau"), selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);

					fill_2d(string("control_met_singlemu_pretau_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singlemu_pretau_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singlemu_pretau_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singlemu_pretau_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);

					fill_1d(string("singlemu_pretauselection_nRawJets"),        10, 0,10,  jets.size(),              weight);
					fill_1d(string("singlemu_pretauselection_njets"),           10, 0,10,  selJets.size(),           weight);
					fill_1d(string("singlemu_pretauselection_njetsNoLep"),      10, 0,10,  selJetsNoLep.size(),      weight);
					fill_1d(string("singlemu_pretauselection_njetsNoLepNoTau"), 10, 0,10,  selJetsNoLepNoTau.size(), weight);
					fill_1d(string("singlemu_pretauselection_nRawTaus"),        10, 0,10,  taus.size(),              weight);
					fill_1d(string("singlemu_pretauselection_ntaus"),           10, 0,10,  selTaus.size(),           weight);
					fill_1d(string("singlemu_pretauselection_ntausNoLep"),      10, 0,10,  selTausNoLep.size(),      weight);
					fill_1d(string("singlemu_pretauselection_ntausNoLepNoJet"), 10, 0,10,  selTausNoLepNoJet.size(), weight);

					// calculate jet-to-tau fake rate per all jets and save the sum
					double jet_to_tau_no_fake_prob = 1.0;
					double jet_to_tau_no_fake_prob1_q = 1.0; // done with only histo 1
					double jet_to_tau_no_fake_prob1_w = 1.0; // done with only histo 1
					double jet_to_tau_no_fake_prob2_q = 1.0; // only histo 2
					double jet_to_tau_no_fake_prob2_w = 1.0; // only histo 2

					double jet_to_tau_fake_rate = 1.0;
					double jet_to_tau_fake_rate1_q = 1.0; // done with only histo 1
					double jet_to_tau_fake_rate1_w = 1.0; // done with only histo 1
					double jet_to_tau_fake_rate2_q = 1.0; // only histo 2
					double jet_to_tau_fake_rate2_w = 1.0; // only histo 2

					double jet_to_tau_no_fake_rate_elmu = 1.0;
					double jet_to_tau_fake_rate_elmu = 1.0;

					double jet_to_tau_no_fake_rate_mumu = 1.0;
					double jet_to_tau_fake_rate_mumu = 1.0;

					// using selJetsNoLep jets
					if (debug)
						{
						cout << "(mu) passed pre-tau selection" << "\n";
						cout << "N NoLep jets for pre-tau fakerate = " << selJetsNoLep.size() << "\n";
						cout << "N NoLepNoTau jets for pre-tau fakerate = " << selJetsNoLepNoTau.size() << "\n";

						cout << "FakeRates of NoLepNoTau jets:\n";
						for(size_t n=0; n<selJetsNoLepNoTau.size(); ++n)
							{
							cout << (1. - jetToTauFakeRate_vanila(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]), debug));
							cout << " ";
							cout << (1. - jetToTauFakeRate_vanila(tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]), debug));
							cout << "\n";
							}
						}

					//for(size_t n=0; n<selJetsNoLep.size(); ++n)
					// new fake rates are calculated
					// from all 20GeV jets
					// and all 30GeV, if there are not 2 of them... (yeah..)
					for(size_t n=0; n<selJetsNoLep.size(); ++n)
						// TODO: selJetsNoLepNoTau ???
						{
						pat::Jet& jet = selJetsNoLep[n];
						if (isMC)
							{
							int partID = jet.partonFlavour();
							//int hadrID = jet.hadronFlavour(); // TODO: test later
							//const reco::GenParticle & jet_origin = jet.genParton();
							//const reco::GenParticle* jet_origin = jet.genParton();
							//const reco::LeafCandidate* jet_origin = jet.genParton();
							//cout << "BBB " << abs( jet_origin.pdgId() ) << "\n";
							//fill_particle_ids(string("smu_pretau_jet_origin_ids"), 1, weight);

							fill_1d(string("smu_pretau_jet_origin_ids"), 100, 0, 100,  partID, weight);
							}

						// FAKERATES
						/*
						//jet_to_tau_no_fake_prob  *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob1_q *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob1_w *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_q *= (1. - jetToTauFakeRate(tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_w *= (1. - jetToTauFakeRate(tau_fake_rate2_jets_histo_w, tau_fake_rate2_taus_histo_w, tau_fake_rate2_jets_histo_w, tau_fake_rate2_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));

						jet_to_tau_no_fake_rate_elmu *= (1. - jetToTauFakeRate(tau_fake_rate_jets_histo_elmu, tau_fake_rate_taus_histo_elmu, tau_fake_rate_jets_histo_elmu, tau_fake_rate_taus_histo_elmu, 1.0, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_mumu *= (1. - jetToTauFakeRate(tau_fake_rate_jets_histo_mumu, tau_fake_rate_taus_histo_mumu, tau_fake_rate_jets_histo_mumu, tau_fake_rate_taus_histo_mumu, 1.0, jet.pt(), jet.eta(), jet_radius(jet), debug));
						*/

						//double jetToTauFakeRate_Projections(
						//		FakeRateProjections & tau_fake_rate_jets_histo,
						//		FakeRateProjections & tau_fake_rate_taus_histo,
						//		Double_t tau_fake_rate_histo1_fraction,
						//		Double_t jet_pt, Double_t jet_eta, Double_t jet_radius,
						//		bool debug)

						jet_to_tau_no_fake_prob1_q *= (1. - jetToTauFakeRate_Projections(frates_qcd_jets_proj, frates_qcd_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_w *= (1. - jetToTauFakeRate_Projections(frates_wjets_jets_proj, frates_wjets_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_elmu *= (1. - jetToTauFakeRate_Projections(frates_elmu_jets_proj, frates_elmu_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_mumu *= (1. - jetToTauFakeRate_Projections(frates_mumu_jets_proj, frates_mumu_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						}



					jet_to_tau_fake_rate  = 1.0 - jet_to_tau_no_fake_prob;
					jet_to_tau_fake_rate1_q = 1.0 - jet_to_tau_no_fake_prob1_q; // done with only histo 1
					jet_to_tau_fake_rate1_w = 1.0 - jet_to_tau_no_fake_prob1_w; // done with only histo 1
					jet_to_tau_fake_rate2_q = 1.0 - jet_to_tau_no_fake_prob2_q; // only histo 2
					jet_to_tau_fake_rate2_w = 1.0 - jet_to_tau_no_fake_prob2_w; // only histo 2

					jet_to_tau_fake_rate_elmu = 1.0 - jet_to_tau_no_fake_rate_elmu;
					jet_to_tau_fake_rate_mumu = 1.0 - jet_to_tau_no_fake_rate_mumu;

					if (debug)
						{
						cout << "no-fake probs: " << jet_to_tau_no_fake_prob1_q << " " << jet_to_tau_no_fake_prob << " " << jet_to_tau_no_fake_prob2_q << "\n";
						cout << "fakerates: " << jet_to_tau_fake_rate1_q << " " << jet_to_tau_fake_rate << " " << jet_to_tau_fake_rate2_q << "\n";
						}

					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 1,  weight * (jet_to_tau_fake_rate  < 1. ? jet_to_tau_fake_rate  : 1.));
					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 2,  weight * (jet_to_tau_fake_rate1_q  < 1. ? jet_to_tau_fake_rate1_q  : 1.));
					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 3,  weight * (jet_to_tau_fake_rate1_w  < 1. ? jet_to_tau_fake_rate1_w  : 1.));
					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 4,  weight * (jet_to_tau_fake_rate2_q  < 1. ? jet_to_tau_fake_rate2_q  : 1.));
					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 5,  weight * (jet_to_tau_fake_rate2_w  < 1. ? jet_to_tau_fake_rate2_w  : 1.));

					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 6,  weight * (jet_to_tau_fake_rate_elmu  < 1. ? jet_to_tau_fake_rate_elmu  : 1.));
					fill_1d(string("singlemu_pretauselection_jettotaufakerates"), 10, 0,10, 7,  weight * (jet_to_tau_fake_rate_mumu  < 1. ? jet_to_tau_fake_rate_mumu  : 1.));
					}

				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection) {
					// post-tau selection
					fill_1d(string("smu_passtau_selection_njets"), 20, 0, 20,   jets.size(), weight);
					fill_1d(string("smu_passtau_selection_nselJets"), 20, 0, 20,   selJets.size(), weight);
					fill_1d(string("smu_passtau_selection_nselJetsNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
					fill_1d(string("smu_passtau_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

					fill_1d(string("smu_passtau_selection_ntaus"), 20, 0, 20,   taus.size(), weight);
					fill_1d(string("smu_passtau_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
					fill_1d(string("smu_passtau_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
					fill_1d(string("smu_passtau_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

					fill_1d( string("weightflow_control_mu_passtau"),  300, -3, 3, weight, weight);

					fill_2d(string("control_met_singlemu_passtau_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singlemu_passtau_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singlemu_passtau_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singlemu_passtau_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
					}


				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS) {
					fill_1d(string("smu_passos_selection_ntaus"), 20, 0, 20,   taus.size(), weight);
					fill_1d(string("smu_passos_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
					fill_1d(string("smu_passos_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
					fill_1d(string("smu_passos_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

					fill_1d( string("singlemu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("singlemu_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("singlemu_selection_njets"), 10, 0, 10, selJetsNoLepNoTau.size(), weight);
					fill_1d( string("singlemu_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);

					// pts
					fill_1d( string("singlemu_selection_muon_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("singlemu_selection_tau_pt"),  200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("singlemu_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[0].pt(), weight);
					fill_1d( string("singlemu_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[1].pt(), weight);
					fill_1d( string("singlemu_selection_bjet_pt"), 200, 0, 400, selBJets[0].pt(), weight);
					// fill_1d( string("singlemu_selection_met_pt"), 200, 0, 400, n_met.pt(), weight);
					fill_1d( string("singlemu_selection_met_pt"),  200, 0, 400, met.pt(), weight);

					// energy
					fill_1d( string("singlemu_selection_muon_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("singlemu_selection_tau_energy"),  200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("singlemu_selection_jet1_energy"), 200, 0, 400, selJetsNoLepNoTau[0].energy(), weight);
					fill_1d( string("singlemu_selection_jet2_energy"), 200, 0, 400, selJetsNoLepNoTau[1].energy(), weight);
					fill_1d( string("singlemu_selection_bjet_energy"), 200, 0, 400, selBJets[0].energy(), weight);
					// fill_1d( string("singlemu_selection_met_energy"), 200, 0, 400,  n_met.energy(), weight);
					fill_1d( string("singlemu_selection_met_energy"),  200, 0, 400, met.energy(), weight);

					// eta
					fill_1d( string("singlemu_selection_muon_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("singlemu_selection_tau_eta"),  300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("singlemu_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[0].eta(), weight);
					fill_1d( string("singlemu_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[1].eta(), weight);
					fill_1d( string("singlemu_selection_bjet_eta"), 300, -3, 3, selBJets[0].eta(), weight);
					// fill_1d( string("singlemu_selection_met_eta"), 300, -3, 3, n_met.eta(), weight);
					fill_1d( string("singlemu_selection_met_eta"),  300, -3, 3, met.eta(), weight);

					weightflow_control_mu_selection += weight;
					fill_1d( string("weightflow_control_mu_selection"),  300, -3, 3, weight, weight);
					}

				// the 5 geometrical for loops...
				// TODO: it's awekward

				bool s1 = false;
				bool s2 = false;
				bool s3 = false;
				bool s4 = false;
				bool s5 = false;
				for (int i1 = 0; i1 < 2; i1++) // consider selection1 or not
					{
					s1 = !s1;
					for (int i2 = 0; i2 < 2; i2++)
					{
					s2 = !s2;
					for (int i3 = 0; i3 < 2; i3++)
					{
					s3 = !s3;
					for (int i4 = 0; i4 < 2; i4++)
					{
					s4 = !s4;
					for (int i5 = 0; i5 < 2; i5++)
						{
						s5 = !s5;
						// if a selection is not considered
						// it is passed
						bool pass = (s1 ? passJetSelection : true);
						pass &= (s2 ? passMetSelection : true);
						pass &= (s3 ? passBtagsSelection : true);
						pass &= (s4 ? passTauSelection : true);
						pass &= (s5 ? passOS : true);

						if (!pass) continue;

						unsigned int multi = 0;
						multi += (s1 ? 1 : 0);
						multi += (s2 ? 2 : 0);
						multi += (s3 ? 4 : 0);
						multi += (s4 ? 8 : 0);
						multi += (s5 ? 16 : 0);

						//fill_1d(string("tauIDSFs_in_weightflow_mu_") + to_string(multi), 200, 0., 2.,   weight_tauIDsf, 1);
						// + to_string(multi)
						fill_1d(string("weightflow_mu"), 300, 0, 300,   31 + multi, weight );
						fill_1d(string("weightflow_mu_without_tauIDSFs"), 300, 0, 300,   31 + multi, weight_without_tauIDsf);
						fill_1d(string("eventflow_mu"), 300, 0, 300,   31 + multi, 1);

						fill_1i(string("weightflow_mu_I"), 300, 0, 300,   31 + multi, weight );
						fill_1i(string("weightflow_mu_without_tauIDSFs_I"), 300, 0, 300,   31 + multi, weight_without_tauIDsf);
						fill_1i(string("eventflow_mu_I"), 300, 0, 300,   31 + multi, 1);
						}
					}
					}
					}
					}
				}

			if (isSingleE)
				{
				fill_2d(string("control_met_singleel_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
				for (int i=0; i<selJetsNoLep.size(); i++)
					fill_2d(string("control_jet_singleel_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
				for (int i=0; i<selTausNoLep.size(); i++)
					fill_2d(string("control_tau_singleel_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
				for (int i=0; i<selLeptons.size(); i++)
					fill_2d(string("control_lep_singleel_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);

				if (passJetSelection) {
					record_jets_fakerate_distrs(string("singleel_"), string("control_passjet"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					fill_2d(string("control_met_singleel_passjet_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singleel_passjet_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singleel_passjet_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singleel_passjet_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
					}

				//
				// the 5 geometrical for loops...

				bool s1 = false;
				bool s2 = false;
				bool s3 = false;
				bool s4 = false;
				bool s5 = false;
				for (int i1 = 0; i1 < 2; i1++) // consider selection1 or not
					{
					s1 = !s1;
					for (int i2 = 0; i2 < 2; i2++)
					{
					s2 = !s2;
					for (int i3 = 0; i3 < 2; i3++)
					{
					s3 = !s3;
					for (int i4 = 0; i4 < 2; i4++)
					{
					s4 = !s4;
					for (int i5 = 0; i5 < 2; i5++)
						{
						s5 = !s5;
						// if a selection is not considered
						// it is passed
						bool pass = (s1 ? passJetSelection : true);
						pass &= (s2 ? passMetSelection : true);
						pass &= (s3 ? passBtagsSelection : true);
						pass &= (s4 ? passTauSelection : true);
						pass &= (s5 ? passOS : true);

						if (!pass) continue;

						unsigned int multi = 0;
						multi += (s1 ? 1 : 0);
						multi += (s2 ? 2 : 0);
						multi += (s3 ? 4 : 0);
						multi += (s4 ? 8 : 0);
						multi += (s5 ? 16 : 0);

						//fill_1d(string("tauIDSFs_in_weightflow_el_") + to_string(multi), 200, 0., 2.,   weight_tauIDsf, 1);
						// + to_string(multi)
						fill_1d(string("weightflow_el"), 300, 0, 300,   31 + multi, weight );
						fill_1d(string("weightflow_el_without_tauIDSFs"), 300, 0, 300,   31 + multi, weight_without_tauIDsf);
						fill_1d(string("eventflow_el"), 300, 0, 300,   31 + multi, 1);

						fill_1i(string("weightflow_el_I"), 300, 0, 300,   31 + multi, weight );
						fill_1i(string("weightflow_el_without_tauIDSFs_I"), 300, 0, 300,   31 + multi, weight_without_tauIDsf);
						fill_1i(string("eventflow_el_I"), 300, 0, 300,   31 + multi, 1);
						}
					}
					}
					}
					}

				if (passJetSelection && passMetSelection && passBtagsSelection)
					{
					// pre-tau selection
					//int record_jets_fakerate_distrs(string & channel, string & selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, double event_weight, bool isMC)
					// control
					record_jets_fakerate_distrs(string("singleel_"), string("control_pretau"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					// and our selection fake rates (counting jet origins)
					record_jets_fakerate_distrs(string("singleel_"), string("pretau"), selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);

					fill_2d(string("control_met_singleel_pretau_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singleel_pretau_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singleel_pretau_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singleel_pretau_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);

					//
					//fill_1d( string("singleel_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);
					fill_1d(string("singleel_pretauselection_nRawJets"),        10, 0,10,  jets.size(),              weight);
					fill_1d(string("singleel_pretauselection_njets"),           10, 0,10,  selJets.size(),           weight);
					fill_1d(string("singleel_pretauselection_njetsNoLep"),      10, 0,10,  selJetsNoLep.size(),      weight);
					fill_1d(string("singleel_pretauselection_njetsNoLepNoTau"), 10, 0,10,  selJetsNoLepNoTau.size(), weight);
					fill_1d(string("singleel_pretauselection_nRawTaus"),        10, 0,10,  taus.size(),              weight);
					fill_1d(string("singleel_pretauselection_ntaus"),           10, 0,10,  selTaus.size(),           weight);
					fill_1d(string("singleel_pretauselection_ntausNoLep"),      10, 0,10,  selTausNoLep.size(),      weight);
					fill_1d(string("singleel_pretauselection_ntausNoLepNoJet"), 10, 0,10,  selTausNoLepNoJet.size(), weight);

					// calculate jet-to-tau fake rate per all jets and save the sum
					double jet_to_tau_no_fake_prob = 1.0;
					double jet_to_tau_no_fake_prob1_q = 1.0; // done with only histo 1
					double jet_to_tau_no_fake_prob1_w = 1.0; // done with only histo 1
					double jet_to_tau_no_fake_prob2_q = 1.0; // only histo 2
					double jet_to_tau_no_fake_prob2_w = 1.0; // only histo 2

					double jet_to_tau_fake_rate = 1.0;
					double jet_to_tau_fake_rate1_q = 1.0; // done with only histo 1
					double jet_to_tau_fake_rate1_w = 1.0; // done with only histo 1
					double jet_to_tau_fake_rate2_q = 1.0; // only histo 2
					double jet_to_tau_fake_rate2_w = 1.0; // only histo 2

					double jet_to_tau_no_fake_rate_elmu = 1.0;
					double jet_to_tau_fake_rate_elmu = 1.0;

					double jet_to_tau_no_fake_rate_mumu = 1.0;
					double jet_to_tau_fake_rate_mumu = 1.0;

					// using selJetsNoLep jets
					if (debug)
						{
						cout << "(el) passed pre-tau selection" << "\n";
						cout << "N jets for pre-tau fakerate = " << selJetsNoLep.size() << "\n";
						cout << "N NoLepNoTau jets for pre-tau fakerate = " << selJetsNoLepNoTau.size() << "\n";

						cout << "FakeRates of NoLepNoTau jets:\n";
						for(size_t n=0; n<selJetsNoLepNoTau.size(); ++n)
							{
							cout << (1. - jetToTauFakeRate_vanila(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]), debug));
							cout << " ";
							cout << (1. - jetToTauFakeRate_vanila(tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate_histo1_fraction, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]), debug));
							cout << "\n";
							}
						}

					//for(size_t n=0; n<selJetsNoLep.size(); ++n)
					for(size_t n=0; n<selJetsNoLep.size(); ++n)
						{
						pat::Jet& jet = selJetsNoLep[n];
						if (debug) cout << n << ":\n";
						// jet_to_tau_fake_rate += jetToTauFakeRate(tau_fake_rate_jets_histo, tau_fake_rate_taus_histo, selJetsNoLep[n].pt(), selJetsNoLep[n].eta(), jet_radius(selJetsNoLep[n]));
						if (isMC)
							{
							int partID = jet.partonFlavour();
							//fill_particle_ids(string("sel_pretau_jet_origin_ids"), selJetsNoLep[n].genParton()->pdgId(), weight);
							//fill_particle_ids(string("sel_pretau_jet_origin_ids"), partID, weight);
							fill_1d(string("sel_pretau_jet_origin_ids"), 100, 0, 100,  partID, weight);
							}

						/*
						jet_to_tau_no_fake_prob  *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob1_q *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob1_w *= (1. - jetToTauFakeRate(tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate1_jets_histo_w, tau_fake_rate1_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_q *= (1. - jetToTauFakeRate(tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate2_jets_histo_q, tau_fake_rate2_taus_histo_q, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_w *= (1. - jetToTauFakeRate(tau_fake_rate2_jets_histo_w, tau_fake_rate2_taus_histo_w, tau_fake_rate2_jets_histo_w, tau_fake_rate2_taus_histo_w, tau_fake_rate_histo1_fraction, jet.pt(), jet.eta(), jet_radius(jet), debug));

						jet_to_tau_no_fake_rate_elmu *= (1. - jetToTauFakeRate(tau_fake_rate_jets_histo_elmu, tau_fake_rate_taus_histo_elmu, tau_fake_rate_jets_histo_elmu, tau_fake_rate_taus_histo_elmu, 1.0, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_mumu *= (1. - jetToTauFakeRate(tau_fake_rate_jets_histo_mumu, tau_fake_rate_taus_histo_mumu, tau_fake_rate_jets_histo_mumu, tau_fake_rate_taus_histo_mumu, 1.0, jet.pt(), jet.eta(), jet_radius(jet), debug));
						*/

						jet_to_tau_no_fake_prob1_q *= (1. - jetToTauFakeRate_Projections(frates_qcd_jets_proj, frates_qcd_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_prob2_w *= (1. - jetToTauFakeRate_Projections(frates_wjets_jets_proj, frates_wjets_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_elmu *= (1. - jetToTauFakeRate_Projections(frates_elmu_jets_proj, frates_elmu_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						jet_to_tau_no_fake_rate_mumu *= (1. - jetToTauFakeRate_Projections(frates_mumu_jets_proj, frates_mumu_taus_proj, 1, jet.pt(), jet.eta(), jet_radius(jet), debug));
						}



					jet_to_tau_fake_rate  = 1.0 - jet_to_tau_no_fake_prob;
					jet_to_tau_fake_rate1_q = 1.0 - jet_to_tau_no_fake_prob1_q; // done with only histo 1
					jet_to_tau_fake_rate1_w = 1.0 - jet_to_tau_no_fake_prob1_w; // done with only histo 1
					jet_to_tau_fake_rate2_q = 1.0 - jet_to_tau_no_fake_prob2_q; // only histo 2
					jet_to_tau_fake_rate2_w = 1.0 - jet_to_tau_no_fake_prob2_w; // only histo 2

					jet_to_tau_fake_rate_elmu = 1.0 - jet_to_tau_no_fake_rate_elmu;
					jet_to_tau_fake_rate_mumu = 1.0 - jet_to_tau_no_fake_rate_mumu;

					if (debug)
						{
						cout << "fakerates: " << jet_to_tau_fake_rate1_q << " " << jet_to_tau_fake_rate << " " << jet_to_tau_fake_rate2_q << "\n";
						}

					/*
					increment(string("singleel_pretauselection_jettotaufakerate"),  weight * (jet_to_tau_fake_rate  < 1. ? jet_to_tau_fake_rate  : 1.));
					increment(string("singleel_pretauselection_jettotaufakerate1_q"), weight * (jet_to_tau_fake_rate1_q < 1. ? jet_to_tau_fake_rate1_q : 1.));
					increment(string("singleel_pretauselection_jettotaufakerate1_w"), weight * (jet_to_tau_fake_rate1_w < 1. ? jet_to_tau_fake_rate1_w : 1.));
					increment(string("singleel_pretauselection_jettotaufakerate2_q"), weight * (jet_to_tau_fake_rate2_q < 1. ? jet_to_tau_fake_rate2_q : 1.));
					increment(string("singleel_pretauselection_jettotaufakerate2_w"), weight * (jet_to_tau_fake_rate2_w < 1. ? jet_to_tau_fake_rate2_w : 1.));
					*/
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 1,  weight * (jet_to_tau_fake_rate  < 1. ? jet_to_tau_fake_rate  : 1.));
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 2,  weight * (jet_to_tau_fake_rate1_q  < 1. ? jet_to_tau_fake_rate1_q  : 1.));
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 3,  weight * (jet_to_tau_fake_rate1_w  < 1. ? jet_to_tau_fake_rate1_w  : 1.));
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 4,  weight * (jet_to_tau_fake_rate2_q  < 1. ? jet_to_tau_fake_rate2_q  : 1.));
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 5,  weight * (jet_to_tau_fake_rate2_w  < 1. ? jet_to_tau_fake_rate2_w  : 1.));

					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 6,  weight * (jet_to_tau_fake_rate_elmu  < 1. ? jet_to_tau_fake_rate_elmu  : 1.));
					fill_1d(string("singleel_pretauselection_jettotaufakerates"), 10, 0,10, 7,  weight * (jet_to_tau_fake_rate_mumu  < 1. ? jet_to_tau_fake_rate_mumu  : 1.));
					}

				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection) {
					// post-tau selection
					fill_1d(string("sel_passtau_selection_njets"), 20, 0, 20,   jets.size(), weight);
					fill_1d(string("sel_passtau_selection_nselJets"), 20, 0, 20,   selJets.size(), weight);
					fill_1d(string("sel_passtau_selection_nselJetsNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
					fill_1d(string("sel_passtau_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

					fill_1d(string("sel_passtau_selection_ntaus"), 20, 0, 20,   taus.size(), weight);
					fill_1d(string("sel_passtau_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
					fill_1d(string("sel_passtau_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
					fill_1d(string("sel_passtau_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

					fill_1d( string("weightflow_control_el_passtau"),  300, -3, 3, weight, weight);

					fill_2d(string("control_met_singleel_passtau_pt_eta"), 250, 0., 500., 200, -4., 4., met.pt(), met.eta(), weight);
					for (int i=0; i<selJetsNoLep.size(); i++)
						fill_2d(string("control_jet_singleel_passtau_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selJetsNoLep[i].pt(), selJetsNoLep[i].eta(), weight);
					for (int i=0; i<selTausNoLep.size(); i++)
						fill_2d(string("control_tau_singleel_passtau_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., selTausNoLep[i].pt(), selTausNoLep[i].eta(), weight);
					for (int i=0; i<selLeptons.size(); i++)
						fill_2d(string("control_lep_singleel_passtau_selLeptons_pt_eta"), 250, 0., 500., 200, -4., 4., selLeptons[i].pt(), selLeptons[i].eta(), weight);
					}

				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS) {
					fill_1d(string("sel_passos_selection_ntaus"), 20, 0, 20,   taus.size(), weight);
					fill_1d(string("sel_passos_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
					fill_1d(string("sel_passos_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
					fill_1d(string("sel_passos_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

					fill_1d( string("singleel_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("singleel_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("singleel_selection_njets"), 10, 0, 10, selJetsNoLepNoTau.size(), weight);
					fill_1d( string("singleel_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);

					// pts
					fill_1d( string("singleel_selection_electron_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("singleel_selection_tau_pt"),  200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("singleel_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[0].pt(), weight);
					fill_1d( string("singleel_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[1].pt(), weight);
					fill_1d( string("singleel_selection_bjet_pt"), 200, 0, 400, selBJets[0].pt(), weight);
					// fill_1d( string("singleel_selection_met_pt"), 200, 0, 400, n_met.pt(), weight);
					fill_1d( string("singleel_selection_met_pt"),  200, 0, 400, met.pt(), weight);

					// energy
					fill_1d( string("singleel_selection_electron_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("singleel_selection_tau_energy"),  200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("singleel_selection_jet1_energy"), 200, 0, 400, selJetsNoLepNoTau[0].energy(), weight);
					fill_1d( string("singleel_selection_jet2_energy"), 200, 0, 400, selJetsNoLepNoTau[1].energy(), weight);
					fill_1d( string("singleel_selection_bjet_energy"), 200, 0, 400, selBJets[0].energy(), weight);
					// fill_1d( string("singleel_selection_met_energy"), 200, 0, 400, n_met.energy(), weight);
					fill_1d( string("singleel_selection_met_energy"),  200, 0, 400, met.energy(), weight);

					// eta
					fill_1d( string("singleel_selection_electron_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("singleel_selection_tau_eta"),  300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("singleel_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[0].eta(), weight);
					fill_1d( string("singleel_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[1].eta(), weight);
					fill_1d( string("singleel_selection_bjet_eta"), 300, -3, 3, selBJets[0].eta(), weight);
					// fill_1d( string("singleel_selection_met_eta"), 300, -3, 3, n_met.eta(), weight);
					fill_1d( string("singleel_selection_met_eta"),  300, -3, 3, met.eta(), weight);

					weightflow_control_el_selection += weight;
					fill_1d( string("weightflow_control_el_selection"),  300, -3, 3, weight, weight);
					}
				}
			}


		// --------------------------------------------- DILEPTON CHANNELS

		if (isDoubleE || isDoubleMu || isEMu)
			{
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
			bool passMetSelection(met.pt()>40.);
			bool passOS(selLeptons[0].pdgId() * selLeptons[1].pdgId() < 0 );
			// bool passBtagsSelection(selBJets.size()>1);
			bool passBtagsSelection(n_bjets>0);

			// MULTISELECT
			unsigned int multisel = 0;
			multisel += (passMllVeto ? 1 : 0);
			multisel += (passJetSelection ? 2 : 0);
			multisel += (passMetSelection ? 4 : 0);
			multisel += (passOS ? 8 : 0);
			multisel += (passBtagsSelection ? 16 : 0);

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
					}

				if (passJetSelection && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("elel_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("elel_"), string("passbtagfinal"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					//fill_1d( string("elel_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("elel_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					//fill_1d( string("elel_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("elel_selection_njets"), 10, 0, 10, selJetsNoLepNoTau.size(), weight);
					fill_1d( string("elel_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);

					// pt
					fill_1d( string("elel_selection_el1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("elel_selection_el2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("elel_selection_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("elel_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[0].pt(), weight);
					fill_1d( string("elel_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[1].pt(), weight);
					fill_1d( string("elel_selection_bjet_pt"), 200, 0, 400, selBJets[0].pt(), weight);
					// fill_1d( string("elel_selection_met_pt"), 200, 0, 400, n_met.pt(), weight);
					fill_1d( string("elel_selection_met_pt"), 200, 0, 400, met.pt(), weight);

					// energies
					fill_1d( string("elel_selection_el1_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("elel_selection_el2_energy"), 200, 0, 400, selLeptons[1].energy(), weight);
					// fill_1d( string("elel_selection_tau_energy"), 200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("elel_selection_jet1_energy"), 200, 0, 400, selJetsNoLepNoTau[0].energy(), weight);
					fill_1d( string("elel_selection_jet2_energy"), 200, 0, 400, selJetsNoLepNoTau[1].energy(), weight);
					fill_1d( string("elel_selection_bjet_energy"), 200, 0, 400, selBJets[0].energy(), weight);
					// fill_1d( string("elel_selection_met_energy"), 200, 0, 400, n_met.energy(), weight);
					fill_1d( string("elel_selection_met_energy"), 200, 0, 400, met.energy(), weight);

					//etas:
					fill_1d( string("elel_selection_el1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("elel_selection_el2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("elel_selection_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("elel_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[0].eta(), weight);
					fill_1d( string("elel_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[1].eta(), weight);
					fill_1d( string("elel_selection_bjet_eta"), 300, -3, 3, selBJets[0].eta(), weight);
					// fill_1d( string("elel_selection_met_eta"), 300, -3, 3, n_met.eta(), weight);
					fill_1d( string("elel_selection_met_eta"), 300, -3, 3, met.eta(), weight);

					weightflow_control_elel_selection += weight;
					fill_1d( string("weightflow_control_elel_selection"),  300, -3, 3, weight, weight);
					}

				// the 5 geometrical for loops...
				// TODO: it's awekward

				bool s1 = false;
				bool s2 = false;
				bool s3 = false;
				bool s4 = false;
				bool s5 = false;
				for (int i1 = 0; i1 < 2; i1++) // consider selection1 or not
					{
					s1 = !s1;
					for (int i2 = 0; i2 < 2; i2++)
					{
					s2 = !s2;
					for (int i3 = 0; i3 < 2; i3++)
					{
					s3 = !s3;
					for (int i4 = 0; i4 < 2; i4++)
					{
					s4 = !s4;
					for (int i5 = 0; i5 < 2; i5++)
						{
						s5 = !s5;
						// if a selection is not considered
						// it is passed
						bool pass = (s1 ? passMllVeto : true);
						pass &= (s2 ? passJetSelection : true);
						pass &= (s3 ? passMetSelection : true);
						pass &= (s4 ? passOS : true);
						pass &= (s5 ? passBtagsSelection : true);

						if (!pass) continue;

						unsigned int multi = 0;
						multi += (s1 ? 1 : 0);
						multi += (s2 ? 2 : 0);
						multi += (s3 ? 4 : 0);
						multi += (s4 ? 8 : 0);
						multi += (s5 ? 16 : 0);

						// + to_string(multi)
						fill_1d(string("weightflow_elel"), 300, 0, 300,   31 + multi, weight);
						fill_1d(string("eventflow_elel"), 300, 0, 300,   31 + multi, 1);
						}
					}
					}
					}
					}
				}

			if (isDoubleMu)
				{


				if (passJetSelection)
					{
					record_jets_fakerate_distrs(string("mumu_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passJetSelection && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("mumu_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("mumu_"), string("passbtagfinal"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					//fill_1d( string("mumu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					fill_1d( string("mumu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					//fill_1d( string("mumu_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("mumu_selection_njets"), 10, 0, 10, selJetsNoLepNoTau.size(), weight);
					fill_1d( string("mumu_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);

					// pt
					fill_1d( string("mumu_selection_mu1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("mumu_selection_mu2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("mumu_selection_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("mumu_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[0].pt(), weight);
					fill_1d( string("mumu_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[1].pt(), weight);
					fill_1d( string("mumu_selection_bjet_pt"), 200, 0, 400, selBJets[0].pt(), weight);
					// fill_1d( string("mumu_selection_met_pt"), 200, 0, 400, n_met.pt(), weight);
					fill_1d( string("mumu_selection_met_pt"), 200, 0, 400, met.pt(), weight);

					// energies
					fill_1d( string("mumu_selection_mu1_energy"), 200, 0, 200, selLeptons[0].energy(), weight);
					fill_1d( string("mumu_selection_mu2_energy"), 200, 0, 200, selLeptons[1].energy(), weight);
					// fill_1d( string("mumu_selection_tau_energy"), 200, 0, 200, selTausNoLep[0].energy(), weight);
					fill_1d( string("mumu_selection_jet1_energy"), 200, 0, 200, selJetsNoLepNoTau[0].energy(), weight);
					fill_1d( string("mumu_selection_jet2_energy"), 200, 0, 200, selJetsNoLepNoTau[1].energy(), weight);
					fill_1d( string("mumu_selection_bjet_energy"), 200, 0, 200, selBJets[0].energy(), weight);
					// fill_1d( string("mumu_selection_met_energy"), 200, 0, 200, n_met.energy(), weight);
					fill_1d( string("mumu_selection_met_energy"), 200, 0, 200, met.energy(), weight);

					//etas:
					fill_1d( string("mumu_selection_mu1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("mumu_selection_mu2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("mumu_selection_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("mumu_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[0].eta(), weight);
					fill_1d( string("mumu_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[1].eta(), weight);
					fill_1d( string("mumu_selection_bjet_eta"), 300, -3, 3, selBJets[0].eta(), weight);
					// fill_1d( string("mumu_selection_met_eta"), 300, -3, 3, n_met.eta(), weight);
					fill_1d( string("mumu_selection_met_eta"), 300, -3, 3, met.eta(), weight);

					weightflow_control_mumu_selection += weight;
					fill_1d( string("weightflow_control_mumu_selection"),  300, -3, 3, weight, weight);
					}

				// the 5 geometrical for loops...
				// TODO: it's awekward

				bool s1 = false;
				bool s2 = false;
				bool s3 = false;
				bool s4 = false;
				bool s5 = false;
				for (int i1 = 0; i1 < 2; i1++) // consider selection1 or not
					{
					s1 = !s1;
					for (int i2 = 0; i2 < 2; i2++)
					{
					s2 = !s2;
					for (int i3 = 0; i3 < 2; i3++)
					{
					s3 = !s3;
					for (int i4 = 0; i4 < 2; i4++)
					{
					s4 = !s4;
					for (int i5 = 0; i5 < 2; i5++)
						{
						s5 = !s5;
						// if a selection is not considered
						// it is passed
						bool pass = (s1 ? passMllVeto : true);
						pass &= (s2 ? passJetSelection : true);
						pass &= (s3 ? passMetSelection : true);
						pass &= (s4 ? passOS : true);
						pass &= (s5 ? passBtagsSelection : true);

						if (!pass) continue;

						unsigned int multi = 0;
						multi += (s1 ? 1 : 0);
						multi += (s2 ? 2 : 0);
						multi += (s3 ? 4 : 0);
						multi += (s4 ? 8 : 0);
						multi += (s5 ? 16 : 0);

						// + to_string(multi)
						fill_1d(string("weightflow_mumu"), 300, 0, 300,   31 + multi, weight);
						fill_1d(string("eventflow_mumu"), 300, 0, 300,   31 + multi, 1);
						}
					}
					}
					}
					}
				}

			if (isEMu)
				{

				if (passJetSelection)
					{
					record_jets_fakerate_distrs(string("elmu_"), string("passjets"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passJetSelection && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("elmu_"), string("passjetsNbtag"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					record_jets_fakerate_distrs(string("elmu_"), string("passbtagfinal"), selJets_JetTauFakeRate_NoLep, selTaus_JetTauFakeRate_NoLep, visible_gen_taus, weight, isMC);

					fill_1d( string("elmu_selection_nleps"), 10, 0, 10, selLeptons.size(), weight);
					//fill_1d( string("elmu_selection_ntaus"), 10, 0, 10, selTausNoLep.size(), weight);
					fill_1d( string("elmu_selection_njets"), 10, 0, 10, selJetsNoLepNoTau.size(), weight);
					fill_1d( string("elmu_selection_nbjets"), 10, 0, 10, selBJets.size(), weight);

					// pt
					fill_1d( string("elmu_selection_lep1_pt"), 200, 0, 400, selLeptons[0].pt(), weight);
					fill_1d( string("elmu_selection_lep2_pt"), 200, 0, 400, selLeptons[1].pt(), weight);
					// fill_1d( string("elmu_selection_tau_pt"), 200, 0, 400, selTausNoLep[0].pt(), weight);
					fill_1d( string("elmu_selection_jet1_pt"), 200, 0, 400, selJetsNoLepNoTau[0].pt(), weight);
					fill_1d( string("elmu_selection_jet2_pt"), 200, 0, 400, selJetsNoLepNoTau[1].pt(), weight);
					fill_1d( string("elmu_selection_bjet_pt"), 200, 0, 400, selBJets[0].pt(), weight);
					// fill_1d( string("elmu_selection_met_pt"), 200, 0, 400, n_met.pt(), weight);
					fill_1d( string("elmu_selection_met_pt"), 200, 0, 400, met.pt(), weight);

					// energies
					fill_1d( string("elmu_selection_lep1_energy"), 200, 0, 400, selLeptons[0].energy(), weight);
					fill_1d( string("elmu_selection_lep2_energy"), 200, 0, 400, selLeptons[1].energy(), weight);
					// fill_1d( string("elmu_selection_tau_energy"), 200, 0, 400, selTausNoLep[0].energy(), weight);
					fill_1d( string("elmu_selection_jet1_energy"), 200, 0, 400, selJetsNoLepNoTau[0].energy(), weight);
					fill_1d( string("elmu_selection_jet2_energy"), 200, 0, 400, selJetsNoLepNoTau[1].energy(), weight);
					fill_1d( string("elmu_selection_bjet_energy"), 200, 0, 400, selBJets[0].energy(), weight);
					// fill_1d( string("elmu_selection_met_energy"), 200, 0, 400, n_met.energy(), weight);
					fill_1d( string("elmu_selection_met_energy"), 200, 0, 400, met.energy(), weight);

					//etas:
					fill_1d( string("elmu_selection_lep1_eta"), 300, -3, 3, selLeptons[0].eta(), weight);
					fill_1d( string("elmu_selection_lep2_eta"), 300, -3, 3, selLeptons[1].eta(), weight);
					// fill_1d( string("elmu_selection_tau_eta"), 300, -3, 3, selTausNoLep[0].eta(), weight);
					fill_1d( string("elmu_selection_jet1_eta"), 300, -3, 3, selJetsNoLepNoTau[0].eta(), weight);
					fill_1d( string("elmu_selection_jet2_eta"), 300, -3, 3, selJetsNoLepNoTau[1].eta(), weight);
					fill_1d( string("elmu_selection_bjet_eta"), 300, -3, 3, selBJets[0].eta(), weight);
					// fill_1d( string("elmu_selection_met_eta"), 300, -3, 3, n_met.eta(), weight);
					fill_1d( string("elmu_selection_met_eta"), 300, -3, 3, met.eta(), weight);

					weightflow_control_elmu_selection += weight;
					fill_1d( string("weightflow_control_elmu_selection"),  300, -3, 3, weight, weight);
					}

				// the 5 geometrical for loops...
				// TODO: it's awekward

				bool s1 = false;
				bool s2 = false;
				bool s3 = false;
				bool s4 = false;
				bool s5 = false;
				for (int i1 = 0; i1 < 2; i1++) // consider selection1 or not
					{
					s1 = !s1;
					for (int i2 = 0; i2 < 2; i2++)
					{
					s2 = !s2;
					for (int i3 = 0; i3 < 2; i3++)
					{
					s3 = !s3;
					for (int i4 = 0; i4 < 2; i4++)
					{
					s4 = !s4;
					for (int i5 = 0; i5 < 2; i5++)
						{
						s5 = !s5;
						// if a selection is not considered
						// it is passed
						bool pass = (s1 ? passMllVeto : true);
						pass &= (s2 ? passJetSelection : true);
						pass &= (s3 ? passMetSelection : true);
						pass &= (s4 ? passOS : true);
						pass &= (s5 ? passBtagsSelection : true);

						if (!pass) continue;

						unsigned int multi = 0;
						multi += (s1 ? 1 : 0);
						multi += (s2 ? 2 : 0);
						multi += (s3 ? 4 : 0);
						multi += (s4 ? 8 : 0);
						multi += (s5 ? 16 : 0);

						// + to_string(multi)
						fill_1d(string("weightflow_elmu"), 300, 0, 300,   31 + multi, weight);
						fill_1d(string("eventflow_elmu"), 300, 0, 300,   31 + multi, 1);
						}
					}
					}
					}
					}

				//increment(string("weightflow_emu_") + to_string(multisel), weight);
				// increment(string("weightflow_up_emu_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_emu_") + to_string(multisel), weight_down);

				/*
				fill_pt_e( string("doubleelmu_channel_met_pt"), met.pt(), weight);
				if (passMllVeto && passJetSelection)
					{
					fill_pt_e( string("doubleelmu_jetsel_met_pt"), met.pt(), weight);
					}
				if (passMllVeto && passJetSelection && passOS && passBtagsSelection)
					{
					fill_pt_e( string("doubleelmu_allbutmetsel_met_pt"), met.pt(), weight);
					}

				fill_n( string("elmu_channel_n_jets"), n_jets, weight);
				fill_n( string("elmu_channel_n_bjets"), n_bjets, weight);
				fill_n( string("elmu_channel_n_taus"), n_taus, weight);

				if (passMllVeto && passJetSelection && passMetSelection && passOS)
					{
					fill_n( string("elmu_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("elmu_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("elmu_prebselpoint_n_taus"), n_taus, weight);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					increment( string("weightflow_events_passed_doubleemu_selection"), 1 );

					increment( string("weightflow_weight_passed_doubleemu_selection"), weight );

					fill_n( string("elmu_selection_n_jets"), n_jets, weight);
					fill_n( string("elmu_selection_n_bjets"), n_bjets, weight);
					fill_n( string("elmu_selection_n_taus"), n_taus, weight);

					fill_pt_e( string("elmu_selection_el_pt"), selElectrons[0].pt(), weight);
					fill_pt_e( string("elmu_selection_mu_pt"), selMuons[0].pt(), weight);
					// fill_pt_e( string("elmu_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("elmu_selection_jet1_pt"), selJetsNoLep[0].pt(), weight);
					fill_pt_e( string("elmu_selection_jet2_pt"), selJetsNoLep[1].pt(), weight);
					fill_pt_e( string("elmu_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("elmu_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("elmu_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("elmu_selection_el_energy"), selElectrons[0].energy(), weight);
					fill_pt_e( string("elmu_selection_mu_energy"), selMuons[0].energy(), weight);
					// fill_pt_e( string("elmu_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("elmu_selection_jet1_energy"), selJetsNoLep[0].energy(), weight);
					fill_pt_e( string("elmu_selection_jet2_energy"), selJetsNoLep[1].energy(), weight);
					fill_pt_e( string("elmu_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("elmu_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("elmu_selection_met_energy"), met.energy(), weight);

					// etas
					fill_eta( string("elmu_selection_el_eta"), selElectrons[0].eta(), weight);
					fill_eta( string("elmu_selection_mu_eta"), selMuons[0].eta(), weight);
					// fill_eta( string("elmu_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("elmu_selection_jet1_eta"), selJetsNoLep[0].eta(), weight);
					fill_eta( string("elmu_selection_jet2_eta"), selJetsNoLep[1].eta(), weight);
					fill_eta( string("elmu_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("elmu_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("elmu_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_elmuselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_elmuselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_elmuselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_elmuselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_elmuselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_elmuselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_elmuselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_elmuselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				*/
				}
			}

		/* debugging
		if (multisel > MULTISEL_SIZE)
			{
			printf("in event %d too large multisel: %d", iev, multisel);
			break;
			}
		*/

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


if(saveSummaryTree)
	{
	TDirectory* cwd = gDirectory;
	summaryFile->cd();
	summaryTree->Write();
	summaryFile->Close();
	delete summaryFile;
	cwd->cd();
	}


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
	TFile* out_f = TFile::Open (outdir + TString(string("/") + dtag_s + channel + string("_") + job_num + string(".root")), "CREATE");
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

printf ("New output results saved in %s\n", (outdir.Data() + string("/") + dtag_s + string("_<channel>") + string("_") + job_num + string(".root")).c_str());


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


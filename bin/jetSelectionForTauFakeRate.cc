//
// Oleksii Toldaiev, <oleksii.toldaiev@gmail.com>
//
// ttbar to leptons event selection
// CMSSW version 80X (8_0_X, 8_0_5 for example)
//

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

// for trigger matching:
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

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
#include "TH1I.h"
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

#include "UserCode/ttbar-leptons-80X/interface/recordFuncs.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingBJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingMuons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingElectrons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingTaus.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingDRCleaning.h"

using namespace std;

namespace utils
	{
	namespace cmssw
		{
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

		/*
		FactorizedJetCorrector* getJetCorrector(TString baseDir, TString pf, bool isMC)
			{
			gSystem->ExpandPathName(baseDir);
			//TString pf(isMC ? "MC" : "DATA");
			// TString pf("Fall15_25nsV2_");
			pf += (isMC ? "MC" : "DATA");
			
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
		*/

		std::vector<double> smearJER(double pt, double eta, double genPt)
			{
			std::vector<double> toReturn(3,pt);
			if(genPt<=0) return toReturn;

			// These are from MacroUtils
			// FIXME: These are the 8 TeV values.
			//
			// eta=fabs(eta);
			// double ptSF(1.0), ptSF_err(0.06);
			// if(eta<0.8)                  { ptSF=1.061; ptSF_err=sqrt(pow(0.012,2)+pow(0.023,2)); }
			// else if(eta>=0.8 && eta<1.3) { ptSF=1.088; ptSF_err=sqrt(pow(0.012,2)+pow(0.029,2)); }
			// else if(eta>=1.3 && eta<1.9) { ptSF=1.106; ptSF_err=sqrt(pow(0.017,2)+pow(0.030,2)); }
			// else if(eta>=1.9 && eta<2.5) { ptSF=1.126; ptSF_err=sqrt(pow(0.035,2)+pow(0.094,2)); }
			// else if(eta>=2.5 && eta<3.0) { ptSF=1.343; ptSF_err=sqrt(pow(0.127,2)+pow(0.123,2)); }
			// else if(eta>=3.0 && eta<3.2) { ptSF=1.303; ptSF_err=sqrt(pow(0.127,2)+pow(1.303,2)); }
			// else if(eta>=3.2 && eta<5.0) { ptSF=1.320; ptSF_err=sqrt(pow(0.127,2)+pow(1.320,2)); }

			// TODO: 13TeV table for CMSSW_76X
			// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
			// following https://github.com/pfs/TopLJets2015/blob/master/TopAnalysis/src/CommonTools.cc
			eta=fabs(eta);
			double ptSF(1.0), ptSF_err(0.06);
			if      (eta<0.5) { ptSF=1.095; ptSF_err=0.018; }
			else if (eta<0.8) { ptSF=1.120; ptSF_err=0.028; }
			else if (eta<1.1) { ptSF=1.097; ptSF_err=0.017; }
			else if (eta<1.3) { ptSF=1.103; ptSF_err=0.033; }
			else if (eta<1.7) { ptSF=1.118; ptSF_err=0.014; }
			else if (eta<1.9) { ptSF=1.100; ptSF_err=0.033; }
			else if (eta<2.1) { ptSF=1.162; ptSF_err=0.044; }
			else if (eta<2.3) { ptSF=1.160; ptSF_err=0.048; }
			else if (eta<2.5) { ptSF=1.161; ptSF_err=0.060; }
			else if (eta<2.8) { ptSF=1.209; ptSF_err=0.059; }
			else if (eta<3.0) { ptSF=1.564; ptSF_err=0.321; }
			else if (eta<3.2) { ptSF=1.384; ptSF_err=0.033; }
			else if (eta<5.0) { ptSF=1.216; ptSF_err=0.050; }

			toReturn[0]=TMath::Max(0., (genPt+ptSF*(pt-genPt))/pt );
			toReturn[1]=TMath::Max(0., (genPt+(ptSF+ptSF_err)*(pt-genPt))/pt );
			toReturn[2]=TMath::Max(0., (genPt+(ptSF-ptSF_err)*(pt-genPt))/pt );
			return toReturn;
			}


		/* Using new stuff above
		std::vector<double> smearJER(double pt, double eta, double genPt)
			{
			std::vector<double> toReturn(3,pt);
			if(genPt<=0) return toReturn;
			
			// FIXME: These are the 8 TeV values.
			//
			eta=fabs(eta);
			double ptSF(1.0), ptSF_err(0.06);
			if(eta<0.5)                  { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
			else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
			else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
			else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
			else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
			
			toReturn[0]=TMath::Max(0.,(genPt+ptSF*(pt-genPt)));
			toReturn[1]=TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt)));
			toReturn[2]=TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt)));
			return toReturn;
			}

		//
		//std::vector<double> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
		std::vector<float> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
			{
			jecUnc->setJetEta(eta);
			jecUnc->setJetPt(pt);
			double relShift=fabs(jecUnc->getUncertainty(true));
			//std::vector<double> toRet;
			std::vector<float> toRet;
			toRet.push_back((1.0+relShift)*pt);
			toRet.push_back((1.0-relShift)*pt);
			return toRet;
			}
		*/

		/* using the one in src/MacroUtils.cc
		void updateJEC(pat::JetCollection &jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC)
			{
			for(size_t ijet=0; ijet<jets.size(); ijet++)
				{
				pat::Jet jet = jets[ijet];
				//correct JES
				LorentzVector rawJet = jet.correctedP4("Uncorrected");
				//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
				//LorentzVector rawJet(jet*toRawSF);
				jesCor->setJetEta(rawJet.eta());
				jesCor->setJetPt(rawJet.pt());
				jesCor->setJetA(jet.jetArea());
				jesCor->setRho(rho);
				jesCor->setNPV(nvtx);
				double newJECSF=jesCor->getCorrection();
				rawJet *= newJECSF;
				//jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
				jet.setP4(rawJet);

				//smear JER
				double newJERSF(1.0);
				if(isMC)
					{
					const reco::GenJet* genJet=jet.genJet();
					double genjetpt( genJet ? genJet->pt(): 0.);
					std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
					newJERSF=smearJER[0]/jet.pt();
					rawJet *= newJERSF;
					//jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
					jet.setP4(rawJet);
					// FIXME: change the way this is stored (to not storing it)
					// //set the JER up/down alternatives 
					// jets[ijet].setVal("jerup",   smearJER[1] );
					// jets[ijet].setVal("jerdown", smearJER[2] );
					}

				// FIXME: change the way this is stored (to not storing it)
				////set the JES up/down pT alternatives
				//std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
				//jets[ijet].setVal("jesup",    ptUnc[0] );
				//jets[ijet].setVal("jesdown",  ptUnc[1] );

				// FIXME: this is not to be re-set. Check that this is a desired non-feature.
				// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine. 
				//to get the raw jet again
				//jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));  
				}
			}
		*/

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

//string mc_decay("");
// Some MC datasets are inclusive, but analysis needs info on separate channels from them
// thus the events are traversed and the channel is found
// currently it is done only for TTbar channel (isTTbarMC)
//
// the sub-channel of MC is paired together with the distr_name string into the key of control distrs
// it is then printed out with dtag


std::map<std::pair <string,string>, TH1D> th1d_distr_control;
std::map<string, TH1D> th1d_distr_control_headers;

std::map<std::pair <string,string>, TH2D> th2d_distr_control;
std::map<string, TH2D> th2d_distr_control_headers;

std::map<std::pair <string,string>, TH3F> th3f_distr_control;
std::map<string, TH3F> th3f_distr_control_headers;

std::map<std::pair <string,string>, TH1I> th1i_distr_control;
std::map<string, TH1I> th1i_distr_control_headers;

std::map<std::pair <string,string>, double> weight_flow_control;


/*
std::map<string, std::map<string, TH1D>> th1d_distr_maps_control;
std::map<string, TH1D> th1d_distr_maps_control_headers;

std::map<string, std::map<string, TH1I>> th1i_distr_maps_control;
std::map<string, TH1I> th1i_distr_maps_control_headers;

std::map<string, std::map<string, TH2D>> th2d_distr_maps_control;
std::map<string, TH2D> th2d_distr_maps_control_headers;

std::map<string, std::map<string, TH3D>> th3d_distr_maps_control;
std::map<string, TH3D> th3d_distr_maps_control_headers;
*/







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

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		// th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;ID", 600, -300., 300.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		th1d_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;ID", 600, -300., 300.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	// th1i_distr_control[key].Fill(value, weight);
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("p_id")) == th1d_distr_control_headers.end() )
		{
		// th1d_distr_control_headers[string("p_id")] = TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.);
		th1d_distr_control_headers.insert( std::make_pair(string("p_id"), TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.)));
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



int fill_pt_pt(string control_point_name, double pt1, double pt2, double weight)
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
		th2d_distr_control.insert( std::make_pair(key, TH2D((mc_decay + control_point_name).c_str(), ";;Pt/E(GeV)", 400, 0., 400., 400, 0., 400.)));
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th2d_distr_control[key].Fill(pt1, pt2, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th2d_distr_control[control_point_name].Integral() << endl;

	if (th2d_distr_control_headers.find(string("pt_pt")) == th2d_distr_control_headers.end() )
		{
		// th2d_distr_control_headers[string("pt_e")] = TH2D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th2d_distr_control_headers.insert( std::make_pair(string("pt_e"), TH2D("Header of Pt vs Pt (E vs E) distributions", ";;Pt/E(GeV)", 400, 0., 400., 400, 0., 400.)));
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

	for(std::map<string, TH1D>::iterator it = th1d_distr_control_headers.begin(); it != th1d_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1D * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}

	for(std::map<string, TH1I>::iterator it = th1i_distr_control_headers.begin(); it != th1i_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1I * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}

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

bool conf_record_electrons =  runProcess.getParameter<bool> ("conf_record_electrons");
bool conf_record_muons     =  runProcess.getParameter<bool> ("conf_record_muons");
bool conf_record_taus_ID   =  runProcess.getParameter<bool> ("conf_record_taus_ID");
bool conf_record_taus_kino =  runProcess.getParameter<bool> ("conf_record_taus_kino");

cout << "recording stuff:" << endl;
cout << "(el-s) " << conf_record_electrons << '\t';
cout << "(mu-s) " << conf_record_muons << '\t';
cout << "(tau-ID-s) " << conf_record_taus_ID << '\t';
cout << "(tau-kino-s) " << conf_record_taus_kino << '\t' << endl;

string jetHLT = runProcess.getParameter<std::string>("jetHLT"),
	muHLT_MC1   = runProcess.getParameter<std::string>("muHLT_MC1"), muHLT_MC2   = runProcess.getParameter<std::string>("muHLT_MC2"),
	muHLT_Data1 = runProcess.getParameter<std::string>("muHLT_Data1"), muHLT_Data2 = runProcess.getParameter<std::string>("muHLT_Data2");

cout << "Triggers:" << endl;
cout << muHLT_MC1 << '\t' << muHLT_MC2 << '\t' << muHLT_Data1 << '\t' << muHLT_Data2 << endl;
cout << jetHLT << endl;

// Kino cuts
double jettaufr_jet_kino_cuts_pt          = runProcess.getParameter<double>("jettaufr_jet_kino_cuts_pt");
double jettaufr_jet_kino_cuts_eta         = runProcess.getParameter<double>("jettaufr_jet_kino_cuts_eta");
double jettaufr_tau_kino_cuts_pt          = runProcess.getParameter<double>("jettaufr_tau_kino_cuts_pt");
double jettaufr_tau_kino_cuts_eta         = runProcess.getParameter<double>("jettaufr_tau_kino_cuts_eta");

cout << "Kino cuts" << endl;
cout << "jets: (pt)\t" << jettaufr_jet_kino_cuts_pt << "\t(eta)" << jettaufr_jet_kino_cuts_eta << endl;
cout << "taus: (pt)\t" << jettaufr_tau_kino_cuts_pt << "\t(eta)" << jettaufr_tau_kino_cuts_eta << endl;

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
TString outUrl = runProcess.getParameter<std::string>("outfile");
TString outdir = runProcess.getParameter<std::string>("outdir");

cout << "Output directory: " << outdir << endl;

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

bool isNoHLT = dtag.Contains("noHLT");

bool isJetHT = (!isMC) && dtag.Contains("JetHT");
bool isSingleMuon = (!isMC) && dtag.Contains("SingleMuon");

bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );

// Summer16_23Sep2016BCDV4_DATA_        Summer16_23Sep2016EFV4_DATA_        Summer16_23Sep2016GV4_DATA_        Summer16_23Sep2016HV4_DATA_
bool period_BCD = !isMC && (dtag.Contains("2016B") || dtag.Contains("2016C") || dtag.Contains("2016D"));
bool period_EF  = !isMC && (dtag.Contains("2016E") || dtag.Contains("2016F"));
bool period_G   = !isMC && (dtag.Contains("2016G"));
bool period_H   = !isMC && (dtag.Contains("2016H"));

TString outTxtUrl = outUrl + ".txt";
FILE *outTxtFile = NULL;
if (!isMC) outTxtFile = fopen (outTxtUrl.Data (), "w");
printf ("TextFile URL = %s\n", outTxtUrl.Data ());

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


//MC normalization (to 1/pb)
if(debug) cout << "DEBUG: xsec: " << xsec << endl;

// ------------------------------------- jet energy scale and uncertainties 
TString jecDir = runProcess.getParameter < std::string > ("jecDir");
gSystem->ExpandPathName (jecDir);

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


// jet parameters: pt, eta, radius
// in the datasets jet eta varies within [-2.5, 2.5]
//TH3F* jets_distr = (TH3F*) new TH3F("jets_distr", ";;", 500, 0, 500, 400, -3, 3, 100, 0, 2);
//TH3F* tau_jets_distr = (TH3F*) new TH3F("tau_jets_distr", ";;", 500, 0, 500, 400, -3, 3, 100, 0, 2);

//// good bins 1
//Float_t bins_pt[11] = { 0, 29, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
//Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
//Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
//	0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges

//TH3F* wjets_jets_distr      = (TH3F*) new TH3F("wjets_jets_distr",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);
//TH3F* wjets_tau_jets_distr  = (TH3F*) new TH3F("wjets_tau_jets_distr",  ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);
//TH3F* wjets_distinct_tau_distr = (TH3F*) new TH3F("wjets_distinct_tau_distr",  ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);

// TH3F* qcd_jets_distr      = (TH3F*) new TH3F("qcd_jets_distr",        ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);
// TH3F* qcd_tau_jets_distr  = (TH3F*) new TH3F("qcd_tau_jets_distr",    ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);
// TH3F* qcd_distinct_tau_distr = (TH3F*) new TH3F("qcd_distinct_tau_distr",  ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);

TH1D* wjets_taujet_distance = (TH1D*) new TH1D("wjets_taujet_distance", ";Distance [phi-eta];Events",            100, 0.,  2.  );
TH1I* wjets_jet_origin      = (TH1I*) new TH1I("wjets_jet_origin",    ";PDg ID;Particles",            100, 0,  100  );
TH1I* wjets_taujet_origin   = (TH1I*) new TH1I("wjets_taujet_origin", ";PDg ID;Particles",            100, 0,  100  );
TH1D* qcd_taujet_distance = (TH1D*) new TH1D("qcd_taujet_distance",   ";Distance [phi-eta];Events",            100, 0.,  2.  );
TH1I* qcd_jet_origin      = (TH1I*) new TH1I("qcd_jet_origin",    ";PDG ID;N Particles",            100, 0,  100  );
TH1I* qcd_taujet_origin   = (TH1I*) new TH1I("qcd_taujet_origin", ";PDG ID;N Particles",            100, 0,  100  );

/* For reference, some PDG IDs:
 * QUARKS
 * d  1
 * u  2
 * s  3
 * c  4
 * b  5
 * t  6
 * b' 7
 * t' 8
 * g 21
 * gamma 22
 * Z     23
 * W     24
 * h     25
 * e, ve     11, 12
 * mu, vmu   13, 14
 * tau, vtau 15, 16
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


for(size_t f=0; f<urls.size();++f)
	{
	//fprintf(csv_out, "Processing file: %s\n", urls[f].c_str());
	cout << "Processing file: " << urls[f].c_str() << "\n";
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());


	int treeStep (ev.size()/50);
	unsigned int iev = 0;

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
					//
					// in genParticles tree such a decay looks like:
					// 36: 15 2        40.6392,1.48436,2.06287,1.777
					//  <-  15 23;
					//          |-> 16 1; -211 1; 
					// -- no intermediate between tau, neutrinos and visible part
					LorentzVector vis_ds(0,0,0,0);
					for (int j = 0; j < n_daughters; ++j)
						{
						const reco::Candidate * d = p.daughter(j);
						unsigned int d_id = fabs(d->pdgId());
						if (d_id == 12 || d_id == 14 || d_id == 16) continue;
						if (d->status() == 1) vis_ds += d->p4();
						// it should be useless -- all subdaughters are photons from neutral pions
						// and they are visible..
						// but still:
						for (int u = 0; u < d->numberOfDaughters(); ++u)
							{
							const reco::Candidate * d2 = d->daughter(u);
							unsigned int d_id = fabs(d2->pdgId());
							if (d_id == 12 || d_id == 14 || d_id == 16) continue;
							vis_ds += d2->p4();
							}
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

		/* no need for ttbar channels in jet-tau fake rate study
		// extract the TTbar decay in this file
		// by traversing separate t quarks
		//if (debug) {
		if (isTTbarMC && genHandle.isValid()) {
			//string mc_decay(""); // move to all job parameters
			// every found t-quark decay branch will be added as a substring to this string (el, mu, tau, q)
			// TODO: what to do if there are > t-decays than 1 ttbar pair in an event? (naive traverse approach also doesn't work?)
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
						mc_decay += string("bar");
					}
				}
			}
			// so, mc_decay will be populated with strings matching t decays
			// hopefully, in TTbar sample only ttbar decays are present and it is not ambigous
		}

		if (debug) {
			cout << "MC suffix " << mc_decay << " is found\n";
			}

		if (!mc_decay.empty()) mc_decay = string("_") + mc_decay;
		*/

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
		/* jobs still fail here?
		if (debug && (iev % treeStep == 0))
			{
			printf (".");
			//if(!debug) fflush (stdout); // Otherwise debug messages are flushed
			}
		*/

		edm::EventBase const & myEvent = ev;

		// ---------------------------------- MC shaping to data
		// MC is weighted according to distributions of a bunch of data properties

		// NLO -1 corrections
		double weight_Gen(1.);

		// there are also the (dissabled now, since NLO samples are used) HT-binned and pthat-binned stitching of LO and NLO

		// ---------------------------------- pileup weight
		double weight_PU         (1.0);
		double weight_PU_up      (1.0);
		double weight_PU_down    (1.0);

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


		// ---------------------------------- RHO variable for PU control
		// rho = sum of energy flow per dR angle
		// there are 4 of them in 80X MINIAODs:
		// fixedGridRhoFastjetAll
		// fixedGridRhoFastjetCentral
		// fixedGridRhoFastjetCentralNeutral
		// fixedGridRhoFastjetCentralChargedPileUp
		// people in PU Jet ID example use fixedGridRhoFastjetAll
		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;

		//
		// DERIVE WEIGHTS TO APPLY TO SAMPLE
		//

		// weight init 1
		fill_1d(string("weightflow"), 300, 0, 300,   1, weight);

		// ---------------------------------- these are weird NLO -1 weights
		// TODO: figure out how exactly they correct for NLO
		// Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
		if(isNLOMC)
			{

			fwlite::Handle<GenEventInfoProduct> evt;
			evt.getByLabel(ev, "generator");
			if(evt.isValid())
				{
				weight_Gen = (evt->weight() > 0 ) ? 1. : -1. ;
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

		weight *= weight_Gen;
		weight_up *= weight_Gen;
		weight_down *= weight_Gen;
		rawWeight *=weight_Gen;

		// weight Gen 2
		fill_1d(string("weightflow"), 300, 0, 300,   2, weight);

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
		for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
			{
			if(utils::isGoodVertex(vtx[ivtx]))
				{
				if(nGoodPV==0) goodPV=vtx[ivtx];
				nGoodPV++;
				}
			}

		// ----------------------------------------- Apply PILEUP reweighting
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

		weight *= weight_PU;
		weight_up *= weight_PU_up;
		weight_down *= weight_PU_down;

		// weight PU 3
		fill_1d(string("weightflow"), 300, 0, 300,   3, weight);

		// pu distrs
		fill_1d( string("pileup_beforetrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
		fill_1d( string("pileup_beforetrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);
		// nGoodPV = vtx.size() now

		// pu distrs
		fill_1d( string("pileup_beforetrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_beforetrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
		// nGoodPV = vtx.size() now

		// RHO distrs
		fill_1d( string("rho_beforetrig_rawWeight"), 100, 0, 100, rho, rawWeight);
		fill_1d( string("rho_beforetrig_weight"),    100, 0, 100, rho, weight);


		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		//num_inters = 1;
		if (num_inters>99) num_inters = 99;
		if (nGoodPV>100) nGoodPV = 99;


		// -------------------------------------   END of SECTION1, MC weights

		fill_1d(string("weightflow"), 300, 0, 300,   10, weight);

		// -------------------------------------   Basic event selection

		// ---------------------- Orthogonalize Run2015B PromptReco+17Jul15 mix
		// let's remove Run2015B
		// if(isRun2015B)
		// {
		// if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
		// }

		// it's not needed with the latest versions of RunB rereconstruction
		
		// -------------------------------------------------- Skip bad LUMI
		// 80X, v2
		if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock())) continue; 
		// Notice: it is the first continue in the event loop
		// there is no sum_weights_pass_lumi -- lumi is for data only..

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		//increment( string("weightflow_weight_passed_lumi"), weight ); // should not matter
		//increment( string("weightflow_weight_up_passed_lumi"), weight_up ); // should not matter
		//increment( string("weightflow_weight_down_passed_lumi"), weight_down ); // should not matter

		// passlumi 5
		fill_1d(string("weightflow"), 300, 0, 300,   11, weight);
		/*
		fill_1d(string("weightflow_mu"), 300, 0, 300,   5, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   5, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 5, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 5, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 5, weight);
		*/

		// --------------------------------------------- HLT TRIGGER
		// ---------------- and require compatibilitiy of the event with the PD

		// it seems (from ConfDB /dev/CMSSW/CMSSW_8_0_25/...) PFJet140 should not be prescaled
		// and PFJet450 is also not prescaled
		/* passing from cfg.py
		string jetHLT("HLT_PFJet140_v"), //jetHLT("HLT_PFJet40_v"), // jetHLT("HLT_AK4PFJet30_v*"),
			muHLT_MC1("HLT_IsoMu24_v4"), muHLT_MC2("HLT_IsoTkMu24_v4"),
			muHLT_Data1("HLT_IsoMu24_v*"), muHLT_Data2("HLT_IsoTkMu24_v*");
		*/

		// Spring16 MC has these couple issues with HLT
		// some sets don't have it at all -- noHLT (madgraf QCD, DY)
		// some have it, but under HLT2 process -- reHLT (aMCatNLO WJets, powheg TTbar)
		edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT2");

		// in QCD selection the trigger jet is removed:
		// find matching trig objects (dR < 0.3)
		// if there is only 1 matching trig object with the name ~ HLT_PFJet40*
		// skip the corresponding jet
		// otherwise use all jets
		// TODO: check this scheme, if it doesn't remove the 40GeV bump -- try removing closest matching jet
		if (debug) cout << "Trigger objects from MiniAOD workbook" << endl;
		//#include <InputTag.h>
		fwlite::Handle< vector<pat::TriggerObjectStandAlone> > triggerObjectsHandle;
		//myEvent.getByLabel( edm::InputTag("selectedPatTrigger"), triggerObjects );
		triggerObjectsHandle.getByLabel( ev, "selectedPatTrigger" );

		if (!triggerObjectsHandle.isValid())
			{
			cout << "!triggerObjectsHandle.isValid()" << endl;
			}
		vector<pat::TriggerObjectStandAlone> trig_objs = *triggerObjectsHandle;

		// get TriggerNames (needed to "unpack" the trigger path names of TriggerObject
		edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
		edm::InputTag * trigResultsTag;
		//edm::TriggerResults trigger_results = ev.triggerResults ("HLT");

		if (!tr.isValid ()){
			//cout << "Trigger HLT2 is not valid, trying HLT" << endl;
			tr = ev.triggerResultsByName ("HLT");
			trigResultsTag = new edm::InputTag("TriggerResults","","HLT"); //make sure have correct process on MC
			if (!tr.isValid ()){
				cout << "Trigger HLT is not valid, exiting" << endl;
				return false;
				}
			}
		else
			{
			trigResultsTag = new edm::InputTag("TriggerResults","","HLT2"); //make sure have correct process on MC
			}

		ev.getByLabel(*trigResultsTag, trigResults);
		const edm::TriggerNames& trigNames = ev.triggerNames(*trigResults);   

		if(debug){
			cout << "Printing HLT trigger list" << endl;
			//cout << "-- Commented out --" << endl;
			//for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames) cout << *trnames << endl;
			cout << "----------- End of trigger list ----------" << endl;
			//return 0;
			}

		// Need either to simulate the HLT (https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#How_to_easily_emulate_HLT_paths) to match triggers.
		// trigger for Jet-heavy events in analysis note AN-2012-489: HLT_Jet30_v1

		// TODO: find the corresponding trigger in 2016

		bool JetHLTTrigger = false;
		bool MuonHLTTrigger = false;

		if (isNoHLT)
			{
			JetHLTTrigger = true;
			MuonHLTTrigger = true;
			}
		else if (tr.isValid())
			{
			// bool hltTrigger( utils::passTriggerPatterns(tr, "HLT_Jet30_v*") );
			//bool hltTrigger( utils::passTriggerPatterns(tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*") );
			// bool hltTrigger( utils::passTriggerPatterns(tr, "HLT_PFJet40_v*", "HLT_IsoMu22*", "HLT_IsoTkMu22*") );
			JetHLTTrigger = utils::passTriggerPatterns(tr, jetHLT + "*"); // NOTICE the pattern-matching wild card *
			MuonHLTTrigger = (isMC ?
				//utils::passTriggerPatterns(tr, "HLT_IsoMu22*", "HLT_IsoTkMu22*");
				utils::passTriggerPatterns (tr, muHLT_MC1, muHLT_MC2) : // tecommended inn ttbar trig for reHLT
				//utils::passTriggerPatterns (tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*") :
				utils::passTriggerPatterns (tr, muHLT_Data1, muHLT_Data2));
			}

		// if (!(JetHLTTrigger || MuonHLTTrigger)) continue; // No orthogonalization -- run on only 1 trigger type of datasets

		/*
		 * making HLTs not exclusive:
		 * events from SingleMuon dataset passing Muon HLTs and selections are needed
		 * .....  from JetHT dataset passing Jet HLT and selections
		 * --- if it is JetHT check if they pass HLT_Jet
		 *              SingleMuon --- HLT_Muon
		if (JetHLTTrigger && MuonHLTTrigger)
			{
			hlt_channel = string("HLTjetmu_");
			}
		else if (JetHLTTrigger)
			{
			hlt_channel = string("HLTjet_");
			}
		else if (MuonHLTTrigger)
			{
			hlt_channel = string("HLTmu_");
			}
		else continue;
		*/

		if (!(JetHLTTrigger || MuonHLTTrigger)) continue;

		vector<string> hlt_channels;
		if (JetHLTTrigger)  hlt_channels.push_back("HLTjet_");
		if (MuonHLTTrigger) hlt_channels.push_back("HLTmu_");

		// pu distrs
		fill_1d( string("pileup_passtrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
		fill_1d( string("pileup_passtrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);
		// nGoodPV = vtx.size() now

		// pu distrs
		fill_1d( string("pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);
		// nGoodPV = vtx.size() now

		// RHO distrs
		fill_1d( string("rho_passtrig_rawWeight"), 100, 0, 100, rho, rawWeight);
		fill_1d( string("rho_passtrig_weight"),    100, 0, 100, rho, weight);

		// TODO: ----------------------------- HLT efficiency scale factors
		// one should run it on the fired trigger objects,
		// I run it on selection candidates now
		// which is done below, when the candidates are selected

		// double HLT_efficiency_sf = 1.0;

		//HLT_efficiency_sf *= eTrigger  ? eHLT_sf[] : 1 ;
		//HLT_efficiency_sf *= muTrigger ? muHLT_SF[] : 1 ;




		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		// passtrig 12
		fill_1d(string("weightflow"), 300, 0, 300,   12, weight);

		//increment( string("weightflow_weight_passed_trig"), weight ); // should not matter
		//increment( string("weightflow_weight_up_passed_trig"), weight_up ); // should not matter
		//increment( string("weightflow_weight_down_passed_trig"), weight_down ); // should not matter

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

		// ------------------------------------------------- Apply MET FILTERS
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
				cout << "Printing PAT/RECO trigger list for MET filters here" << endl;
				for(edm::TriggerNames::Strings::const_iterator trnames = metFilters.triggerNames().begin(); trnames!=metFilters.triggerNames().end(); ++trnames)
					cout << *trnames << endl;
				cout << "----------- End of trigger list ----------" << endl;
				//return 0;
			}
			//std::vector<std::string>& patterns("Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_CSCTightHalo2015Filter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter");
			if (!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*", "Flag_CSCTightHalo2015Filter*", "Flag_EcalDeadCellTriggerPrimitiveFilter*"))
				continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_goodVertices")) continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter")) continue;
		}
		

		// passmetfilters 13
		fill_1d(string("weightflow"), 300, 0, 300,   13, weight);

		if(debug)
			{
			cout << "met filters applied here" << endl;
			}

		// ------------------------- END of SECTION 2, event cuts

		fill_1d(string("weightflow"), 300, 0, 300,   20, weight);

		// ------------------------- event physics and the corresponding selection

		//------------------------- load all the objects we will need to access

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

		pat::PhotonCollection photons;
		fwlite::Handle<pat::PhotonCollection> photonsHandle;
		photonsHandle.getByLabel(ev, "slimmedPhotons");
		if(photonsHandle.isValid() ) photons = *photonsHandle;

		pat::METCollection mets;
		fwlite::Handle<pat::METCollection> metsHandle;
		metsHandle.getByLabel(ev, "slimmedMETs");
		if(metsHandle.isValid() ) mets = *metsHandle;
		pat::MET met = mets[0];
		// LorentzVector met = mets[0].p4 ();


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

		// triggered objects:
		//const TriggerObjectStandAloneCollection pat::PATObject< ObjectType >::triggerObjectMatchesByPath
		// returns a collection of trigger objects for a pat::Jet?
		if (debug && JetHLTTrigger)
			{
			cout << "JetHLTTrigger matched" << endl;
			for (int i=0; i<jets.size(); i++)
				{
				//const pat::TriggerObjectStandAloneCollection trig_objects = jets[i].triggerObjectMatchesByPath("HLT_PFJet40_v*");
				const pat::TriggerObjectStandAloneCollection trig_objects = jets[i].triggerObjectMatches();
				cout << "trig objects matching jet " << i;
				for (int u=0; u<trig_objects.size(); u++)
					cout << " " << trig_objects[u].pdgId();
				cout << endl;
				}
			}
		// this doesn't return any objects for any jet in JetHT matched events

		// triggered objects for muons:
		//const TriggerObjectStandAloneCollection pat::PATObject< ObjectType >::triggerObjectMatchesByPath
		// returns a collection of trigger objects for a pat::Jet?
		if (debug && MuonHLTTrigger)
			{
			cout << "MuonHLTTrigger matched" << endl;
			for (int i=0; i<muons.size(); i++)
				{
				//const pat::TriggerObjectStandAloneCollection trig_objects = muons[i].triggerObjectMatchesByPath("HLT_IsoMu24_v*");
				const pat::TriggerObjectStandAloneCollection trig_objects = muons[i].triggerObjectMatches();
				cout << "trig objects matching muon " << i;
				for (int u=0; u<trig_objects.size(); u++)
					cout << " " << trig_objects[u].pdgId();
				cout << endl;
				}
			}
		// this doesn't return any objects for any jet in JetHT matched events



		/*
		// in pat exercises they use this:
		fwlite::Handle< pat::TriggerEvent > triggerEventHandle;
		triggerEventHandle.getByLabel(ev, "patTriggerEvent" );
		pat::TriggerEvent triggerEvent = *triggerEventHandle;
		if (debug) cout << "TriggerEvent beamMode = " << triggerEvent.beamMode() << endl;
		*/

		//edm::Handle< pat::TriggerEvent > triggerEventHandle;
		//ev.getByLabel( "patTriggerEvent", triggerEvent );
		//iEvent.getByLabel( "patTriggerEvent", triggerEvent );
		// and go manually over the matching procedure
		/*
		pat::TriggerAlgorithmRefVector::const_iterator itrBit = algoBits.begin();   
		int nbit = 0;
		while ( itrBit != algoBits.end())
			{
			// put here the logic of your analysis (compare and check decision)
			// consult the documentation for the TriggerAlgorithm class        
			itrBit++;    
			nbit++;
			}
		if (debug)
			{
			cout << "bla bits " << nbit << endl;
			}
		*/

		// CONTROLINFO
		// Control values for raw particles:


		//
		//
		// BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS. Whatever that means
		//
		//
		


		if(debug){
			cout << "got objects from the event, starting the analysis" << endl;
			}

		//
		// LEPTON ANALYSIS
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
		pat::ElectronCollection selElectrons;
		unsigned int nVetoE(0);

		processElectrons_ID_ISO_Kinematics(electrons, goodPV, rho, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
			35., 2.4, 15., 2.5, selElectrons, elDiff, nVetoE, conf_record_electrons, debug);

		if(debug){
			cout << "processed electrons" << endl;
			}

		// ---------------------------------- MUONS SELECTION
		LorentzVector muDiff(0., 0., 0., 0.);
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
			30., 2.4, 10., 2.5, selMuons, muDiff, nVetoMu, conf_record_muons, debug);

		if(debug){
			cout << "processed muons" << endl;
			}





		// Finally, merge leptons for cross-cleaning with taus and jets:

		std::vector<patUtils::GenericLepton> selLeptons;
		for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
		for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
		std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);


		if(debug){
			cout << "merged selected leptons" << endl;
			}




		// ------------------------------------------ TAUS SELECTION

		//int processTaus_ID_ISO(pat::TauCollection& taus, double weight, // input
		//	string& tauID_decayMode, string& tauID,               // config/cuts
		//	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
		//	pat::TauCollection& selTaus,                          // output
		//	bool record, bool debug) // more output

		pat::TauCollection IDtaus, selTaus;

		processTaus_ID_ISO(taus, weight, tau_decayMode, tau_ID, tau_againstMuon, tau_againstElectron, IDtaus, conf_record_taus_ID, debug); 

		//int processTaus_Kinematics(pat::TauCollection& taus,          // input
		//	double weight,
		//	double pt_cut, double eta_cut,
		//	pat::TauCollection& selTaus,                          // output
		//	bool record, bool debug) // more output

		processTaus_Kinematics(IDtaus, weight, jettaufr_tau_kino_cuts_pt, jettaufr_tau_kino_cuts_eta, selTaus, conf_record_taus_kino, debug);

		if(debug){
			cout << "selected taus [individual]" << endl;
			}



		// ------------------------------------------ select the taus cleaned from leptons


		//int crossClean_in_dR(pat::TauCollection& selTaus, std::vector<patUtils::GenericLepton>& leptons,
		//	float min_dR,
		//	pat::TauCollection& selTausNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		pat::TauCollection selTausNoMuons, selTausNoElectrons, selTausNoLep;
		crossClean_in_dR(selTaus, selElectrons, 0.4, selTausNoElectrons, weight, string("selTausNoElectrons"), true, debug);
		crossClean_in_dR(selTaus, selMuons,     0.4, selTausNoMuons,     weight, string("selTausNoMuons"),     true, debug);
		crossClean_in_dR(selTaus, selLeptons,   0.4, selTausNoLep,       weight, string("selTausNoLep"),       true, debug);

		// actually these will be recorded in 2d pt-eta distr-s of crossClean-ing
		// check the cleaning of taus from leptons per eta (the problem with high fake rate in endcups):
		for (size_t i = 0; i < selTaus.size(); ++i)
			fill_1d(string("control_selTaus_eta"), 128, -3.0, 3.0, selTaus[i].eta(), weight);
		for (size_t i = 0; i < selTausNoElectrons.size(); ++i)
			fill_1d(string("control_selTausNoElectrons_eta"), 128, -3.0, 3.0, selTausNoElectrons[i].eta(), weight);
		for (size_t i = 0; i < selTausNoMuons.size(); ++i)
			fill_1d(string("control_selTausNoMuons_eta"),     128, -3.0, 3.0, selTausNoMuons[i].eta(), weight);
		for (size_t i = 0; i < selTausNoLep.size(); ++i)
			fill_1d(string("control_selTausNoLep_eta"),       128, -3.0, 3.0, selTausNoLep[i].eta(), weight);

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
		fill_1d(string("weight_tauIDsf"), 200, 0., 2.,   weight_tauIDsf, 1);;

		weight_without_tauIDsf = weight;
		// weight *= weight_tauIDsf;// apply tau weight at selection

		if(debug){
			cout << "processed taus" << " N selTausNoLep = " << selTausNoLep.size() << endl;
			cout << "ID SF = " << weight_tauIDsf << endl;
			}




		//
		// ----------------------------------------------- JET/MET ANALYSIS
		//

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

		//int processJets_CorrectJES_SmearJERnJES_ID_ISO_Kinematics(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
		//	bool isMC, double weight,
		//	double rho, unsigned int nGoodPV,
		//	FactorizedJetCorrector *jesCor,
		//	JetCorrectionUncertainty *totalJESUnc,
		//	double dR_max, // for jet matching in jet corrections smearing for MC
		//	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
		//	string& jetID,
		//	string& jetPUID,
		//	double pt_cut, double eta_cut,
		//	TRandom *r3,   // the randomizer for the smearing
		//	LorentzVector& full_jet_corr, pat::JetCollection& selJets,                          // output
		//	bool record, bool debug) // more output

		LorentzVector full_jet_corr(0., 0., 0., 0.);
		pat::JetCollection IDjets;
		//string jetID("Loose");
		//string jetPUID("MediumPU");

		processJets_CorrectJES_SmearJERnJES_ID_ISO(jets, genJets, isMC, weight, rho, nGoodPV, jesCor, totalJESUnc, 0.4/2,
			jet_resolution_in_pt, jet_resolution_sf_per_eta, jet_m_systematic_variation, jetID, jetPUID, with_PU, r3, full_jet_corr, IDjets, true, debug);


		fill_3d(string("control_jet_full_jet_corr_pX_pY_pZ"), 10, -100., 100., 10, -100., 100., 10, -100., 100.,  full_jet_corr.X(), full_jet_corr.Y(), full_jet_corr.Z(), weight);
		// 1000 bins

		fill_2d(string("control_jet_full_jet_corr_pX_pY"), 100, 0., 100., 100, 0., 100.,  full_jet_corr.X(), full_jet_corr.Y(), weight);
		fill_1d(string("control_jet_full_jet_corr_pZ"),    100, 0., 100., full_jet_corr.Z(), weight);
		// 10 000 and 100 bins

		met.setP4(met.p4() - full_jet_corr); // just return the full correction and propagate in place

		fill_1d(string("control_met_slimmedMETs_fulljetcorrs_pt"), 200, 0., 200., met.pt(), weight);
		//fill_2d(string("control_met_slimmedMETs_fulljetcorrs_pt"), 200, 0., 200., 200, -4., 4., met.pt(), met.eta(), weight);

		//int processJets_Kinematics(pat::JetCollection& jets, // input
		//	//bool isMC, double weight,
		//	double pt_cut, double eta_cut,
		//	pat::JetCollection& selJets,                 // output
		//	bool record, bool debug) // more output

		pat::JetCollection selJets;
		processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jettaufr_jet_kino_cuts_pt, jettaufr_jet_kino_cuts_eta, selJets, true, debug);

		// ---------------------------- Clean jet collections from selected leptons
		// TODO: add gamma-cleaning as well?

		//int crossClean_in_dR(pat::JetCollection& selJets, std::vector<patUtils::GenericLepton>& selLeptons,
		//	float min_dR,
		//	pat::JetCollection& selJetsNoLep, // output
		//	string control_name,
		//	bool record, bool debug) // more output

		pat::JetCollection selJetsNoLep;
		crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), true, debug);



		// ---------------------------- Clean jet collection from selected taus

		//int crossClean_in_dR(pat::JetCollection& selJets, pat::TauCollection& selTaus,
		//	float min_dR,
		//	pat::JetCollection& selJetsOut, // output
		//	string control_name,
		//	bool record, bool debug) // more output
		pat::JetCollection selJetsNoLepNoTau;
		crossClean_in_dR(selJetsNoLep, selTausNoLep, 0.4, selJetsNoLepNoTau, weight, string("selJetsNoLepNoTau"), true, debug);

		if(debug){
			cout << "processed jets" << endl;
			}



		// check distances between jets:
		for (int i=0; i<selJetsNoLep.size(); i++)
			{
			for (int j=i+1; j<selJetsNoLep.size(); j++)
				{
				double delta_R = reco::deltaR(selJetsNoLep[i], selJetsNoLep[j]);
				fill_1d(string("control_selJetsNoLep_deltaR"), 60, 0., 6., delta_R, weight);
				fill_2d(string("control_selJetsNoLep_deltaR_eta1"), 60, 0., 6., 20, -3., 3., delta_R, selJetsNoLep[i].eta(), weight);
				}
			}

		// check distances between taus:
		for (int j=0; j<selTausNoLep.size(); j++)
			{
			for (int i=j+1; i<selTausNoLep.size(); i++)
				{
				double delta_R = reco::deltaR(selTausNoLep[i], selTausNoLep[j]);
				fill_1d(string("control_selTausNoLep_deltaR"), 60, 0., 6., delta_R, weight);
				fill_2d(string("control_selTausNoLep_deltaR_eta1"), 60, 0., 6., 20, -3., 3., delta_R, selTausNoLep[j].eta(), weight);
				}
			}

		// check distances between jets and taus:
		for (int j=0; j<selTausNoLep.size(); j++)
			{
			unsigned int number_of_matches = 0;
			for (int i=0; i<selJetsNoLep.size(); i++)
				{
				double delta_R = reco::deltaR(selJetsNoLep[i], selTausNoLep[j]);
				fill_1d(string("control_selJetsNoLep_selTausNoLep_deltaR"), 60, 0., 6., delta_R, weight);
				fill_2d(string("control_selJetsNoLep_selTausNoLep_deltaR_tau_eta"), 60, 0., 6., 20, -3., 3., delta_R, selTausNoLep[j].eta(), weight);
				if (delta_R < tau_fake_distance)
					number_of_matches += 1;
				}
			fill_1d(string("control_selTausNoLep_selJetsNoLep_number_of_deltaR_matches"), 10, 0., 10., number_of_matches, weight);
			fill_2d(string("control_selTausNoLep_selJetsNoLep_number_of_deltaR_matches_eta"), 10, 0., 10., 20, -3., 3., number_of_matches, selTausNoLep[j].eta(), weight);
			}

		/*
		// --------------------------- B-TAGGED JETS
		pat::JetCollection selBJets;

		//int processBJets_BTag(pat::JetCollection& jets, bool isMC, double weight, // input
		//	BTagCalibrationReader& btagCal, BTagSFUtil& btsfutil,
		//	string& b_tagger_label, float b_tag_WP,
		//	float beff, leff,
		//	pat::JetCollection& selBJets,                          // output
		//	bool record, bool debug) // more output

		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
		string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		float btag_WP = 0.8484; // medium
		processBJets_BTag(selJetsNoLepNoTau, isMC, weight, btagCal, btsfutil, btagger_label, btag_WP, 0.747, 0.13, selBJets, true, debug);

		if(debug){
			cout << "processed b-tagged jets" << endl;
			}
		*/





		// ------------------------------------------ the taus far from jets
		pat::TauCollection selTausNoLepNoJet;
		crossClean_in_dR(selTausNoLep, selJetsNoLep, 0.4, selTausNoLepNoJet, weight, string("selTausNoLepNoJet"), false, debug);




		// ------------------------- END of SECTION 3, selection of particles, SF weights
		fill_1d(string("weightflow"), 300, 0, 300,   30, weight);




		// -------------------------------------------------- All particles are selected
		// now the channel/selection can be done

		if(debug){
			cout << "all particle-objects are processed, checking channel selection" << endl;
			}

		unsigned int n_leptons = selLeptons.size();

		unsigned int n_taus = selTausNoLep.size();
		unsigned int n_jets = selJetsNoLep.size();
		//unsigned int n_bjets = selBJets.size();

		// Select jets for tau-fake-rate study, save the counts

		// W+jets selection:
		bool selection_W_jets = (n_leptons > 0) && (n_jets > 0);

		// And all-jets selection:
		bool selection_QCD_jets = (n_jets > 1);


		// select jets:
		//   - don't consider the jet if it is the only one matching HLT jet
		// compare with taus
		// save fake-rate

		//TODO: remove HLT jet from the probe jets, if there is only 1 HLT jet
		/*
		HLT_jet = ;
		int hlt_bias_jet = -1;
		unsigned int n_hlt_jets = 0;

		for (size_t ijet = 0; ijet < selJetsNoLepNoTau.size(); ++ijet)
			{
			pat::Jet& jet = selJetsNoLepNoTau[ijet];


			if (reco::deltaR(jet, HLT_jet) < 0.1)
				{
				// 
				if (n_hlt_jets == 0)
					{
					hlt_bias_jet = ijet;
					}
				else if (n_hlt_jets == 1)
					{
					hlt_bias_jet = -1;
					}
				n_hlt_jets += 1;
				}
			}

		if (hlt_bias_jet > 0)
			{
			probeJets = // selected_jets without the bias one
			}
		else
			{
			probeJets = // selected_jets
			}
		*/

		// Check among 2 channels of selection:
		//     QCD-channel -- 2 or more jets
		//     W+jets      -- only 1 iso muon, one or more jets

		bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;

		bool isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && clean_lep_conditions;
		//bool isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && clean_lep_conditions;

		bool Wjets_selection = clean_lep_conditions && isSingleMu && (selJetsNoLep.size() > 0);
		bool W1jet_selection = clean_lep_conditions && isSingleMu && (selJetsNoLep.size() == 1);
		//bool QCD_selection  = selJetsNoLep.size() > 1;
		// Trying get HT > 100 jets --- to match the smallest HT bin I have in QCD MC
		// HT = sum of jet pT-s
		// my jets are 20 GeV
		// ...
		// let's do 2 separate requirements and & them
		double event_ht = 0.;
		for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
			{
			event_ht += selJetsNoLep[ijet].pt();
			}
		bool QCD_selection  = (selJetsNoLep.size() > 1) && (event_ht > 120.);

		if (Wjets_selection) for (int i = 0; i<hlt_channels.size(); i++)
			{
			string hlt_channel = hlt_channels[i];
			string selection("wjets");

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);

			// RHO distrs
			fill_1d(hlt_channel + selection + string("_rho_passtrig_rawWeight"), 100, 0, 100, rho, rawWeight);
			fill_1d(hlt_channel + selection + string("_rho_passtrig_weight"),    100, 0, 100, rho, weight);

			// WJets, HLT channels
			if (JetHLTTrigger && MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   40, weight);
				}
			else if (JetHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   41, weight);
				}
			else if (MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   42, weight);
				}
			else continue;

			fill_1d(string(hlt_channel + "wjets_selection_ntaus"), 20, 0, 20,   jets.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJets"), 20, 0, 20,   selJets.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJetsNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

			fill_1d(string(hlt_channel + "wjets_selection_njets"), 20, 0, 20,   taus.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

			record_jets_fakerate_distrs(hlt_channel, selection, selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);
			record_jets_fakerate_distrs_large_bins(hlt_channel, selection, selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);
			}

		if (W1jet_selection) for (int i = 0; i<hlt_channels.size(); i++)
			{
			string hlt_channel = hlt_channels[i];
			string selection("w1jet");

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);

			// RHO distrs
			fill_1d(hlt_channel + selection + string("_rho_passtrig_rawWeight"), 100, 0, 100, rho, rawWeight);
			fill_1d(hlt_channel + selection + string("_rho_passtrig_weight"),    100, 0, 100, rho, weight);

			// WJets, HLT channels
			if (JetHLTTrigger && MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   50, weight);
				}
			else if (JetHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   51, weight);
				}
			else if (MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   52, weight);
				}
			else continue;

			fill_1d(string(hlt_channel + selection +"_selection_ntaus"), 20, 0, 20,   jets.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselJets"), 20, 0, 20,   selJets.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselJetsNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

			fill_1d(string(hlt_channel + selection +"_selection_njets"), 20, 0, 20,   taus.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
			fill_1d(string(hlt_channel + selection +"_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

			record_jets_fakerate_distrs(hlt_channel, selection, selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);
			record_jets_fakerate_distrs_large_bins(hlt_channel, selection, selJetsNoLep, selTausNoLep, visible_gen_taus, weight, isMC);
			}


		if (QCD_selection) for (int i = 0; i<hlt_channels.size(); i++)
			{
			string hlt_channel = hlt_channels[i];
			string selection("qcd");

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_rawWeight"), 100, 0, 100, vtx.size(), rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nvtx_weight"),    100, 0, 100, vtx.size(), weight);

			// pu distrs
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_rawWeight"), 100, 0, 100, nGoodPV, rawWeight);
			fill_1d(hlt_channel + selection + string("_pileup_passtrig_nGoodPV_weight"),    100, 0, 100, nGoodPV, weight);

			// RHO distrs
			fill_1d(hlt_channel + selection + string("_rho_passtrig_rawWeight"), 100, 0, 100, rho, rawWeight);
			fill_1d(hlt_channel + selection + string("_rho_passtrig_weight"),    100, 0, 100, rho, weight);

			// In QCD selection we need to skip the biasing triggered jet
			/*
			// to unpack the TriggerNames in the triggerobject:
			pat::TriggerObjectStandAlone& obj = trig_objs[closest_trigger_object_i];
			obj.unpackPathNames(trigNames);
			// maybe one can loop and unpack all names beforehand?
			*/

			if (debug)
				cout << "selecting " << jetHLT << " trigger objects" << endl;

			// select only our_hlt_trigger_objects
			std::vector<pat::TriggerObjectStandAlone> our_hlt_trigger_objects;
			for (size_t i = 0; i < trig_objs.size(); i++)
				{
				pat::TriggerObjectStandAlone& obj = trig_objs[i];
				obj.unpackPathNames(trigNames);

				bool is_our_hlt = false;
				for (unsigned h = 0; h < obj.pathNames().size(); ++h)
					{
					// HLT_PFJet40_v
					is_our_hlt |= (obj.pathNames()[h].find(jetHLT) != std::string::npos);
					}

				if (is_our_hlt)
					our_hlt_trigger_objects.push_back(obj);
				}

			// if a jet matches to our_hlt_trigger_objects add it to separate collection
			// afterwards there is a procedure removing the "leading jet" from this collection
			// thus the "leading triggering jet" is removed
			// --- it's done due to presence of many triggering jets (matching to our jetHLT)
			//     and no 1 jet can be removed simply due to it matching to trigger
			pat::JetCollection probeJets;
			pat::JetCollection probeJets_our_hlt;
			std::vector<pat::TriggerObjectStandAlone> our_hlt_probeJets;
			//for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
				{
				pat::Jet& jet = selJetsNoLep[ijet];

				// find trigger object closest to the jet
				double minDRtj (9999.);
				int closest_trigger_object_i = 0;
				for (size_t i = 0; i < our_hlt_trigger_objects.size(); i++)
					{
					double jet_trig_distance = TMath::Min(minDRtj, reco::deltaR (jet, our_hlt_trigger_objects[i]));
					if (jet_trig_distance<minDRtj)
						{
						closest_trigger_object_i = i;
						minDRtj = jet_trig_distance;
						}
					}

				// jet matches to trigger object, if dR < 0.3
				// if a jet doesn't match to any object -- add it to probeJets
				if (minDRtj > 0.3)
					{
					probeJets.push_back(jet);
					continue;
					}
				// in case there are no our trigger objects
				// this will catch all jets

				// the jet matches to an object of our trigger within minimum dR
				// save the jet
				// and also save the matching object -- in case it is needed for "leading" jet procedure
				pat::TriggerObjectStandAlone& obj = our_hlt_trigger_objects[closest_trigger_object_i];
				probeJets_our_hlt.push_back(jet);
				our_hlt_probeJets.push_back(obj);

				// TODO: remove this outdated comment
				// jets that match to trigger objects
				// are probeJets, if the trigger objects don't correspond to the chosen HLT trigger
				// jets matching to the chosen HLT trigger will bias the distributions
				// (with HLT_PFJet40_* the jet and tau distrs have a bump at about 40GeV)
				// thus the jet corresponding to the HLT should be removed
				// now the procedure is:
				//   skip the jet if it is only 1 matching to HLT
				//   if many jets match to HLT -- add them all to probeJets
				}

			// if there is only 1 jet matching to our HLT -- skip it
			// otherwise -- add all of them to the probeJets
			// DIDNT WORK
			/*
			if (probeJets_our_hlt.size() > 1)
				for (int i = 0; i<probeJets_our_hlt.size(); i++) probeJets.push_back(probeJets_our_hlt[i]);
			*/
			//std::sort (probeJets_our_hlt.begin(),  probeJets_our_hlt.end(),  utils::sort_CandidatesByPt);
			//don't sort -- need max pt of matched trigger object
			// now the highest-pT is first

			// anyway, test it:
			/*
			if (debug && probeJets_our_hlt.size() > 0)
				{
				cout << "QCD selection, our-HLT jets " << probeJets_our_hlt.size();
				cout << ", largest pT, smallest pT [check]" << endl;
				//cout << probeJets_our_hlt[0].pt() << "\t" << probeJets_our_hlt[probeJets_our_hlt.size()-1].pt() << endl;
				cout << probeJets_our_hlt[0] << "\t" << probeJets_our_hlt[probeJets_our_hlt.size()-1] << endl;
				}

			// skip the smallest-pT jet of our-HLT jets
			if (probeJets_our_hlt.size() > 1)
				{
				for (int i = 0; i<probeJets_our_hlt.size()-1; i++) probeJets.push_back(probeJets_our_hlt[i]);
				}
			// also record the pT and eta distr of the skipped jet:
			if (probeJets_our_hlt.size() > 0)
				{
				fill_1d(hlt_channel + string("skiped_jet_pt"),  200, 0, 200,   probeJets_our_hlt[probeJets_our_hlt.size()-1].pt(),  weight);
				fill_1d(hlt_channel + string("skiped_jet_eta"), 200, 0, 200,   probeJets_our_hlt[probeJets_our_hlt.size()-1].eta(), weight);
				}
			*/
			// it didn't work -- it removes the whole slope at low pT with small part of the bump

			// let's try removing 2-nd smallest jet
			// [worked pretty well, but not perfect]
			/*
			if (probeJets_our_hlt.size() > 1)
				{
				// all jets up to one-before-last
				for (int i = 0; i<probeJets_our_hlt.size()-2; i++) probeJets.push_back(probeJets_our_hlt[i]);
				// last jet
				probeJets.push_back(probeJets_our_hlt[probeJets_our_hlt.size()-1]);
				}
			// else {if there is only 1 our_hlt jet -- skip it}
			// also record the pT and eta distr of the skipped jet:
			if (probeJets_our_hlt.size() > 1)
				{
				fill_1d(hlt_channel + string("skiped_jet_pt"),  200, 0, 200,   probeJets_our_hlt[probeJets_our_hlt.size()-2].pt(),  weight);
				fill_1d(hlt_channel + string("skiped_jet_eta"), 200, 0, 200,   probeJets_our_hlt[probeJets_our_hlt.size()-2].eta(), weight);
				}
			*/

			// skip leading pt jet of our-hlt jets
			if (probeJets_our_hlt.size() == 1)
				{
				fill_1d(hlt_channel + selection + string("_skipped_jet_pt"),  200, 0, 200,   probeJets_our_hlt[0].pt(),  weight);
				fill_1d(hlt_channel + selection + string("_skipped_jet_energy"),  200, 0, 200,   probeJets_our_hlt[0].energy(),  weight);
				fill_1d(hlt_channel + selection + string("_skipped_jet_eta"), 200, -2.5, 2.5,   probeJets_our_hlt[0].eta(), weight);
				}

			if (probeJets_our_hlt.size() > 1)
				{
				// find "leading" trig object (the one to remove):
				int leading_jet_n = 0;

				/*
				double leading_jet = our_hlt_probeJets[0].pt();
				for (int i = 1; i<our_hlt_probeJets.size(); i++)
					{
					if (our_hlt_probeJets[i].pt() > leading_jet)
						{
						leading_jet = our_hlt_probeJets[i].pt();
						leading_jet_n = i;
						}
					}
				*/

				// returning back to removing highest-reco-pT jet
				double leading_jet = probeJets_our_hlt[0].pt();
				for (int i = 1; i<probeJets_our_hlt.size(); i++)
					{
					if (probeJets_our_hlt[i].pt() > leading_jet)
						{
						leading_jet = probeJets_our_hlt[i].pt();
						leading_jet_n = i;
						}
					}

				// removing highest-reco-E jet
				// worse then highest-reco-pT
				/*
				double leading_jet = probeJets_our_hlt[0].energy();
				for (int i = 1; i<probeJets_our_hlt.size(); i++)
					{
					if (probeJets_our_hlt[i].energy() > leading_jet)
						{
						leading_jet = probeJets_our_hlt[i].energy();
						leading_jet_n = i;
						}
					}
				*/

				// and skip the jet corresponding to this object
				// record the pT and eta distr of the skipped jet:
				fill_1d(hlt_channel + selection + string("_skipped_jet_pt"),  200, 0, 200,   probeJets_our_hlt[leading_jet_n].pt(),  weight);
				fill_1d(hlt_channel + selection + string("_skipped_jet_energy"),  200, 0, 200,   probeJets_our_hlt[leading_jet_n].energy(),  weight);
				fill_1d(hlt_channel + selection + string("_skipped_jet_eta"), 200, -2.5, 2.5,   probeJets_our_hlt[leading_jet_n].eta(), weight);
				// all jets, except the leading jet (considered the main trigger jet)
				for (int i = 0; i<probeJets_our_hlt.size(); i++)
					{
					if (i != leading_jet_n) probeJets.push_back(probeJets_our_hlt[i]);
					}
				}


			// QCD, HLT channels
			if (JetHLTTrigger && MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   60, weight);
				}
			else if (JetHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   61, weight);
				}
			else if (MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   62, weight);
				}
			else continue;

			fill_1d(string(hlt_channel + "qcd_selection_njets"), 20, 0, 20,   jets.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJets"), 20, 0, 20,   selJets.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJetsNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nprobeJets"), 20, 0, 20,   probeJets.size(), weight);

			fill_1d(string(hlt_channel + "qcd_selection_ntaus"), 20, 0, 20,   taus.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTaus"), 20, 0, 20,   selTaus.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTausNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTausNoLepNoJet"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

			//record_jets_fakerate_distrs(hlt_channel, selection, selJetsNoLep, selTausNoLep, weight, isMC);
			record_jets_fakerate_distrs(hlt_channel, selection, probeJets, selTausNoLep, visible_gen_taus, weight, isMC);
			record_jets_fakerate_distrs_large_bins(hlt_channel, selection, probeJets, selTausNoLep, visible_gen_taus, weight, isMC);
			//int record_jets_fakerate_distrs(string & channel, string & selection, pat::JetCollection & selJets, pat::TauCollection & selTaus)
			}

		/*
		measure:
		tau PT VS tau-jet PT
		N "distinct taus" for jet-match with NoLep jets, all jets, ID jets, Nolep jets with lower PT reuirement

		the "jets - taus" selJetsNoLepNoTau are corrected with this tau fakerate
		*/

		} // End single file event loop

	delete file;
	} // End loop on files

printf("Done processing the job of files\n");

printf("End of (file loop) the job.\n");

// Controls distributions of processed particles


// CONTROLINFO

/*
FILE *csv_out;
string FileName = ((outUrl.ReplaceAll(".root",""))+".csv").Data();
csv_out = fopen(FileName.c_str(), "w");

fprintf(csv_out, "New output (sums per whole job!):\n");

//printout_counters(csv_out, string(isMC ? "MC,": "Data,") + dtag_s + string(",") + job_num);
//printout_distrs(csv_out, string(isMC ? "MC,": "Data,") + dtag_s + string(",") + job_num);
//printout_counters(csv_out, job_def);
//printout_distrs(csv_out, job_def);

// So, each job output contains for each value:
// value_name,MC/Data,dtag,job_num,value
// and the table for distr:
// distr_name,MC/Data,dtag,job_num,bin,contents
// simple grep/ack by first coloumn produces a ready data.frame for R
// which should be summarised for each dataset -- easy with ddply
// the only problem is:
// FIXME: how to do MULTISELECT nicely here?

fprintf(csv_out, "End of job output.\n\n");

fclose(csv_out);
*/

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
//########     SAVING HISTO TO FILE     ########
//##############################################
//save control plots to file
printf ("Results save in %s\n", outUrl.Data());

/*
// re-enabled the ROOT output
// for resubmit option of the job script
TFile *ofile = TFile::Open (outUrl + ".root", "recreate");

//singlelep_ttbar_initialevents->Write();
//singlelep_ttbar_preselectedevents->Write();

//qcd_jets_distr->Write();
//qcd_tau_jets_distr->Write();
//wjets_jets_distr->Write();
//wjets_tau_jets_distr->Write();
//wjets_distinct_tau_distr->Write();
//qcd_distinct_tau_distr->Write();

qcd_taujet_origin->Write();
qcd_jet_origin->Write();
wjets_taujet_origin->Write();
wjets_jet_origin->Write();

qcd_taujet_distance->Write();
wjets_taujet_distance->Write();

//for(std::map<string, TH1D>::iterator it = th1d_distr_control_headers.begin(); it != th1d_distr_control_headers.end(); ++it)
//	{
//	string name = it->first;
//	TH1D * header_distr = & it->second;
//	// Header:
//	fprintf(out, "header,%s", name.c_str());
//	for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
//	fprintf(out, "\n");
//	}

// write all jet distr histos

for(std::map<std::pair <string,string>, TH3F>::iterator it = th3f_distr_control.begin(); it != th3f_distr_control.end(); ++it)
	{
	const std::pair <string,string> *key = &it->first;
	string mc_decay_suffix = key->first;
	string name = key->second;

	TH3F * distr = & it->second;

	distr->Write();
	}

//th2d_distr_control
//fill_pt_pt(string control_point_name, double pt1, double pt2, double weight)
for(std::map<std::pair <string,string>, TH2D>::iterator it = th2d_distr_control.begin(); it != th2d_distr_control.end(); ++it)
	{
	const std::pair <string,string> *key = &it->first;
	string mc_decay_suffix = key->first;
	string name = key->second;

	TH2D * distr = & it->second;

	distr->Write();
	}



ofile->Close();
*/

for(std::map<string, std::map<string, TH1D>>::iterator it = th1d_distr_maps_control.begin(); it != th1d_distr_maps_control.end(); ++it)
	{
	// const std::pair <string,string> *key = &it->first;
	string channel = it->first;

	//outUrl.Data() is dtag_jobnum
	// use them separately, take from: dtag_s, job_num
	// TFile* out_f = TFile::Open (TString(outUrl.Data() + string("_") + channel + string(".root")), "CREATE");

	TString filename = outdir + TString(string("/") + dtag_s + channel + string("_") + job_num + string(".root"));
	cout << "opening " << filename << endl;
	TFile* out_f = TFile::Open(filename, "CREATE");
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
	cout << "closed " << endl;
	}

FILE *csv_out;
string FileName = ((outUrl.ReplaceAll(".root",""))+".job_done").Data();
csv_out = fopen(FileName.c_str(), "w");

if (outTxtFile) fclose (outTxtFile);

// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}


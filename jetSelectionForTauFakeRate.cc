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
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
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

bool passPFJetID(std::string label, pat::Jet jet)
	{
	// Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

	bool passID(false); 

	// float rawJetEn(jet.correctedJet("Uncorrected").energy() );

	double eta=jet.eta();
 
 	// from the twiki:
	float NHF = jet.neutralHadronEnergyFraction();
	float NEMF = jet.neutralEmEnergyFraction();
	float CHF = jet.chargedHadronEnergyFraction();
	float MUF = jet.muonEnergyFraction();
	float CEMF = jet.chargedEmEnergyFraction();
	float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	float NumNeutralParticles =jet.neutralMultiplicity();
	float CHM = jet.chargedMultiplicity(); 
	// TODO: check if these change corresponding to jet corrections? Or apply corrections after passing the selection?

	// float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
	// float nef( jet.neutralEmEnergy()/rawJetEn );
	// float cef( jet.chargedEmEnergy()/rawJetEn );
	// float chf( jet.chargedHadronEnergy()/rawJetEn );
	// float nch    = jet.chargedMultiplicity();
	// float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	// float muf(jet.muonEnergy()/rawJetEn); 

	if (label=="Loose")
		// passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
		// the same, just with new names:
		passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;
	if (label=="Tight")
		// passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
		passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;

	// there is also:
	//tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=3.0

	// And the abs(eta)>3.0 part, but we never consider such jets, so... Meh!

	return passID; 
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

string mc_decay("");
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

// good bins 1
// Float_t bins_pt[11] = { 0, 29, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
Float_t bins_pt[11] = { 0, 30, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
Float_t bins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
int n_bins_eta = 7;
Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges

//TH3F* wjets_jets_distr      = (TH3F*) new TH3F("wjets_jets_distr",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad);

std::map<string, std::map<string, TH1D>> th1d_distr_maps_control;
std::map<string, TH1D> th1d_distr_maps_control_headers;

std::map<string, std::map<string, TH1I>> th1i_distr_maps_control;
std::map<string, TH1I> th1i_distr_maps_control_headers;

std::map<string, std::map<string, TH2D>> th2d_distr_maps_control;
std::map<string, TH2D> th2d_distr_maps_control_headers;

std::map<string, std::map<string, TH3D>> th3d_distr_maps_control;
std::map<string, TH3D> th3d_distr_maps_control_headers;





int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th1d_distr_maps_control
	if (th1d_distr_maps_control.find(mc_decay) == th1d_distr_maps_control.end() )
		{
		//
		th1d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH1D>()));
		}

	if (th1d_distr_maps_control[mc_decay].find(control_point_name) == th1d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th1d_distr_maps_control[mc_decay][control_point_name] = TH1D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_maps_control[mc_decay][control_point_name].Fill(value, weight);

	if (th1d_distr_maps_control_headers.find(control_point_name) == th1d_distr_maps_control_headers.end() )
		{
		// th1d_distr_maps_control_headers[control_point_name] = TH1D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th1d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH1D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup)));
		}

	// return success:
	return 0;
	}



int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th1i_distr_maps_control
	if (th1i_distr_maps_control.find(mc_decay) == th1i_distr_maps_control.end() )
		{
		//
		th1i_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH1I>()));
		}

	if (th1i_distr_maps_control[mc_decay].find(control_point_name) == th1i_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th1i_distr_maps_control[mc_decay][control_point_name] = TH1I((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1i_distr_maps_control[mc_decay][control_point_name].Fill(value, weight);

	if (th1i_distr_maps_control_headers.find(control_point_name) == th1i_distr_maps_control_headers.end() )
		{
		// th1i_distr_maps_control_headers[control_point_name] = TH1I("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th1i_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH1I((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup)));
		}

	// return success:
	return 0;
	}



int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th2d_distr_maps_control
	if (th2d_distr_maps_control.find(mc_decay) == th2d_distr_maps_control.end() )
		{
		//
		th2d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH2D>()));
		}

	if (th2d_distr_maps_control[mc_decay].find(control_point_name) == th2d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th2d_distr_maps_control[mc_decay][control_point_name] = TH2D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th2d_distr_maps_control[mc_decay][control_point_name].Fill(x, y, weight);

	if (th2d_distr_maps_control_headers.find(control_point_name) == th2d_distr_maps_control_headers.end() )
		{
		// th2d_distr_maps_control_headers[control_point_name] = TH2D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th2d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH2D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup)));
		}

	// return success:
	return 0;
	}


int fill_3d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, double x, double y, double z, double weight)
	{
	// check if control point has been initialized
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
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup );
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(x, y, z, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup)));
		}

	// return success:
	return 0;
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

	/*
	if (th3f_distr_control.find(key) == th3f_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		// th3f_distr_control[key] = TH1D(control_point_name.c_str(), ";;ID", 600, -300., 300.);
		//cout << "creating " << mc_decay << " - " << control_point_name << endl;
		// th3f_distr_control.insert( std::make_pair(key, TH1D((mc_decay + control_point_name).c_str(), ";;ID", 600, -300., 300.)));
		th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		//cout << "creating " << control_point_name << endl;
		}
	*/

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad);
		// th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}


	// fill the distribution:
	// th1i_distr_control[key].Fill(value, weight);
	// th3f_distr_control[key].Fill(value, weight);
	// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
	// th3f_distr_control[key].Fill(pt, eta, radius, weight);
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, radius, weight);

	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	/*
	if (th3f_distr_control_headers.find(string("j_distr")) == th3f_distr_control_headers.end() )
		{
		// th3f_distr_control_headers[string("p_id")] = TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.)));
		th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		}
	*/
	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
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
TRandom *r3 = new TRandom3();

// configure the process
const edm::ParameterSet & runProcess = edm::readPSetsFrom (argv[1])->getParameter < edm::ParameterSet > ("runProcess");

bool debug           = runProcess.getParameter<bool>  ("debug");
bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
bool isMC            = runProcess.getParameter<bool>  ("isMC");
double xsec          = runProcess.getParameter<double>("xsec");
int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
TString dtag         = runProcess.getParameter<std::string>("dtag");
string dtag_s        = runProcess.getParameter<std::string>("dtag");
string job_num       = runProcess.getParameter<std::string>("job_num");

JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

	
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

bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );
	
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
// v1
// getJetCorrector(TString baseDir, TString pf, bool isMC)
FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector (jecDir, "/Spring16_25nsV6_", isMC);
//TString pf(isMC ? "MC" : "DATA");
//JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/"+pf+"_Uncertainty_AK4PFchs.txt").Data ());
//JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/" + (isMC ? "MC" : "DATA") + "_Uncertainty_AK4PFchs.txt").Data ());
// JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/Fall15_25nsV2_" + (isMC ? "MC" : "DATA") + "_Uncertainty_AK4PFchs.txt").Data ());
JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/Spring16_25nsV6_" + (isMC ? "MC" : "DATA") + "_Uncertainty_AK4PFchs.txt").Data ());
// JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/MC_Uncertainty_AK4PFchs.txt").Data ());

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
// TODO: move all these numbers to where they are applied??
// btagMedium is used twice in the code
// merge those tagging procedures and eliminated the variable?

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

// new btag calibration
// TODO: check callibration readers in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
// and latest standalone callibrator:
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h

// Setup calibration readers
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
// BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_76X.csv");
BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_ichep_80X.csv");
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
                             "central",             // central sys type
                             {"up", "down"});      // other sys types
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_B,      // btag flavour
            "mujets");               // measurement type
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_C,      // btag flavour
            "mujets");              // measurement type
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



// ----------------------------
// So here we got all the parameters from the config














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
			if(iev == 5)
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
		if (iev % treeStep == 0)
			{
			printf (".");
			if(!debug) fflush (stdout); // Otherwise debug messages are flushed
			}

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

		weight *= weight_PU;
		weight_up *= weight_PU_up;
		weight_down *= weight_PU_down;

		// weight PU 3
		fill_1d(string("weightflow"), 300, 0, 300,   3, weight);

		// pu distrs
		fill_1d( string("pileup_passtrig_rawWeight_pergoodpv"), 100, 0, 100, nGoodPV, rawWeight);
		fill_1d( string("pileup_passtrig_weight_pergoodpv"),    100, 0, 100, nGoodPV, weight);

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		//num_inters = 1;
		if (num_inters>99) num_inters = 99;
		if (nGoodPV>100) nGoodPV = 99;
		//if (num_inters<0)  num_inters = 0;
		if (weight_Gen<0)
			{
			increment( string("negative_events"), 1 );
			fill_pu( string("pileup_negative_weight_pernuminters"), num_inters, weight);
			fill_pu( string("pileup_negative_rawweight_pernuminters"), num_inters, rawWeight);
			}
		else
			{
			increment( string("positive_events"), 1 );
			fill_pu( string("pileup_positive_weight_pernuminters"), num_inters, weight);
			fill_pu( string("pileup_positive_rawweight_pernuminters"), num_inters, rawWeight);
			}

		// -------------------------------------   Basic event selection

		// ---------------------- Orthogonalize Run2015B PromptReco+17Jul15 mix
		// let's remove Run2015B
		// if(isRun2015B)
		// {
		// if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
		// }

		// it's not needed with the latest versions of RunB rereconstruction
		
		// -------------------------------------------------- Skip bad lumi
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
		fill_1d(string("weightflow"), 300, 0, 300,   5, weight);
		/*
		fill_1d(string("weightflow_mu"), 300, 0, 300,   5, weight);
		fill_1d(string("weightflow_el"), 300, 0, 300,   5, weight);
		fill_1d(string("weightflow_elel"), 300, 0, 300, 5, weight);
		fill_1d(string("weightflow_elmu"), 300, 0, 300, 5, weight);
		fill_1d(string("weightflow_mumu"), 300, 0, 300, 5, weight);
		*/

		// --------------------------------------------- apply trigger
		// ---------------- and require compatibilitiy of the event with the PD
		edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");

		if (!tr.isValid ()){
			//cout << "Trigger HLT is not valid, trying HLT2" << endl;
			tr = ev.triggerResultsByName ("HLT2");
			if (!tr.isValid ()){
				cout << "Trigger HLT2 is not valid, exiting" << endl;
				return false;
				}
			}

		if(debug){
			cout << "Printing HLT trigger list" << endl;
			//cout << "-- Commented out --" << endl;
			for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames)
				cout << *trnames << endl;
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
			JetHLTTrigger = utils::passTriggerPatterns(tr, "HLT_PFJet40_v*");
			MuonHLTTrigger = (isMC ?
				//utils::passTriggerPatterns(tr, "HLT_IsoMu22*", "HLT_IsoTkMu22*");
				utils::passTriggerPatterns (tr, "HLT_IsoMu24_v2", "HLT_IsoTkMu24_v2") : // tecommended inn ttbar trig for reHLT
				//utils::passTriggerPatterns (tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*") :
				utils::passTriggerPatterns (tr, "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*"));
			}

		// if (!(JetHLTTrigger || MuonHLTTrigger)) continue; // No orthogonalization -- run on only 1 trigger type of datasets

		string hlt_channel("");

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

		/*
		bool eTrigger    ( isMC ?
			//utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v2") :
			utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v*") :
			utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v*") );
		bool muTrigger   ( isMC ?
			//utils::passTriggerPatterns (tr, "HLT_IsoMu22_v3", "HLT_IsoTkMu22_v3") :
			utils::passTriggerPatterns (tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*") :
			utils::passTriggerPatterns (tr, "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*")
			// utils::passTriggerPatterns (tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*") // the efficiency scale factor are available for these only
			// utils::passTriggerPatterns (tr, "HLT_IsoMu18_v*", "HLT_IsoTkMu18_v*")
			// utils::passTriggerPatterns (tr, "HLT_IsoMu18_v*")
			);

		// if data and SingleElectron dataset and both triggers -- skip event
		if (!debug) {
			if (!isMC && isSingleElectronDataset && eTrigger && muTrigger) continue;
		
			if (!(eTrigger || muTrigger)) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
		}
		*/

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

		// passtrig 6
		fill_1d(string("weightflow"), 300, 0, 300,   6, weight);

		//increment( string("weightflow_weight_passed_trig"), weight ); // should not matter
		//increment( string("weightflow_weight_up_passed_trig"), weight_up ); // should not matter
		//increment( string("weightflow_weight_down_passed_trig"), weight_down ); // should not matter

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

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
		

		// passmetfilters 7
		fill_1d(string("weightflow"), 300, 0, 300,   7, weight);

		if(debug)
			{
			cout << "met filters applied here" << endl;
			}


		// ------------------------- event physics and the corresponding selection

		//------------------------- load all the objects we will need to access

		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;

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
		LorentzVector elDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton>
		pat::ElectronCollection selElectrons;
		unsigned int nVetoE(0);

		for(unsigned int count_idiso_electrons = 0, n=0; n<electrons.size (); ++n)
			{
			pat::Electron& electron = electrons[n];

			fill_2d(string("control_el_slimmedelectrons_pt_eta"), 250, 0., 500., 200, -3., 3., electron.pt(), electron.eta(), weight);
			fill_1d(string("control_el_slimmedelectrons_phi"), 128, -3.2, 3.2, electron.phi(), weight);
			//fill_1d(string("control_tau_slimmedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);

			bool 
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);
			bool passSigma(false), passSigmaVeto(false);
			// from passId( pat::Photon .. ) of PatUtils:
			//bool elevto = photon.hasPixelSeed();  //LQ  REACTIVATED FOR TIGHT ID, OTHERWISE MANY ELECtRONS pass the photon Id
			// and then, in Tight ID, they include:
			// && !elevto 
			//
			// So, it would be nice to add it to Tight Electron ID:
			//bool hasPixelSeed = electron.hasPixelSeed();
			// but electrons don't have this method
			// will have to cross-clean with photons or etc

			// removing all electrons close to tight Photons
			// actually, people do it the other way around, testing v9.5
			// nope, it doesn't fix anything
			/*
			double minDRlg(9999.);
			for(size_t i=0; i<selPhotons.size(); i++)
				{
				minDRlg = TMath::Min(minDRlg, deltaR(electron.p4(), selPhotons[i].p4()));
				}
			if(minDRlg<0.1) continue;
			*/

			int lid(electron.pdgId()); // should always be 11

			//apply electron corrections
			/* no lepcorrs in 13.6
			if(abs(lid)==11)
				{
				elDiff -= electron.p4();
				ElectronEnCorrector.calibrate(electron, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
				//electron = patUtils::GenericLepton(electron.el); //recreate the generic electron to be sure that the p4 is ok
				elDiff += electron.p4();
				}
			*/

			// TODO: probably, should make separate collections for each step, for corrected particles, then -- passed ID etc
			// fill_pt_e( string("all_electrons_corrected_pt"), electron.pt(), weight);
			// if (n < 2)
			// 	{
			// 	fill_pt_e( string("top2pt_electrons_corrected_pt"), electron.pt(), weight);
			// 	}

			//no need for charge info any longer
			//lid = abs(lid);
			//TString lepStr(lid == 13 ? "mu" : "e");
			// should always be 11!
					
			// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
			// no need to mess with photon ID // double minDRlg(9999.);
			// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
			// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[n].p4(),selPhotons[ipho].p4()));
			// no need to mess with photon ID // if(minDRlg<0.1) continue;


			// ------------------------- electron IDs
			//Cut based identification
			//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
			// full5x5_sigmaIetaIeta is forgotten in our passId for electrons
			// getting high MC/Data in electrons now -- maybe due to photons,
			// checking sigmaIetaIeta, then photon rejection (hasPixelSeed()) and cross-cleaning with photons
			float eta = std::abs(electron.superCluster()->position().eta());
			float sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();
			if (eta <= 1.479) // barrel, newer selection is precise?
				{
				passSigma =     sigmaIetaIeta < 0.00998; // Tight WP
				passSigmaVeto = sigmaIetaIeta < 0.0115;  // Veto WP
				}
			else if (eta > 1.479) // endcap
				{
				passSigma =     sigmaIetaIeta < 0.0292; // Tight WP
				passSigmaVeto = sigmaIetaIeta < 0.037;  // Veto WP
				}
			passId     = patUtils::passId(electron, goodPV, patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut) && passSigma;
			passVetoId = patUtils::passId(electron, goodPV, patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut) && passSigmaVeto;

			// ------------------------- electron isolation

			passIso     = patUtils::passIso(electron, patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, rho);
			passVetoIso = patUtils::passIso(electron, patUtils::llvvElecIso::Loose, patUtils::CutVersion::ICHEP16Cut, rho);


			// ---------------------------- kinematics
			//double leta(fabs(lid==11 ? lepton.el.superCluster()->eta() : lepton.eta()));
			double leta( electron.superCluster()->eta() );

			// ---------------------- Main lepton kin
			//if(electron.pt() < 30.)                passKin = false;
			if(electron.pt() < 35.)                passKin = false;
			// if(leta > 2.1)                         passKin = false;
			if(leta > 2.4)                         passKin = false;
			if(leta > 1.4442 && leta < 1.5660)     passKin = false; // Crack veto

			// ---------------------- Veto lepton kin
			// if (electron.pt () < 20)            passVetoKin = false;
			if (electron.pt () < 15)            passVetoKin = false;
			// if (leta > 2.1)                     passVetoKin = false;
			if (leta > 2.5)                     passVetoKin = false;
			if (leta > 1.4442 && leta < 1.5660) passVetoKin = false; // Crack veto


			if     (passKin     && passId     && passIso)
				{
				selElectrons.push_back(electron);
				fill_2d(string("control_el_selElectrons_pt_eta"), 250, 0., 500., 200, -3., 3., electron.pt(), electron.eta(), weight);
				fill_1d(string("control_el_selElectrons_phi"), 128, -3.2, 3.2, electron.phi(), weight);
				}
			else if(passVetoKin && passVetoId && passVetoIso) nVetoE++;

			}

		// TODO: there should be no need to sort selected electrons here again -- they are in order of Pt
		std::sort(selElectrons.begin(),   selElectrons.end(),   utils::sort_CandidatesByPt);


		if(debug){
			cout << "processed electrons" << endl;
			}





		// ---------------------------------- MUONS SELECTION
		LorentzVector muDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton> selLeptons;
		pat::MuonCollection selMuons;
		unsigned int nVetoMu(0);
		// unsigned int count_idiso_muons = 0;

		for(unsigned int count_idiso_muons = 0, n=0; n<muons.size (); ++n)
			{
			// patUtils::GenericLepton& lepton = leptons[n];
			pat::Muon& muon = muons[n];

			fill_2d(string("control_mu_slimmedmions_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
			fill_1d(string("control_mu_slimmedmions_phi"), 128, -3.2, 3.2, muon.phi(), weight);

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
			passId     = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdTight, patUtils::CutVersion::ICHEP16Cut);
			passVetoId = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdLoose, patUtils::CutVersion::ICHEP16Cut);

			// ------------------------- muon isolation
			passIso     = patUtils::passIso(muon, patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
			passVetoIso = patUtils::passIso(muon, patUtils::llvvMuonIso::Loose, patUtils::CutVersion::ICHEP16Cut);

			if (passId && passIso)
				{
				fill_pt_e( string("all_muons_idiso_pt"), muon.pt(), weight);
				if (count_idiso_muons < 2) fill_pt_e( string("top2pt_muons_idiso_pt"), muon.pt(), weight);
				count_idiso_muons += 1;
				}


			// ---------------------------- kinematics
			double leta( muon.eta());

			// ---------------------- Main muon kin
			if(muon.pt() < 30.)   passKin = false;
			// if(leta > 2.1)        passKin = false;
			//if(muon.pt() < 26.)   passKin = false;
			if(leta > 2.4)        passKin = false;

			// ---------------------- Veto muon kin
			//if (muon.pt () < 20)  passVetoKin = false;
			// if (leta > 2.1)       passVetoKin = false;
			if (muon.pt () < 10)  passVetoKin = false;
			if (leta > 2.5)       passVetoKin = false;

			if     (passKin     && passId     && passIso)
				{
				selMuons.push_back(muon);
				fill_2d(string("control_mu_selMuons_pt_eta"), 250, 0., 500., 200, -3., 3., muon.pt(), muon.eta(), weight);
				fill_1d(string("control_mu_selMuons_phi"), 128, -3.2, 3.2, muon.phi(), weight);
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

		if(debug){
			cout << "processed muons" << endl;
			}


		// Check for the single-electron/single-muon channels and save Pt-s if the event is in channel

		// Here we can already assign leptonic channel
		// electron or muon -- tau
		// the selection is:

		// ------------------------------------- Propagate lepton energy scale to MET
		// Propagate now (v13)
		/* no lepton corrections propagation v13.1
		met.setP4(met.p4() - muDiff - elDiff); // TODO: note this also propagates to all MET uncertainties -- does it??
		met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
		*/

		if(debug){
			cout << "NOT propagated lepton corrections to met" << endl;
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
		pat::TauCollection selTaus;
		for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
			{
			pat::Tau& tau = taus[n];

			fill_2d(string("control_tau_slimmedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
			fill_1d(string("control_tau_slimmedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);

			// ---------- IDs
					
			// if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
			// if(tau.emFraction() >=2.) continue;
					
			// Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
			// "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
			//if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
			if (tau.tauID("decayModeFinding")<0.5) continue; // OldDMs is synonim for no <DMs>
			// Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
			// Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM
			// if (tau.tauID ("byMediumCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
			// if (tau.tauID ("byTightCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
			if (tau.tauID ("byMediumIsolationMVArun2v1DBoldDMwLT")<0.5) continue;

			if (tau.tauID ("againstMuonTight3")                          <0.5) continue; // Medium working point not usable. Available values: Loose, Tight
			//if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
			// if (tau.tauID ("againstElectronMediumMVA6")                  <0.5) continue;
			if (tau.tauID ("againstElectronTightMVA6")                  <0.5) continue;

			// maybe use MVA-based discriminators:
			// byVLooseIsolationMVArun2v1DBoldDMwLT
			// byLooseIsolationMVArun2v1DBoldDMwLT
			// byMediumIsolationMVArun2v1DBoldDMwLT
			// byTightIsolationMVArun2v1DBoldDMwLT
			// byVTightIsolationMVArun2v1DBoldDMwLT

			// or the new discriminators:
			// byLooseCombinedIsolationDeltaBetaCorr3HitsdR03
			// byMediumCombinedIsolationDeltaBetaCorr3HitsdR03
			// byTightCombinedIsolationDeltaBetaCorr3HitsdR03
			// -- recommended for multi-object final states (ttH, H->tau-tau)
			// -- not found in noHLT TTbar

			fill_2d(string("control_tau_idedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
			fill_1d(string("control_tau_idedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
			//fill_pt_e( string("all_taus_4discrs_pt"), tau.pt(), weight);

			// Pixel hits cut (will be available out of the box in new MINIAOD production)
			//{
			//int nChHadPixelHits = 0;
			//reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
			//for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++)
			//	{
			//	pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
			//	int pixelHits = packedCand->numberOfPixelHits();
			//	if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
			//	}
			//if(nChHadPixelHits==0) continue;
			//}

			/*
			fill_pt_e( string("all_taus_ided_pt"), tau.pt(), weight);
			if (count_ided_taus<1)
				{
				fill_pt_e( string("top1pt_taus_ided_pt"), tau.pt(), weight);
				count_ided_taus += 1;
				}
			*/

			// --------- Kinematics
			if (tau.pt() < 20. || fabs (tau.eta()) > 2.3) continue;

			selTaus.push_back(tau);
			fill_2d(string("control_tau_selTaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
			fill_1d(string("control_tau_selTaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
			}

		std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);


		if(debug){
			cout << "selected taus [individual]" << endl;
			}


		// ------------------------------------------ select the taus cleaned from leptons

		pat::TauCollection selTausNoLep;
		int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
		for (size_t itau = 0; itau < selTaus.size(); ++itau)
			{
			pat::Tau& tau = selTaus[itau];
			if(debug) cout << "selecting NoLep taus " << itau << endl;

			// cross-cleaning taus with leptons
			bool overlapWithLepton(false);
			for(int l=0; l<(int)selLeptons.size();++l)
				{
				if (reco::deltaR(tau, selLeptons[l])<0.4)
					{ overlapWithLepton=true; break; }
				}
			if (overlapWithLepton) continue;

			selTausNoLep.push_back(tau);
			// so these are the final taus we use in the selection
			fill_2d(string("control_tau_selTausNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
			fill_1d(string("control_tau_selTausNoLep_phi"), 128, -3.2, 3.2, tau.phi(), weight);

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			//increment(string("number_of_tausnolep_found"), 1);

			// for MC find the generated candidate closest to the tau
			// to get the fakeness
			/*
			if (isMC)
				{
				double min_deltaR = 99999.9;
				for(size_t i = 0; i < gen.size(); ++ i)
					{
					const reco::GenParticle & p = gen[i];
					double deltaR_to_p = reco::deltaR(tau, p);
					if (deltaR_to_p < min_deltaR)
						{
						min_deltaR = deltaR_to_p;
						closest_totaunolep_particle_id = p.pdgId();
						}
					//int id = p.pdgId();
					//int st = p.status();
					//int n = p.numberOfDaughters();
					//cout << i << ": " << id << " " << st;
					//if (p.numberOfMothers() != 0) cout << " <- " ;
					//for (int j = 0 ; j < p.numberOfMothers(); ++j)
					//	{
					//	const reco::Candidate * mom = p.mother(j);
					//	cout << " " << mom->pdgId() << " " << mom->status() << ";";
					//	}
					//cout << "\n";
					//if (n>0)
					//	{
					//	cout << "\t|-> " ;
					//	for (int j = 0; j < n; ++j)
					//		{
					//		const reco::Candidate * d = p.daughter( j );
					//		cout << d->pdgId() << " " << d->status() << "; " ;
					//		}
					//	cout << "\n";
					//	}
					}
				fill_particle_ids(string("nearest_particle_around_taunolep"), closest_totaunolep_particle_id, weight);
				// electron-tau fake-rate scale factor
				if (fabs(closest_totaunolep_particle_id)==11)
					{
					increment(string("number_of_tausnolep_from_electron_found"), 1);
					// apply the data/MC fake rate scale factor
					// (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Electron_to_tau_fake_rate)
					// for tight working points it is 1.80 +- 0.23 in barrel (eta < 1.460)
					// and 1.3 +- 0.42 in endcaps (eta > 1.558)
					// repeating 13.4, no electron-to-tau sf
					if (fabs(tau.eta()) < 1.460)
						{
						weight *= 1.8 + r3->Gaus(0, 0.2);  // gaussian +- 0.2
						}
					else if (fabs(tau.eta()) > 1.558)
						{
						weight *= 1.3 + r3->Gaus(0, 0.42); // gaussian +- 0.42
						}
					//
					// TODO: and then I need to renormalize the MC integral??
					}

				if (fabs(closest_totaunolep_particle_id)==13)
					{
					increment(string("number_of_tausnolep_from_muon_found"), 1);
					// TODO: add muon-tau fake rate?
					}
				}
			*/
			}

		if (isMC && selTausNoLep.size() > 0) // 2016 data/MC tau ID efficiency (all discriminators, all pt and eta ranges) = 0.83 +- 0.06
			{
			double pre_weight_tauIDsf = 1;
			for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
				//pre_weight_tauIDsf *= (1 - 0.83 + r3->Gaus(0, 0.06)); // gaussian +- 0.06
				// recommendations update (noticed 8-11-2016, last page update 2016-11-07)
				// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_ID_efficiency
				pre_weight_tauIDsf *= (1 - 0.9 + r3->Gaus(0, 0.1)); // gaussian +- 0.1
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



		if(debug){
			cout << "processed taus" << endl;
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

		//void updateJEC(pat::JetCollection& jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC){

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )


		if(debug) cout << "jet eta pt e, e x y z" << endl;

		// v6, adding jet corrections and b-tagging
		LorentzVector full_jet_corr(0., 0., 0., 0.);
		for(size_t ijet=0; ijet<jets.size(); ijet++)
			{
			// TODO: so does this mean "in place"?
			pat::Jet& jet = jets[ijet];
			//fill_2d(string("slimmedjet_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_2d(string("control_jet_slimmedjet_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_slimmedjet_phi"), 128, -3.2, 3.2, jet.phi(), weight);

			LorentzVector jet_corr(0., 0., 0., 0.);

			LorentzVector jet_initial_momentum = jet.p4();

			if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
			//correct JES
			LorentzVector rawJet = jet.correctedP4("Uncorrected");

			fill_2d(string("control_jet_slimmedjet_uncorrected_pt_eta"), 400, 0., 400., 200, -4., 4., rawJet.pt(), rawJet.eta(), weight);

			if(debug) cout << rawJet.eta() << " " << rawJet.pt() << " " << rawJet.energy() << endl;

			//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
			//LorentzVector rawJet(jet*toRawSF);
			// 13.8_2 trying all jet corrections + new (if it is new) jetID procedure
			jesCor->setJetEta(rawJet.eta());
			jesCor->setJetPt(rawJet.pt());
			jesCor->setJetA(jet.jetArea());
			jesCor->setRho(rho);
			jesCor->setNPV(nGoodPV);
			float jes_correction = jesCor->getCorrection();
			jet.setP4(rawJet*jes_correction);

			fill_2d(string("control_jet_slimmedjet_jescor_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_slimmedjet_jescorrection"), 400, 0., 2., jes_correction, weight);

			if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

			//smear JER
			//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
			//double newJERSF(1.0);
			// 13.8_2 trying all jet corrections + new (if it is new) jetID procedure
			if(isMC)
				{
				const reco::GenJet* genJet=jet.genJet();
				if(genJet)
					{
					double genjetpt( genJet ? genJet->pt(): 0.);                    
					std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
					double jer_smearing = smearJER[0];
					jet.setP4(jet.p4()*jer_smearing);
					fill_1d(string("control_jet_slimmedjet_mc_jersmearing"), 400, 0., 2., jer_smearing, weight);
					
					//printf("jet pt=%f gen pt = %f smearing %f %f %f\n", jet.pt(), genjetpt, smearJER[0], smearJER[1], smearJER[2]);
					// //set the JER up/down alternatives
					jet.addUserFloat("jerup", smearJER[1]);  //kept for backward compatibility
					jet.addUserFloat("jerdown", smearJER[2] ); //kept for backward compatibility
					jet.addUserFloat("_res_jup", smearJER[1]);
					jet.addUserFloat("_res_jdown", smearJER[2] );
					}
				else{
					jet.addUserFloat("jerup", 1.0); //kept for backward compatibility
					jet.addUserFloat("jerdown", 1.0);  //kept for backward compatibility
					jet.addUserFloat("_res_jup", 1.0);
					jet.addUserFloat("_res_jdown", 1.0 );
					}
				fill_2d(string("control_jet_slimmedjet_mc_jercor_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
				}


			// here is the correct3 jet correction point

			if(isMC)
				{
				////set the JES up/down pT alternatives
				std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
				jet.addUserFloat("jesup",    ptUnc[0] );  //kept for backward compatibility
				jet.addUserFloat("jesdown",  ptUnc[1] );  //kept for backward compatibility
				jet.addUserFloat("_scale_jup",    ptUnc[0] );
				jet.addUserFloat("_scale_jdown",  ptUnc[1] );
				}

			// FIXME: this is not to be re-set. Check that this is a desired non-feature.
			// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
			//to get the raw jet again
			//jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));

			// Add the jet momentum correction:
			// jet_cor propagation is on in 13.4
			jet_corr = jet.p4() - jet_initial_momentum;
			full_jet_corr += jet_corr;

			fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_eta"), 400, 0., 400., 200, -4., 4., full_jet_corr.pt(), full_jet_corr.eta(), weight);

			if(debug)
				{
				cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
				cout << "-" << jet_initial_momentum.eta() << " " << jet_initial_momentum.pt() << " " << jet_initial_momentum.energy() << endl;
				}
			}


		std::sort (jets.begin(),  jets.end(),  utils::sort_CandidatesByPt);

		// ----------------------------------- here is the correctF jet correction point
		// Propagate full_jet_corr to MET:
		met.setP4(met.p4() - full_jet_corr);
		// TODO: uncertainties?
		// for leptons they are done as:
		//met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		//met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		//met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		//met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction



		// FIXME: So are these MET corrections?
		if(debug) cout << "Update also MET" << endl;
		// LorentzVector n_met = met.p4();
		// std::vector<LorentzVector> newMet = utils::cmssw::getMETvariations(n_met/*recoMet*/,jets,selLeptons,isMC);
		// FIXME: Must choose a lepton collection. Perhaps loose leptons?
		// n_met = newMet[utils::cmssw::METvariations::NOMINAL];

		//fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), n_met.pt(), weight);
		//if (isSingleMu) fill_pt_e( string("met0_all_leptoncorr_jetcorr_singlemu_pt"), n_met.pt(), weight);
		//if (isSingleE)  fill_pt_e( string("met0_all_leptoncorr_jetcorr_singleel_pt"), n_met.pt(), weight);
		//fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), met.pt(), weight);
		fill_1d(string("control_met_allcorr_pt"), 200, 0., 200., met.pt(), weight);

		if(debug) cout << "Jet Energy Corrections updated" << endl;


		// TODO: should MET corrections be done here?
		// METs with corrections
		// not used now
		//
		//double met_pt_values[7];
		//met_pt_values[0] = n_met.pt();
		//met_pt_values[1] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp).pt();
		//met_pt_values[2] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown).pt();
		//met_pt_values[3] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp).pt();
		//met_pt_values[4] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown).pt();
		//met_pt_values[5] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp).pt();
		//met_pt_values[6] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown).pt();

		// ------------------------------- JETS SELECTION
		// I need different collections because of tau cleaning, but this is needed only for the single lepton channels, so the tau cleaning is performed later.
		pat::JetCollection selJets;
		// TODO: do all jet selection right here
		// now selBJets are not used anywhere
		// selJets pass cross-cleaning with taus later
		// and b-tagging again
		double mindphijmet (9999.);
		for (unsigned int count_ided_jets = 0, ijet = 0; ijet < jets.size(); ++ijet)
			{
			pat::Jet& jet = jets[ijet];

			fill_2d(string("control_jet_jetscorrected_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_jetscorrected_phi"), 128, -3.2, 3.2, jet.phi(), weight);

			// TODO: what do we do here exactly?
			// a loose selection on jets, and then tighten it later?
			// if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue;
			// Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)

			//mc truth for this jet
			//const reco::GenJet * genJet = jet.genJet();
			//TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");
			// TODO: this mctruth for jets it is never used in the code

			//jet id
			bool passPFloose = passPFJetID("Loose", jet); 
			// bool passPFloose = passPFJetID("Tight", jet); 
			//if (label=="Tight")
			// FIXME: check when pileup ID will come out

			if (passPFloose)
				{
				fill_2d(string("control_jet_jetsIDed_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
				fill_1d(string("control_jet_jetsIDed_phi"), 128, -3.2, 3.2, jet.phi(), weight);
				}

			// and now the tighter final selection
			double eta = jet.eta();
			double pt  = jet.pt();
			bool passKino = pt > 30. && fabs(eta) < 2.4;

			// corrections:
			// TODO: are they MC-only?
			/* no smeared  jet pt values for now
			std::vector<double> pt_values;
			if (isMC)
				pt_values = utils::cmssw::smearJES(pt, eta, totalJESUnc);
			else
				{
				pt_values.push_back(pt);
				pt_values.push_back(pt);
				}
			*/
			// vary JesUp   is pt_values[0]
			// vary JesDown is pt_values[1]
			// if (!passPFloose || jet.pt() <30. || fabs(jet.eta()) > 2.5) continue;
			// if (passPFloose && (pt > 30. || pt_values[0] > 30. || pt_values[1] > 30.) && fabs(eta) < 2.5)


			if (passPFloose && passKino)
				{
				selJets.push_back(jet);

				fill_2d(string("control_jet_selJets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
				fill_1d(string("control_jet_selJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);

				// double dphijmet = fabs (deltaPhi (n_met.phi(), jet.phi()));
				double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
				if (dphijmet < mindphijmet) mindphijmet = dphijmet;
				// FIXME: mindphijmet is not used anywhere now
				}
			}

		std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);



		// ---------------------------- Clean jet collections from selected leptons
		// TODO: add gamma-cleaning as well?

		pat::JetCollection selJetsNoLep;
		for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			{
			pat::Jet& jet = selJets[ijet];

			double minDRlj (9999.);

			for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
				minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));

			if (minDRlj < 0.4) continue;

			selJetsNoLep.push_back(jet);

			fill_2d(string("control_jet_selJetsNoLep_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_selJetsNoLep_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}



		// ---------------------------- Clean jet collection from selected taus
		pat::JetCollection selJetsNoLepNoTau;

		for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
			{
			pat::Jet& jet = selJetsNoLep[ijet];

			double minDRtj(9999.);

			for(size_t itau=0; itau < selTausNoLep.size(); ++itau)
				minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTausNoLep[itau]));

			if (minDRtj < 0.4) continue;

			selJetsNoLepNoTau.push_back(jet);

			fill_2d(string("control_jet_selJetsNoLepNoTau_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_selJetsNoLepNoTau_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}


		if(debug){
			cout << "processed jets" << endl;
			}


		// --------------------------- B-TAGGED JETS
		pat::JetCollection selBJets;

		// for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
		for (size_t ijet = 0; ijet < selJetsNoLepNoTau.size(); ++ijet)
			{
			pat::Jet& jet = selJetsNoLep[ijet];

			fill_2d(string("control_jet_bjetCandidates_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_bjetCandidates_phi"), 128, -3.2, 3.2, jet.phi(), weight);

			double eta=jet.eta();

			bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
			//bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.935);
			bool hasCSVtag_BTagUp(false), hasCSVtag_BTagDown(false);

			fill_2d(string("all_b_tagging_candidate_jets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
			if (hasCSVtag) fill_2d(string("all_b_tagging_candidate_jets_pt_eta_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);

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
				// fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

				double sf;
				if (abs(flavId)==5) {
					// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_b_jets_pt_eta_beforesf"), jet.pt(), eta, weight);
					// int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight)
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_b_jets_pt_eta_beforesf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
					// jet_radius(pat::Jet& jet)
					fill_1d(string("btag_radius_flavour_b"), 400, 0., 4.,   jet_radius(jet), weight);

					sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, jet.pt(), 0.);
					// fill_btag_sf(string("btag_sf_flavour_b"), sf, weight);
					// fill_1i(string("weightflow_mu_passmetfilters"), 300, 0, 300,   7, weight);
					fill_1d(string("btag_sf_flavour_b"), 200, 0., 2.,   sf, weight);

					btsfutil.modifyBTagsWithSF(hasCSVtag, sf, beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_b_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_b_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				} else if(abs(flavId)==4) {
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_c_jets_pt_eta_beforesf"), jet.pt(), eta, weight);
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_c_jets_pt_eta_beforesf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
					fill_1d(string("btag_radius_flavour_c"), 400, 0., 4.,   jet_radius(jet), weight);

					// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
					sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, jet.pt(), 0.);
					// fill_btag_sf(string("btag_sf_flavour_c"), sf, weight);
					fill_1d(string("btag_sf_flavour_c"), 200, 0., 2.,   sf, weight);

					btsfutil.modifyBTagsWithSF(hasCSVtag, sf, beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_c_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_c_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				} else {
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_udsg_jets_pt_eta_beforesf"), jet.pt(), eta, weight);
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_udsg_jets_pt_eta_beforesf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
					fill_1d(string("btag_radius_flavour_udsg"), 400, 0., 4.,   jet_radius(jet), weight);

					// hasCSVtag = false; // FIXME: 8.26-27 no b-tagged light jets in MC <<<<<-----------------------------------

					// btsfutil.modifyBTagsWithSF(hasCSVtag, btagCalL.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
					sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, jet.pt(), 0.);
					// fill_btag_sf(string("btag_sf_flavour_udsg"), sf, weight);
					fill_1d(string("btag_sf_flavour_udsg"), 200, 0., 2.,   sf, weight);

					btsfutil.modifyBTagsWithSF(hasCSVtag, sf, leff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
					// if (hasCSVtag) fill_btag_eff(string("mc_all_b_tagged_udsg_jets_pt_eta_aftersf"), jet.pt(), eta, weight);
					if (hasCSVtag) fill_2d(string("mc_all_b_tagged_udsg_jets_pt_eta_aftersf"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
				}
			if (hasCSVtag) fill_2d(string("all_b_tagging_candidate_jets_pt_eta_tagged_after_mc_sfs"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
			}

			if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
				{
				selBJets.push_back(jet);
				fill_2d(string("control_jet_selBJets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
				fill_1d(string("control_jet_selBJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
				}
			}


		if(debug){
			cout << "processed b-tagged jets" << endl;
			}





		// ------------------------------------------ the taus far from jets

		pat::TauCollection selTausNoLepNoJet; // selTausNoLep
		// int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
		for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
			{
			pat::Tau& tau = selTausNoLep[itau];

			// cross-cleaning taus from selected jets
			bool overlapWithJet(false);
			for(int n=0; n<(int)selJetsNoLep.size();++n)
				{
				if (reco::deltaR(tau, selJetsNoLep[n])<0.4)
					{ overlapWithJet=true; break; }
				}
			if (overlapWithJet) continue;

			selTausNoLepNoJet.push_back(tau);
			// so these are the final taus we use in the selection

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			increment(string("weight_of_tausnolepnojet_found"), weight);
			}

		/*
		double selTausNoLep_noAnyJets = 0;
		for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
			{
			pat::Tau& tau = selTausNoLep[itau];

			// cross-cleaning taus with leptons
			bool overlapWithJet(false);
			for(int n=0; n<(int)jets.size();++n)
				{
				if (reco::deltaR(tau, jets[n])<0.4)
					{ overlapWithJet=true; break; }
				}
			if (overlapWithJet) continue;

			selTausNoLep_noAnyJets += weight;
			// selTausNoLepNoJet.push_back(tau);
			// so these are the final taus we use in the selection

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			// increment(string("weight_of_tausnolepnojet_found"), weight);
			}

		double selTausNoLep_noIdJets = 0;
		double selTausNoLep_noSelLowPtJets = 0;

		for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
			{
			pat::Tau& tau = selTausNoLep[itau];

			// cross-cleaning taus with leptons
			bool overlapWithJet(false);
			for(int n=0; n<(int)idJets.size();++n)
				{
				if (reco::deltaR(tau, idJets[n])<0.4)
					{ overlapWithJet=true; break; }
				}
			if (overlapWithJet) continue;

			selTausNoLep_noIdJets += weight;
			// selTausNoLepNoJet.push_back(tau);
			// so these are the final taus we use in the selection

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			// increment(string("weight_of_tausnolepnojet_found"), weight);
			}

		for (size_t itau = 0; itau < selTausNoLep.size(); ++itau)
			{
			pat::Tau& tau = selTausNoLep[itau];

			// cross-cleaning taus with leptons
			bool overlapWithJet(false);
			for(int n=0; n<(int)selLowPtJets.size();++n)
				{
				if (reco::deltaR(tau, selLowPtJets[n])<0.4)
					{ overlapWithJet=true; break; }
				}
			if (overlapWithJet) continue;

			selTausNoLep_noSelLowPtJets += weight;
			// selTausNoLepNoJet.push_back(tau);
			// so these are the final taus we use in the selection

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			// increment(string("weight_of_tausnolepnojet_found"), weight);
			}
		*/


		// -------------------------------------------------- All particles are selected
		// now the channel/selection can be done

		if(debug){
			cout << "all particle-objects are processed, checking channel selection" << endl;
			}

		unsigned int n_leptons = selLeptons.size();

		unsigned int n_taus = selTausNoLep.size();
		unsigned int n_jets = selJetsNoLep.size();
		unsigned int n_bjets = selBJets.size();

		// Select jets for tau-fake-rate study, save the counts

		// W+jets selection:
		bool selection_W_jets = (n_leptons > 0) && (n_jets > 0);

		// And all-jets selection:
		bool selection_QCD_jets = (n_jets > 1);


		pat::JetCollection probeJets;
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

		// for now just take all jets:
		probeJets = selJetsNoLep;

		//TH3F jets_distr;
		//TH3F tau_jets_distr;
		float tau_fake_distance = 0.1; // the distance to tau for a jet to be considered tau's origin

		// loop through probe jets and fill jets_distr
		// loop through jets, find the ones close to taus -- fill tau_jets_distr
		// new standard function for jet distr
		//int fill_3d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, double x, double y, double z, double weight)

		// Check among 2 channels of selection:
		//     QCD-channel -- 2 or more jets
		//     W+jets      -- only 1 iso muon, one or more jets

		bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;

		bool isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && clean_lep_conditions;
		//bool isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && clean_lep_conditions;

		bool Wjets_selection = clean_lep_conditions && isSingleMu && (selJetsNoLep.size() > 0);
		bool QCD_selection  = selJetsNoLep.size() > 1;

		if (Wjets_selection)
			{

			// WJets, HLT channels
			if (JetHLTTrigger && MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   10, weight);
				}
			else if (JetHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   11, weight);
				}
			else if (MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   12, weight);
				}
			else continue;

			fill_1d(string(hlt_channel + "wjets_selection_ntaus"), 20, 0, 20,   jets.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTaus"), 20, 0, 20,   selJets.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTausNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselTausNoLepNoJet"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

			fill_1d(string(hlt_channel + "wjets_selection_njets"), 20, 0, 20,   taus.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJets"), 20, 0, 20,   selTaus.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJetsNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
			fill_1d(string(hlt_channel + "wjets_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

			/*
			// N jets, taus
			increment( hlt_channel + string("wjets_selection_njets"), weight*jets.size() );
			// pat::JetCollection idJets;
			// pat::JetCollection selLowPtJets;
			increment( hlt_channel + string("wjets_selection_nidJets"), weight*idJets.size() );
			increment( hlt_channel + string("wjets_selection_nselLowPtJets"), weight*selLowPtJets.size() );
			increment( hlt_channel + string("wjets_selection_nselJetsNoLep"), weight*selJetsNoLep.size() );
			// selJetsNoLepNoTau
			increment( hlt_channel + string("wjets_selection_nselJetsNoLepNoTau"), weight*selJetsNoLepNoTau.size() );

			increment( hlt_channel + string("wjets_selection_ntaus"), weight*taus.size() );
			increment( hlt_channel + string("wjets_selection_nselTausNoLep"), weight*selTausNoLep.size() );
			// double selTausNoLep_noIdJets = 0;
			// double selTausNoLep_noSelLowPtJets = 0;
			increment( hlt_channel + string("wjets_selection_nselTausNoLep_noIdJets"), selTausNoLep_noIdJets );
			increment( hlt_channel + string("wjets_selection_nselTausNoLep_noSelLowPtJets"), selTausNoLep_noSelLowPtJets );
			increment( hlt_channel + string("wjets_selection_nselTausNoLepNoJet"), weight*selTausNoLepNoJet.size() );
			*/

			// TODO: make only 1 loop through jets,
			// checking closeness to taus along the way

			// loop through all jets,
			// save them
			// if it is MC save also the gluon/quark composition
			for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
				{
				pat::Jet& jet = selJetsNoLep[ijet];
				// wjets_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
				//fill_jet_distr(string control_point_name, double weight, double pt, double eta, double radius)
				fill_jet_distr(hlt_channel + string("wjets_jets_distr"), weight, jet.pt(), jet.eta(), jet_radius(jet));

				// const reco::GenParticle* genParton()
				if (isMC)
					{
					//const reco::GenParticle* jet_origin = jet.genParton();
					// TODO this genParton is apparently not defined in PAT jets
					// using the jet.partonFlavor
					// need to check if it is correct parton ID from generation MC
					// there is also hadronFlavor..
					// the ID should be in:
					// jet_origin->pdgId();
					//wjets_jet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(hlt_channel + string("wjets_jet_origins"), 100, 0, 100,   abs(jet.partonFlavour()), weight);
					//wjets_jet_origin->Fill(abs( jet_origin->pdgId() ));
					// wjets_taujet_origin
					}

				// loop through taus
				// find jets close to tau -- consider them faking the tau
				// save such jets
				// if it is MC save also their gluon/quark origin -- just in case, for interest
				for(size_t itau=0; itau < selTausNoLep.size(); ++itau)
					{
					//for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
					//{
					//double fake_distance = reco::deltaR(selJetsNoLep[ijet], selTausNoLep[itau]);
					double fake_distance = reco::deltaR(jet, selTausNoLep[itau]);
					//wjets_taujet_distance->Fill(fake_distance);
					fill_1d(hlt_channel + string("wjets_taujet_distance"), 100, 0, 2,   fake_distance, weight);

					if (fake_distance <= tau_fake_distance)
						{
						// the tau is fake by this jet -- save distr
						// wjets_tau_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
						fill_jet_distr(hlt_channel + string("wjets_tau_jets_distr"), weight, jet.pt(), jet.eta(), jet_radius(jet));
						//fill_pt_pt(hlt_channel + string("wjets_tau_taujet_pts"), jet.pt(), selTausNoLep[itau].pt(), weight);
						fill_2d(hlt_channel + string("wjets_tau_taujet_pts"), 400, 0., 400., 400, 0., 400, jet.pt(), selTausNoLep[itau].pt(), weight);

						// N tau-jets
						//increment( hlt_channel + string("wjets_selection_ntaujet"), weight );

						// const reco::GenParticle* genParton()
						if (isMC)
							{
							//const reco::GenParticle* jet_origin = selJetsNoLep[ijet].genParton();
							// the ID should be in:
							// jet_origin->pdgId();
							//wjets_taujet_origin->Fill(abs( jet.partonFlavour() ));
							fill_1d(hlt_channel + string("wjets_taujet_origins"), 100, 0, 100,   abs(jet.partonFlavour()), weight);
							//wjets_taujet_origin->Fill(abs( jet_origin->pdgId() ));
							// wjets_jet_origin
							}
						continue;
						}
					//}
					}
				}

				// and the fakerate of "distinct taus"
				/* TODO: apparently pat:Tau-s may not have phiphi/etaeta momenta for the radius
				for(size_t itau=0; itau < selTausNoLepNoJet.size(); ++itau)
					{
					//
					pat::Tau& tau = selTausNoLepNoJet[itau];
					// wjets_distinct_tau_distr->Fill(tau.pt(), tau.eta(), jet_radius(tau));
					fill_jet_distr(hlt_channel + string("wjets_distinct_tau_distr"), weight, tau.pt(), tau.eta(), jet_radius(tau));

					// there is no partonFlavour origin for taus..........
					}
				*/
			}

		if (QCD_selection)
			{

			// QCD, HLT channels
			if (JetHLTTrigger && MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   20, weight);
				}
			else if (JetHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   21, weight);
				}
			else if (MuonHLTTrigger)
				{
				fill_1d(string("weightflow"), 300, 0, 300,   22, weight);
				}
			else continue;

			fill_1d(string(hlt_channel + "qcd_selection_ntaus"), 20, 0, 20,   jets.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTaus"), 20, 0, 20,   selJets.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTausNoLep"), 20, 0, 20,   selJetsNoLep.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselTausNoLepNoJet"), 20, 0, 20,   selJetsNoLepNoTau.size(), weight);

			fill_1d(string(hlt_channel + "qcd_selection_njets"), 20, 0, 20,   taus.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJets"), 20, 0, 20,   selTaus.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJetsNoLep"), 20, 0, 20,   selTausNoLep.size(), weight);
			fill_1d(string(hlt_channel + "qcd_selection_nselJetsNoLepNoTau"), 20, 0, 20,   selTausNoLepNoJet.size(), weight);

			/*
			increment( hlt_channel + string("weightflow_qcd_selection"), weight );
			// N jets, taus
			increment( hlt_channel + string("qcd_selection_njets"), weight*jets.size() );
			increment( hlt_channel + string("qcd_selection_nidJets"), weight*idJets.size() );
			increment( hlt_channel + string("qcd_selection_nselLowPtJets"), weight*selLowPtJets.size() );
			increment( hlt_channel + string("qcd_selection_nselJetsNoLep"), weight*selJetsNoLep.size() );
			increment( hlt_channel + string("qcd_selection_nselJetsNoLepNoTau"), weight*selJetsNoLepNoTau.size() );

			increment( hlt_channel + string("qcd_selection_ntaus"), weight*taus.size() );
			increment( hlt_channel + string("qcd_selection_nselTausNoLep"), weight*selTausNoLep.size() );
			//
			increment( hlt_channel + string("qcd_selection_nselTausNoLep_noIdJets"), selTausNoLep_noIdJets );
			increment( hlt_channel + string("qcd_selection_nselTausNoLep_noSelLowPtJets"), selTausNoLep_noSelLowPtJets );
			increment( hlt_channel + string("qcd_selection_nselTausNoLepNoJet"), weight*selTausNoLepNoJet.size() );
			*/


			for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
				{
				pat::Jet& jet = selJetsNoLep[ijet];
				// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
				// selection jets
				fill_jet_distr(hlt_channel + string("qcd_jets_distr"), weight, jet.pt(), jet.eta(), jet_radius(jet));
				//fill_3d(hlt_channel + string("qcd_jets_distr"), 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad, 300, 0, 300,   20, weight);

				// const reco::GenParticle* genParton()
				// jet parton origin
				if (isMC)
					{
					//const reco::GenParticle* jet_origin = jet.genParton();
					// the ID should be in:
					// jet_origin->pdgId();
					//qcd_jet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(hlt_channel + string("qcd_jet_origins"), 100, 0, 100,   abs(jet.partonFlavour()), weight);
					//qcd_jet_origin->Fill(abs( jet_origin->pdgId() ));
					// qcd_taujet_origin
					}

				for(size_t itau=0; itau < selTausNoLep.size(); ++itau)
					{
					// selection taus (fake taus)
					//for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
					//{
					//double fake_distance = reco::deltaR(selJetsNoLep[ijet], selTausNoLep[itau]);
					double fake_distance = reco::deltaR(jet, selTausNoLep[itau]);
					//qcd_taujet_distance->Fill(fake_distance);
					fill_1d(hlt_channel + string("qcd_taujet_distance"), 100, 0, 2,   fake_distance, weight);

					if (fake_distance <= tau_fake_distance)
						{
						// the tau is fake by this jet -- save distr
						//qcd_tau_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
						fill_jet_distr(hlt_channel + string("qcd_tau_jets_distr"), weight, jet.pt(), jet.eta(), jet_radius(jet));
						// fill_pt_pt(string control_point_name, double pt1, double pt2, double weight)
						//fill_pt_pt(hlt_channel + string("qcd_tau_taujet_pts"), jet.pt(), selTausNoLep[itau].pt(), weight);
						fill_2d(hlt_channel + string("qcd_tau_taujet_pts"), 400, 0., 400., 400, 0., 400, jet.pt(), selTausNoLep[itau].pt(), weight);

						// N tau-jets
						//increment( hlt_channel + string("qcd_selection_ntaujet"), weight );

						// const reco::GenParticle* genParton()
						// jet parton origin for faking jet (just in case)
						if (isMC)
							{
							//const reco::GenParticle* jet_origin = selJetsNoLep[ijet].genParton();
							// the ID should be in:
							// jet_origin->pdgId();
							//qcd_taujet_origin->Fill(abs( jet.partonFlavour() ));
							fill_1d(hlt_channel + string("qcd_taujet_origins"), 100, 0, 100,   abs(jet.partonFlavour()), weight);
							//qcd_taujet_origin->Fill(abs( jet_origin->pdgId() ));
							// qcd_jet_origin
							}
						continue;
						}
					//}
					}
				}

				// and the fakerate of "distinct taus"
				/* TODO: apparently pat:Tau-s may not have phiphi/etaeta momenta for the radius
				for(size_t itau=0; itau < selTausNoLepNoJet.size(); ++itau)
					{
					//
					pat::Tau& tau = selTausNoLepNoJet[itau];
					// qcd_distinct_tau_distr->Fill(tau.pt(), tau.eta(), jet_radius(tau));
					fill_jet_distr(hlt_channel + string("qcd_distinct_tau_distr"), weight, tau.pt(), tau.eta(), jet_radius(tau));
					// there is no partonFlavour origin for taus..........
					}
				*/
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


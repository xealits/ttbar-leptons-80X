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
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/ttbar-leptons-80X/interface/ProcessingGenParticles.h"

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

#include "UserCode/ttbar-leptons-80X/interface/recordFuncs.h"

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
			//eta=fabs(eta);
			//double ptSF(1.0), ptSF_err(0.06);
			//if      (eta<0.5) { ptSF=1.095; ptSF_err=0.018; }
			//else if (eta<0.8) { ptSF=1.120; ptSF_err=0.028; }
			//else if (eta<1.1) { ptSF=1.097; ptSF_err=0.017; }
			//else if (eta<1.3) { ptSF=1.103; ptSF_err=0.033; }
			//else if (eta<1.7) { ptSF=1.118; ptSF_err=0.014; }
			//else if (eta<1.9) { ptSF=1.100; ptSF_err=0.033; }
			//else if (eta<2.1) { ptSF=1.162; ptSF_err=0.044; }
			//else if (eta<2.3) { ptSF=1.160; ptSF_err=0.048; }
			//else if (eta<2.5) { ptSF=1.161; ptSF_err=0.060; }
			//else if (eta<2.8) { ptSF=1.209; ptSF_err=0.059; }
			//else if (eta<3.0) { ptSF=1.564; ptSF_err=0.321; }
			//else if (eta<3.2) { ptSF=1.384; ptSF_err=0.033; }
			//else if (eta<5.0) { ptSF=1.216; ptSF_err=0.050; }

			// 2016 80X
			// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
			//abs(eta) region 	0.0–0.5 	0.5-0.8 	0.8–1.1 	1.1-1.3 	1.3–1.7 	1.7 - 1.9 	1.9–2.1 	2.1 - 2.3 	2.3 - 2.5 	2.5–2.8 	2.8-3.0 	3.0-3.2 	3.2-5.0
			//Data/MC SF
			//1.122 +-0.026
			//1.167 +-0.048
			//1.168 +-0.046
			//1.029 +-0.066
			//1.115 +-0.03
			//1.041 +-0.062
			//1.167 +-0.086
			//1.094 +-0.093
			//1.168 +-0.120
			//1.266 +-0.132
			//1.595 +-0.175
			//0.998 +-0.066
			//1.226 +-0.145
			eta=fabs(eta);
			double ptSF(1.0), ptSF_err(0.06);
			if      (eta<0.5) { ptSF=1.122; ptSF_err=0.026; }
			else if (eta<0.8) { ptSF=1.167; ptSF_err=0.048; }
			else if (eta<1.1) { ptSF=1.168; ptSF_err=0.046; }
			else if (eta<1.3) { ptSF=1.029; ptSF_err=0.066; }
			else if (eta<1.7) { ptSF=1.115; ptSF_err=0.03 ; }
			else if (eta<1.9) { ptSF=1.041; ptSF_err=0.062; }
			else if (eta<2.1) { ptSF=1.167; ptSF_err=0.086; }
			else if (eta<2.3) { ptSF=1.094; ptSF_err=0.093; }
			else if (eta<2.5) { ptSF=1.168; ptSF_err=0.120; }
			else if (eta<2.8) { ptSF=1.266; ptSF_err=0.132; }
			else if (eta<3.0) { ptSF=1.595; ptSF_err=0.175; }
			else if (eta<3.2) { ptSF=0.998; ptSF_err=0.066; }
			else if (eta<5.0) { ptSF=1.226; ptSF_err=0.145; }

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

	// at time of 80X all 13TeV (74X, 76X, 80X) recommendations got new specification for |eta| < 2.7
	if (label=="Loose")
		{
		// passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
		// the same, just with new names:
		// passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 
		if (fabs(eta) <= 3.0)
			{
			if (fabs(eta) <= 2.7) passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4);
			else passID = NEMF<0.90 && NumNeutralParticles > 2;
			}
		else passID = NEMF<0.90 && NumNeutralParticles > 10;
		}
	if (label=="Tight")
		{
		// passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
		// passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;
		// the same, only NHF and NEMF are 0.90 at |eta| < 2.7
		if (fabs(eta) <= 3.0)
			{
			if (fabs(eta) <= 2.7) passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4);
			else passID = NEMF<0.90 && NumNeutralParticles > 2;
			}
		else passID = NEMF<0.90 && NumNeutralParticles > 10;
		}

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

string find_W_decay_old(const reco::Candidate * W) {
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
double jetToTauFakeRate(TH3F * tau_fake_rate_jets_histo1, TH3F * tau_fake_rate_taus_histo1, TH3F * tau_fake_rate_jets_histo2, TH3F * tau_fake_rate_taus_histo2, Double_t tau_fake_rate_histo1_fraction, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
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
	// -- taking all combined for 8 TeV
	double a = 0.156;
	double b = -0.00137;
	return TMath::Exp(a + b*x);
	}




//string mc_decay("");
// Some MC datasets are inclusive, but analysis needs info on separate channels from them
// thus the events are traversed and the channel is found
// currently it is done only for TTbar channel (isTTbarMC)
//
// the sub-channel of MC is paired together with the distr_name string into the key of control distrs
// it is then printed out with dtag

// TODO: rearrange the code into particle selection and channel selection
// TODO organize the code with new general record functions and remove these


//void printDaughters(const reco::GenParticle & p, int offset)
void printDaughters(const reco::Candidate * p, int offset)
{
int n_daughters = p->numberOfDaughters();
if (n_daughters>0)
	{
	for (int j = 0; j < n_daughters; ++j) {
		const reco::Candidate * d = p->daughter( j );
		//const reco::GenParticle * d = p.daughter( j );
		for (int o=0; o<offset; o++) cout << ' ';
		cout << d->pdgId() << " (" << d->status() << ")" << " <- " << abs(p->pdgId()) << endl;
		printDaughters(d, offset+1);
		}
	}
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


// Data-driven tau fakerate background

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



//MC normalization (to 1/pb)
if(debug) cout << "DEBUG: xsec: " << xsec << endl;

// ------------------------------------- jet energy scale and uncertainties 
TString jecDir = runProcess.getParameter < std::string > ("jecDir");
gSystem->ExpandPathName (jecDir);
// v1
// getJetCorrector(TString baseDir, TString pf, bool isMC)
//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
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
beff = 0.747;

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
// BTagCalibrationReader btagCal(BTagEntry::OP_TIGHT,  // operating point
                             "central",             // central sys type
                             {"up", "down"});      // other sys types
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_B,      // btag flavour
            "comb");               // measurement type
//          "mujets");               // measurement type
btagCal.load(btagCalib,              // calibration instance
            BTagEntry::FLAV_C,      // btag flavour
            "comb");              // measurement type
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


if(debug)
	{
	cout << "Some input parameters\n";

	cout << "isMC = " << isMC << "\n";
	cout << "isW0JetsSet = " << isW0JetsSet << "\n";

	cout << "isTTbarMC = "    << isTTbarMC << "\n";
	cout << "isPromptReco = " << isPromptReco << "\n";
	cout << "isRun2015B = "   << isRun2015B << "\n";
	cout << "isNLOMC = "      << isNLOMC << "\n";

	cout << "jecDir = "      << jecDir << "\n";
	}











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

	int treeStep (ev.size()/50);

	for (ev.toBegin(); !ev.atEnd(); ++ev)
		{
		iev++;

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

		// tt b and W products
		//
        	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> t_b_p4_s, tb_b_p4_s, t_W_p4_s, tb_W_p4_s;
        	vector<Int_t> t_b_pdgId_s, tb_b_pdgId_s, t_W_pdgId_s, tb_W_pdgId_s;
        	vector<Int_t> t_b_status_s, tb_b_status_s, t_W_status_s, tb_W_status_s;
        	vector<const reco::Candidate*> t_b_parts, tb_b_parts, t_W_parts, tb_W_parts;


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
						const reco::Candidate * b = p.daughter( 1-W_num );
						save_final_states(W, t_W_p4_s, t_W_pdgId_s, t_W_status_s, t_W_parts);
						save_final_states(b, t_b_p4_s, t_b_pdgId_s, t_b_status_s, t_b_parts);
						mc_decay += find_W_decay_old(W);
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
						const reco::Candidate * b = p.daughter( 1-W_num );
						save_final_states(W, tb_W_p4_s, tb_W_pdgId_s, tb_W_status_s, tb_W_parts);
						save_final_states(b, tb_b_p4_s, tb_b_pdgId_s, tb_b_status_s, tb_b_parts);
						mc_decay += find_W_decay_old(W);
						mc_decay += string("bar");
					}
				}
			}
			// so, mc_decay will be populated with strings matching t decays
			// hopefully, in TTbar sample only ttbar decays are present and it is not ambigous
		}

		if (debug) {
			cout << "MC suffix " << mc_decay << " is found\n";

			// printing out t tb final states:
			cout << "t W final states:" << endl;
			for (int i=0; i< t_W_parts.size(); i++)
				cout << "  " << t_W_pdgId_s[i] << '(' << t_W_status_s[i] << ')';
			cout << endl;
			cout << "t b final states:" << endl;
			for (int i=0; i< t_b_parts.size(); i++)
				cout << "  " << t_b_pdgId_s[i] << '(' << t_b_status_s[i] << ')';
			cout << endl;
			cout << "tb W final states:" << endl;
			for (int i=0; i< t_W_parts.size(); i++)
				cout << "  " << t_W_pdgId_s[i] << '(' << t_W_status_s[i] << ')';
			cout << endl;
			cout << "tb b final states:" << endl;
			for (int i=0; i< tb_b_parts.size(); i++)
				cout << "  " << tb_b_pdgId_s[i] << '(' << tb_b_status_s[i] << ')';
			cout << endl;
			}

		//if (!mc_decay.empty()) mc_decay = string("_") + mc_decay;
		//mc_decay = string("_") + mc_decay; // so we'll have "_" or "_mcdecay"

		bool only_tt_lepjets = false;
		if (only_tt_lepjets)
			{
			if (mc_decay != string("elqbar") && mc_decay != string("qelbar") && mc_decay != string("qmubar") && mc_decay != string("muqbar")) continue;
			}

		// recursive print of genTree
		// cmssw crap doesn't print and notes from 2013 are far outdated (include instead of load?)
		// how to print it here?
		// which particles are "original" mothers
		// let's just take two first ones...

		const reco::GenParticle & proton_1 = gen[0];
		const reco::GenParticle & proton_2 = gen[1];

		cout << proton_1.pdgId() << " (" << proton_1.status() << ')' << endl;
		//printDaughters(proton_1, 1);
		int offset = 1;
		int n_daughters = proton_1.numberOfDaughters();
		if (n_daughters>0)
			{
			for (int j = 0; j < n_daughters; ++j) {
				const reco::Candidate * d = proton_1.daughter( j );
				//const reco::GenParticle * d = p.daughter( j );
				for (int o=0; o<offset; o++) cout << ' ';
				cout << d->pdgId() << " (" << d->status() << ")" << endl;
				printDaughters(d, offset+1);
				}
			}

		cout << proton_2.pdgId() << " (" << proton_2.status() << ')' << endl;
		//printDaughters(proton_2, 1);
		n_daughters = proton_2.numberOfDaughters();
		if (n_daughters>0)
			{
			for (int j = 0; j < n_daughters; ++j) {
				const reco::Candidate * d = proton_2.daughter( j );
				//const reco::GenParticle * d = p.daughter( j );
				for (int o=0; o<offset; o++) cout << ' ';
				cout << d->pdgId() << " (" << d->status() << ")" << endl;
				printDaughters(d, offset+1);
				}
			}

		//* List of mother-daughters for all particles
		//* TODO: make it into a separate function

		/*
		 * Prints the nodes of decay tree:
		 * 0 each GenParticle is printed (pdgId & status)
		 * 1 with its' (possible) mothers (pdgId & status)
		 * 2 and (possible) daughters
		 *
		 * N_part: [ mom1_id status;[ mom2...]] >>| particle_string id status |>> [ daught1_id status;[ daught2...]]
		 * 
		 */

		vector<LorentzVector> vis_taus;
		for(size_t i = 0; i < gen.size(); ++ i) {
			const reco::GenParticle & p = gen[i];
			int id = p.pdgId();
			int st = p.status();
			int n_daughters = p.numberOfDaughters();
			//cout << i << ": " ;
			for (int j = 0 ; j < p.numberOfMothers(); ++j) {
				const reco::Candidate * mom = p.mother(j);
				//cout << " " << mom->pdgId() << " " << mom->status() << ";";
				}
			//if (p.numberOfMothers() != 0) cout << "\t >>| \t" ;
			//cout <<  id << " " << st << "\t" << p.pt() << "," << p.eta() << "," << p.phi() << "," << p.mass();
			if (n_daughters>0) {
				//cout << " |>> \t" ;
				for (int j = 0; j < n_daughters; ++j) {
					const reco::Candidate * d = p.daughter( j );
					//cout << d->pdgId() << " " << d->status() << "; " ;
					}
				}
			//cout << endl;

			// if it is a final state tau
			//  the status is 1 or 2
			//  1. (final state, not decays)
			//  2. (decayed or fragmented -- the case for tau)
			if (fabs(id) == 15 && (st == 1 || st == 2))
				{
				//cout << "A final state tau!" << endl;
				//cout << p << endl;
				//cout << "visible daughters:" << endl;
				LorentzVector vis_ds(0,0,0,0);
				for (int j = 0; j < n_daughters; ++j) {
					const reco::Candidate * d = p.daughter(j);
					unsigned int d_id = fabs(d->pdgId());
					if (d_id == 12 || d_id == 14 || d_id == 16) continue;
					//cout << d->pdgId() << " (" << d->status() << ", " << d->numberOfDaughters() << "); " ;
					vis_ds += d->p4();
					}
				vis_taus.push_back(vis_ds);
				//cout << endl << "sum = " << vis_ds << endl;
				//cout << "sum = " << vis_ds.pt() << " " << vis_ds.eta() << endl;
				}
			}


		// JETS, their generation hadron/parton IDs

		pat::JetCollection jets;
		fwlite::Handle<pat::JetCollection>jetsHandle;
		jetsHandle.getByLabel(ev, "slimmedJets");
		if(jetsHandle.isValid() ) jets = *jetsHandle;
		else
			{
			cout << "ooops" << endl;
			continue;
			}

		cout << "partonFlavour, hadronFlavour, dR-match-to-taus of slimmedJets" << endl;
		for (unsigned int count_ided_jets = 0, ijet = 0; ijet < jets.size(); ++ijet)
			{
			pat::Jet& jet = jets[ijet];
			cout << jet.partonFlavour() << "\t" << jet.hadronFlavour() << "\t";

			double minDRtj (9999.);
			unsigned int closest_tau = -1;
			for (size_t i = 0; i < vis_taus.size(); i++)
				{
				double jet_tau_distance = TMath::Min(minDRtj, reco::deltaR (jet, vis_taus[i]));
				if (jet_tau_distance<minDRtj)
					{
					closest_tau = i;
					minDRtj = jet_tau_distance;
					}
				}
			cout << minDRtj << endl;
			}

		/*
		pat::JetCollection selJets;
		for (unsigned int count_ided_jets = 0, ijet = 0; ijet < jets.size(); ++ijet)
			{
			pat::Jet& jet = jets[ijet];

			fill_2d(string("control_jet_jetscorrected_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_jetscorrected_phi"), 128, -3.2, 3.2, jet.phi(), weight);

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

			// Jet Kinematics
			double eta = jet.eta();
			double pt  = jet.pt();
			//bool passKino = pt > 30. && fabs(eta) < 2.4;


			//if (passPFloose && passKino)
			if (passPFloose && (fabs(eta) < 2.4))
				{
				if (pt > 30.)
					{
					// trying new fake rate application:
					// if there are only 2 30GeV jets -- they are not considered for faking taus
					// TODO: and what to do with 20GeV jets? in selection and in fake rates?
					//  --- for now I will just add 20-30 GeV jets for fake rates..
					//      one should actually try without them
					selJets30GeV.push_back(jet);
					selJets.push_back(jet);

					fill_2d(string("control_jet_selJets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
					fill_1d(string("control_jet_selJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
					// double dphijmet = fabs (deltaPhi (n_met.phi(), jet.phi()));
					//double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
					//if (dphijmet < mindphijmet) mindphijmet = dphijmet;
					// FIXME: mindphijmet is not used anywhere now
					}
				else if (pt > 20.)
					{
					selJets20GeV.push_back(jet);
					}
				}
			}
		*/

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

if (outTxtFile) fclose (outTxtFile);

// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}


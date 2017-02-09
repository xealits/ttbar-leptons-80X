
#include "CommonTools/Utils/interface/PtComparator.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
//#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>


#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>

#include <memory>
#include <random>

using namespace JME;


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
	//
	//Moriond:
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF	1.109 +-0.008	1.138 +-0.013	1.114 +-0.013	1.123 +-0.024	1.084 +-0.011	1.082 +-0.035	1.140 +-0.047	1.067 +-0.053	1.177 +-0.041	1.364 +-0.039	1.857 +-0.071	1.328 +-0.022	1.16 +-0.029
	//
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF
	// 1.109 +-0.008
	// 1.138 +-0.013
	// 1.114 +-0.013
	// 1.123 +-0.024
	// 1.084 +-0.011
	// 1.082 +-0.035
	// 1.140 +-0.047
	// 1.067 +-0.053
	// 1.177 +-0.041
	// 1.364 +-0.039
	// 1.857 +-0.071
	// 1.328 +-0.022
	// 1.16  +-0.029
	eta=fabs(eta);
	double ptSF(1.0), ptSF_err(0.06);
	if      (eta<0.5) { ptSF=1.109; ptSF_err=0.008; }
	else if (eta<0.8) { ptSF=1.138; ptSF_err=0.013; }
	else if (eta<1.1) { ptSF=1.114; ptSF_err=0.013; }
	else if (eta<1.3) { ptSF=1.123; ptSF_err=0.024; }
	else if (eta<1.7) { ptSF=1.084; ptSF_err=0.011; }
	else if (eta<1.9) { ptSF=1.082; ptSF_err=0.035; }
	else if (eta<2.1) { ptSF=1.140; ptSF_err=0.047; }
	else if (eta<2.3) { ptSF=1.067; ptSF_err=0.053; }
	else if (eta<2.5) { ptSF=1.177; ptSF_err=0.041; }
	else if (eta<2.8) { ptSF=1.364; ptSF_err=0.039; }
	else if (eta<3.0) { ptSF=1.857; ptSF_err=0.071; }
	else if (eta<3.2) { ptSF=1.328; ptSF_err=0.022; }
	else if (eta<5.0) { ptSF=1.16 ; ptSF_err=0.029; }

	/*
	toReturn[0]=TMath::Max(0., (genPt+ptSF*(pt-genPt))/pt );
	toReturn[1]=TMath::Max(0., (genPt+(ptSF+ptSF_err)*(pt-genPt))/pt );
	toReturn[2]=TMath::Max(0., (genPt+(ptSF-ptSF_err)*(pt-genPt))/pt );
	*/

	// TODO: check new SF-scaling application:
	// from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
	// (called c_JER there)
	toReturn[0]=TMath::Max(0., (1 + (ptSF          - 1)*(pt - genPt)/pt) );
	toReturn[1]=TMath::Max(0., (1 + (ptSF+ptSF_err - 1)*(pt - genPt)/pt) );
	toReturn[2]=TMath::Max(0., (1 + (ptSF-ptSF_err - 1)*(pt - genPt)/pt) );

	return toReturn;
	}

std::vector<double> JER_SF(double eta)
	{
	eta=fabs(eta);
	double ptSF(1.0), ptSF_err(0.06);
	std::vector<double> toReturn(2, ptSF);

	//Moriond:
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF	1.109 +-0.008	1.138 +-0.013	1.114 +-0.013	1.123 +-0.024	1.084 +-0.011	1.082 +-0.035	1.140 +-0.047	1.067 +-0.053	1.177 +-0.041	1.364 +-0.039	1.857 +-0.071	1.328 +-0.022	1.16 +-0.029
	//
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF
	// 1.109 +-0.008
	// 1.138 +-0.013
	// 1.114 +-0.013
	// 1.123 +-0.024
	// 1.084 +-0.011
	// 1.082 +-0.035
	// 1.140 +-0.047
	// 1.067 +-0.053
	// 1.177 +-0.041
	// 1.364 +-0.039
	// 1.857 +-0.071
	// 1.328 +-0.022
	// 1.16  +-0.029
	if      (eta<0.5) { ptSF=1.109; ptSF_err=0.008; }
	else if (eta<0.8) { ptSF=1.138; ptSF_err=0.013; }
	else if (eta<1.1) { ptSF=1.114; ptSF_err=0.013; }
	else if (eta<1.3) { ptSF=1.123; ptSF_err=0.024; }
	else if (eta<1.7) { ptSF=1.084; ptSF_err=0.011; }
	else if (eta<1.9) { ptSF=1.082; ptSF_err=0.035; }
	else if (eta<2.1) { ptSF=1.140; ptSF_err=0.047; }
	else if (eta<2.3) { ptSF=1.067; ptSF_err=0.053; }
	else if (eta<2.5) { ptSF=1.177; ptSF_err=0.041; }
	else if (eta<2.8) { ptSF=1.364; ptSF_err=0.039; }
	else if (eta<3.0) { ptSF=1.857; ptSF_err=0.071; }
	else if (eta<3.2) { ptSF=1.328; ptSF_err=0.022; }
	else if (eta<5.0) { ptSF=1.16 ; ptSF_err=0.029; }

	toReturn[0]=TMath::Max(0., ptSF);
	toReturn[1]=TMath::Max(0., ptSF_err);

	return toReturn;
	}




/* TODO: taken from the latest MacroUtils of llvv
 * need to check it
 */
std::vector<float> smearJES(float pt, float eta, JetCorrectionUncertainty *jecUnc)
	{
	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float relShift=fabs(jecUnc->getUncertainty(true));
	std::vector<float> toRet;
	toRet.push_back((1.0+relShift)*pt);
	toRet.push_back((1.0-relShift)*pt);
	return toRet;
	}



bool passPFJetID(std::string label, pat::Jet jet)
	{
	// Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

	bool looseJetID = false;
	bool tightJetID = false;
	bool passID(false); 

	// float rawJetEn(jet.correctedJet("Uncorrected").energy() );

	double eta = fabs(jet.eta());
 
 	// from the twiki:
 	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	float NHF = jet.neutralHadronEnergyFraction();
	float NEMF = jet.neutralEmEnergyFraction();
	float CHF = jet.chargedHadronEnergyFraction();
	float MUF = jet.muonEnergyFraction();
	float CEMF = jet.chargedEmEnergyFraction(); // chargedEmEnergyFraction (relative to uncorrected jet energy)
	float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	float NumNeutralParticles =jet.neutralMultiplicity();
	float CHM = jet.chargedMultiplicity(); 
	// TODO: check if these change corresponding to jet corrections? Or apply corrections after passing the selection?
	//
	// the doc https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_14/doc/html/d6/d00/classpat_1_1Jet.html
	// says these `...Fraction` are respective uncorrected jets
	// which is exactly what is needed by their note:
	// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.

	// float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
	// float nef( jet.neutralEmEnergy()/rawJetEn );
	// float cef( jet.chargedEmEnergy()/rawJetEn );
	// float chf( jet.chargedHadronEnergy()/rawJetEn );
	// float nch    = jet.chargedMultiplicity();
	// float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	// float muf(jet.muonEnergy()/rawJetEn); 

	// at time of 80X all 13TeV (74X, 76X, 80X) recommendations got new specification for |eta| < 2.7
	/*
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
	*/

	// and latest (at Moriond17) stuff:
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	// quote:
	// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.
	// --- khmm....
	// in MINIAOD the jets are already corrected  --> one gets the uncorrected jet and reapplies the corrections
	// > ... collection of AK4 jets slimmedJets, made from ak4PFJetsCHS ... "These have standard jet energy corrections applied (L1FastJet, L2, L3), and a pT cut at 10 GeV"
	// so now one need to get the uncorrected jet --> 
	if (eta <= 2.7)
		{
		looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ;
		tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ;
		}
	else if (eta <= 3.0)
		{
		looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
		tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
		}
	else
		{
		looseJetID = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
		tightJetID = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
		}

	// there is also:
	//tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=3.0

	// And the abs(eta)>3.0 part, but we never consider such jets, so... Meh!

	if (label == "Loose")
		return looseJetID;
	if (label == "Tight")
		return tightJetID;

	return passID; 
	}




int processJets_CorrectJES_SmearJERnJES_ID_ISO_Kinematics(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
	bool isMC, double weight,
	double rho, unsigned int nGoodPV,
	FactorizedJetCorrector *jesCor,
	JetCorrectionUncertainty *totalJESUnc,
	double dR_max, // for jet matching in jet corrections smearing for MC
	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
	string& jetID,
	double pt_cut, double eta_cut,
	TRandom *r3,   // the randomizer for the smearing
	LorentzVector& full_jet_corr, pat::JetCollection& selJets,                          // output
	bool record, bool debug) // more output

{
// the PF ID (in Moriond17 recommended to be applied before corrections)
/*
JME::JetResolution resolution = JME::JetResolution(jetResolutionFileName);
JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor(jetResolutionSFFileName);
Variation m_systematic_variation = Variation::NOMINAL; // FIXME: it should be in some headers, included before... but remake it somehow
*/


pat::JetCollection IDjets;

for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	// the input jets here are, supposedly, slimmedJets from MINIAODs -- which are already corrected
	// -- apparently, the parameters that are used in PF ID calculation use the uncorrected values stored in the jets
	// so, everything is fine anyway
	pat::Jet& jet = jets[ijet];

	bool passID = passPFJetID(jetID, jet);

	if (passID)
		{
		IDjets.push_back(jet);
		if (record)
			{
			fill_2d(string("control_jet_jetsIDed_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_jetsIDed_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		}
	}

// The jet corrections

// v6, adding jet corrections and b-tagging
//LorentzVector full_jet_corr(0., 0., 0., 0.);
for(size_t ijet=0; ijet<IDjets.size(); ijet++)
	{
	// TODO: so does this mean "in place"?
	pat::Jet& jet = IDjets[ijet];

	// for MC smearing
	//TRandom *r3 = new TRandom3();

	if (record)
		{
		//fill_2d(string("slimmedjet_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_2d(string("control_jet_slimmedjet_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(string("control_jet_slimmedjet_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}

	LorentzVector jet_corr(0., 0., 0., 0.);

	LorentzVector jet_initial_momentum = jet.p4();

	if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
	//correct JES
	LorentzVector rawJet = jet.correctedP4("Uncorrected");

	if (record)
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

	if (record)
		{
		fill_2d(string("control_jet_slimmedjet_jescor_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(string("control_jet_slimmedjet_jescorrection"), 400, 0., 2., jes_correction, weight);
		}

	if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

	//smear JER
	//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
	// if a jet matches well to a genJet -- it can be just scaled and that's fine
	// if id doesn't -- one needs to smear it with a randomizer according to the resolution and other stuff

	// here is the matching of the jet:
	if(isMC)
		{
		if (record)
		       {
		       fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_before"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
		       }
		// the JER SF and resolution for the jet:
		//std::vector<double> jer_sf_pair = JER_SF(jet.eta());
		//double jer_sf = jer_sf_pair[0];
		//double jer_resolution = jer_sf_pair[1]; // TODO: not sure about this -- is the table the same as what their tool returns?
		// getting it with the tool from files:
		double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho}});
		double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, m_systematic_variation);

		// matching to generation jet:
		//const reco::GenJet* genJet=jet.genJet();
		// the PAT doc describes it as "return the matched generated jet"
		// what's the matching procedure?
		// for now I'll do it manually in dR, as shown in Jet POG example
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
		const reco::GenJet* matched_genJet = nullptr;
		//double dR_max = 0.4/2; // 0.4 is the jet cone parameter of AK4 jets, which I use
		// moved it to parameters of the procedure
		for (int i=0; i<genJets.size(); i++)
			{
			reco::GenJet& genJet = genJets[i];
			double dR = reco::deltaR(jet, genJet);

			if (dR > dR_max) continue;

			double dPt = std::abs(genJet.pt() - jet.pt());
			double dPt_max_factor = 3*jet.pt(); // from twiki
			if (dPt > dPt_max_factor * jet_resolution) continue;

			matched_genJet = &genJet;
			}

		if (matched_genJet)
			{ // the scaling is ok
			double dPt = jet.pt() - matched_genJet->pt();
			//double genjetpt( genJet ? genJet->pt(): 0.);                    
			//std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
			// using the local smear:
			//std::vector<double> JER_smearing_factor = smearJER(jet.pt(),jet.eta(),genjetpt);
			//double jer_smearing = JER_smearing_factor[0];
			double jer_smearing = TMath::Max(0., 1. + (jer_sf - 1) * dPt / jet.pt());
			jet.setP4(jet.p4()*jer_smearing); // same as scaleEnergy in the Jet POG example
			// but they also do MIN_ENERGY thing
			// which is static constexpr const double MIN_JET_ENERGY = 1e-2;
			fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_scaling"), 400, 0., 2., jer_smearing, weight);

			if (record)
				{
				fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_scaling_done"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
				}
			}
		else
			{ // the smearing with gaussian randomizer
			// this is the example:
			//double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
			//double smearFactor = 1 + r3->Gaus(0, sigma);
			// this is the twiki:
			double smearFactor = 1 + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf*jer_sf - 1.));
			jet.setP4(jet.p4()*TMath::Max(0., smearFactor));
			fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_stochastic_smearing"), 400, 0., 2., smearFactor, weight);

			if (record)
				{
				fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_smearing_done"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
				}
			}
		}

	/*
	if(isMC)
		{
		const reco::GenJet* genJet=jet.genJet();
		// the PAT doc describes it as "return the matched generated jet"
		// what's the matching procedure?
		// for now I'll do it manually in dR, as shown in Jet POG example

		if(genJet)
			{
			double genjetpt( genJet ? genJet->pt(): 0.);                    
			//std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
			// using the local smear:
			std::vector<double> JER_smearing_factor = smearJER(jet.pt(),jet.eta(),genjetpt);
			double jer_smearing = JER_smearing_factor[0];
			jet.setP4(jet.p4()*jer_smearing);
			fill_1d(string("control_jet_slimmedjet_mc_jersmearing"), 400, 0., 2., jer_smearing, weight);

			//printf("jet pt=%f gen pt = %f smearing %f %f %f\n", jet.pt(), genjetpt, JER_smearing_factor[0], JER_smearing_factor[1], JER_smearing_factor[2]);
			// //set the JER up/down alternatives
			jet.addUserFloat("jerup", JER_smearing_factor[1]);  //kept for backward compatibility
			jet.addUserFloat("jerdown", JER_smearing_factor[2] ); //kept for backward compatibility
			jet.addUserFloat("_res_jup", JER_smearing_factor[1]);
			jet.addUserFloat("_res_jdown", JER_smearing_factor[2] );
			}
		else{
			jet.addUserFloat("jerup", 1.0); //kept for backward compatibility
			jet.addUserFloat("jerdown", 1.0);  //kept for backward compatibility
			jet.addUserFloat("_res_jup", 1.0);
			jet.addUserFloat("_res_jdown", 1.0 );
			}
		if (record)
			fill_2d(string("control_jet_slimmedjet_mc_jercor_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
		}
	*/


	// here is the correct3 jet correction point

	/*
	if(isMC)
		{
		////set the JES up/down pT alternatives
		std::vector<float> ptUnc = utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
		jet.addUserFloat("jesup",    ptUnc[0] );  //kept for backward compatibility
		jet.addUserFloat("jesdown",  ptUnc[1] );  //kept for backward compatibility
		jet.addUserFloat("_scale_jup",    ptUnc[0] );
		jet.addUserFloat("_scale_jdown",  ptUnc[1] );
		}
	*/

	// FIXME: this is not to be re-set. Check that this is a desired non-feature.
	// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
	//to get the raw jet again
	//IDjets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));

	// Add the jet momentum correction:
	// jet_cor propagation is on in 13.4
	jet_corr = jet.p4() - jet_initial_momentum;
	full_jet_corr += jet_corr;

	if (record)
		{
		fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_eta"), 400, 0., 400., 200, -4., 4., full_jet_corr.pt(), full_jet_corr.eta(), weight);
		// so, jes & jer shouldn't change the eta of the jet -- only the energy (pt)
		// record the correction per jet's eta and per original pt
		fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_eta"), 400, 0., 400., 200, -4., 4., full_jet_corr.pt(), jet.eta(), weight);
		fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_pt"), 400, 0., 400., 200, -4., 4.,  full_jet_corr.pt(), jet_initial_momentum.eta(), weight);
		}

	if(debug)
		{
		cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
		cout << "-" << jet_initial_momentum.eta() << " " << jet_initial_momentum.pt() << " " << jet_initial_momentum.energy() << endl;
		}
	}


std::sort (IDjets.begin(),  IDjets.end(),  utils::sort_CandidatesByPt);

// ----------------------------------- here is the correctF jet correction point
// Propagate full_jet_corr to MET:
//met.setP4(met.p4() - full_jet_corr); // just return the full correction and propagate in place
// TODO: uncertainties?
// for leptons they are done as:
//met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
//met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction



// FIXME: So are these MET corrections?
//if(debug) cout << "Update also MET" << endl;
// LorentzVector n_met = met.p4();
// std::vector<LorentzVector> newMet = utils::cmssw::getMETvariations(n_met/*recoMet*/,IDjets,selLeptons,isMC);
// FIXME: Must choose a lepton collection. Perhaps loose leptons?
// n_met = newMet[utils::cmssw::METvariations::NOMINAL];

//fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), n_met.pt(), weight);
//if (isSingleMu) fill_pt_e( string("met0_all_leptoncorr_jetcorr_singlemu_pt"), n_met.pt(), weight);
//if (isSingleE)  fill_pt_e( string("met0_all_leptoncorr_jetcorr_singleel_pt"), n_met.pt(), weight);
//fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), met.pt(), weight);
//fill_1d(string("control_met_allcorr_pt"), 200, 0., 200., met.pt(), weight);

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
// the kinematics only, now (at Moriond17) the IDs are recommended to be applied before corrections

//pat::JetCollection selJets;
//pat::JetCollection selJets20GeV, selJets30GeV;
// TODO: do all jet selection right here
// now selBJets are not used anywhere
// selJets pass cross-cleaning with taus later
// and b-tagging again
double mindphijmet (9999.);
for (unsigned int count_ided_jets = 0, ijet = 0; ijet < IDjets.size(); ++ijet)
	{
	pat::Jet& jet = IDjets[ijet];

	if (record)
		{
		fill_2d(string("control_jet_jetscorrected_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(string("control_jet_jetscorrected_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}

	// TODO: what do we do here exactly?
	// a loose selection on jets, and then tighten it later?
	// if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue;
	// Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)

	//mc truth for this jet
	//const reco::GenJet * genJet = jet.genJet();
	//TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");
	// TODO: this mctruth for jets it is never used in the code

	//jet id (done before corrections)
	//bool passPFloose = passPFJetID(jetID, jet); 
	// bool passPFloose = passPFJetID("Tight", jet); 
	//if (label=="Tight")
	// FIXME: check when pileup ID will come out

	// Jet Kinematics
	double eta = jet.eta();
	double pt  = jet.pt();
	//bool passKino = pt > 30. && fabs(eta) < 2.4;

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


	//if (passPFloose && passKino)
	//if (passPFloose && (fabs(eta) < eta_cut) && pt > pt_cut)
	// now the ID is applied before even corrections
	if ((fabs(eta) < eta_cut) && pt > pt_cut)
		{
		//if (pt > 30.)
			{
			// trying new fake rate application:
			// if there are only 2 30GeV jets -- they are not considered for faking taus
			// TODO: and what to do with 20GeV jets? in selection and in fake rates?
			//  --- for now I will just add 20-30 GeV jets for fake rates..
			//      one should actually try without them
			//selJets30GeV.push_back(jet);
			// drop the 20-30GeV jets for now
			selJets.push_back(jet);

			if (record)
				{
				fill_2d(string("control_jet_selJets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
				fill_1d(string("control_jet_selJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
				}
			// double dphijmet = fabs (deltaPhi (n_met.phi(), jet.phi()));
			//double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
			//if (dphijmet < mindphijmet) mindphijmet = dphijmet;
			// FIXME: mindphijmet is not used anywhere now
			}
		//else if (pt > 20.)
			//{
			//selJets20GeV.push_back(jet);
			//}
		}
	}

std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);

return 0;
}



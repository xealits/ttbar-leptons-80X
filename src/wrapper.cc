#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include <string.h>

int foo()
{
return 5;
}

/* Recoil corrections
 */

static TString recoil_corrections_data_file("/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root");
//TString recoil_corrections_data_file("/HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root");
static RecoilCorrector* recoilPFMetCorrector = new RecoilCorrector(recoil_corrections_data_file);

float met_pt_recoilcor(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return sqrt(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex);
	}

float met_pt_recoilcor_x(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return NT_pfmetcorr_ex;
	}

float met_pt_recoilcor_y(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return NT_pfmetcorr_ey;
	}


// Z pt-mass weight
// zpt_weights_2016.root
TFile* zPtMassWeights_file  = TFile::Open("/UserCode/zpt_weights_2016.root");
TH2D*  zPtMassWeights_histo = (TH2D*) zPtMassWeights_file->Get("zptmass_histo");
TH2D*  zPtMassWeights_histo_err = (TH2D*) zPtMassWeights_file->Get("zptmass_histo_err");

float zPtMass_weight(float genMass, float genPt)
	{
	return zPtMassWeights_histo->GetBinContent(zPtMassWeights_histo->GetXaxis()->FindBin(genMass), zPtMassWeights_histo->GetYaxis()->FindBin(genPt));
	}


/*
double transverse_mass_pts(double v1_x, double v1_y, double v2_x, double v2_y)
	{
	double v1v2 = sqrt((v1_x*v1_x + v1_y*v1_y)*(v2_x*v2_x + v2_y*v2_y));
	return sqrt(2*(v1v2 - (v1_x*v2_x + v1_y*v2_y)));
	}
*/

float MTlep_met_pt_recoilcor(float lep_px, float lep_py,
	float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	//return sqrt(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex);
	double v1v2 = sqrt((lep_px*lep_px + lep_py*lep_py)*(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex));
	return sqrt(2*(v1v2 - (lep_px*NT_pfmetcorr_ex + lep_py*NT_pfmetcorr_ey)));
	}


/*
// b-tagging SF
// for each jet check if it is tagged
// apply SF accordingly
// 
*/

//Float_t NT_jet_id_0;
//Float_t NT_jet_eta_0;
//Float_t NT_jet_phi_0;
//Float_t NT_jet_pt_0;
//Float_t NT_jet_p_0;
//Float_t NT_jet_rad_0;
//Float_t NT_jet_b_discr_0;
//int flavId = (int) NT_jet_hadronFlavour_0;
//Float_t NT_jet_partonFlavour_0;
//double sf = 1.0, eff = 1.0;

/*
struct bTaggingEfficiencyHistograms bEffs = {
	.b_alljet    = bTaggingEfficiencies_b_alljet   ,
	.b_tagged    = bTaggingEfficiencies_b_tagged   ,
	.c_alljet    = bTaggingEfficiencies_c_alljet   ,
	.c_tagged    = bTaggingEfficiencies_c_tagged   ,
	.udsg_alljet = bTaggingEfficiencies_udsg_alljet,
	.udsg_tagged = bTaggingEfficiencies_udsg_tagged
	};
*/

struct bTaggingEfficiencyHistograms {
	TH2F* b_alljet   ;
	TH2F* b_tagged   ;
	TH2F* c_alljet   ;
	TH2F* c_tagged   ;
	TH2F* udsg_alljet;
	TH2F* udsg_tagged;
	};


static double bTagging_b_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.b_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.b_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.b_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_2d(string("btag_eff_retrieved0_flavour_b"), 250, 0., 500., 200, -4., 4., pt, eta, 1);
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_b"), pt, eta, 1);
		return 0; // whatch out -- equality of floats
		}
	return N_tagged/N_alljets;
	}

static double bTagging_c_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.c_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.c_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.c_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_c"), pt, eta, 1);
		return 0;
		}
	return N_tagged/N_alljets;
	}

static double bTagging_udsg_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.udsg_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.udsg_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.udsg_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_udsg"), pt, eta, 1);
		return 0;
		}
	return N_tagged/N_alljets;
	}


/* ------------------------------------------------------------------
 * B-TAGGING SFs
 */
static std::string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
//float b_tag_WP = 0.8484; // medium

/*
void load_btagCal()
	{
	btagCalib = new BTagCalibration("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");
	btagCal.load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_B,      // btag flavour
	//            "comb");              // they say "comb" is better precision, but mujets are independent from ttbar dilepton channels
		"mujets");                //
	btagCal.load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_C,      // btag flavour
		"mujets");               // measurement type
	btagCal.load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_UDSG,   // btag flavour
		"incl");                // measurement type
	}
*/

//gROOT->ProcessLine("#include <vector>");

static double calc_btag_sf_weight(BTagCalibrationReader& btagCal, struct bTaggingEfficiencyHistograms& bEffs, bool hasCSVtag, int flavId, double pt, double eta)
	{
	// int flavId=jet.partonFlavour();
	//int flavId=jet.hadronFlavour();
	// also: patJet->genParton().pdgId()
	// fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

	double sf = 1.0, eff = 1.0;
	/* If the jet is tagged -- weight *= SF of the jet
	 * if not weight *= (1 - eff*SF)/(1 - eff)
	 */

	if (abs(flavId)==5) {
		// get SF for the jet
		sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt, 0.);
		// get eff for the jet
		eff = bTagging_b_jet_efficiency(bEffs, pt, eta);
		}
	else if(abs(flavId)==4) {
		sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt, 0.);
		eff = bTagging_c_jet_efficiency(bEffs, pt, eta);
		}
	else {
		sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt, 0.);
		eff = bTagging_udsg_jet_efficiency(bEffs, pt, eta);
		}

	double jet_weight_factor = 1;
	if (hasCSVtag) // a tagged jet
		{
		jet_weight_factor = sf;
		}
	else // not tagged
		{
		// truncate efficiency to 0 and 0.99
		eff = (eff < 0. ? 0. : (eff > 0.99 ? 0.99 : eff));
		jet_weight_factor = (1 - sf*eff) / (1 - eff);
		}

	return jet_weight_factor;
	}

// these are members of our object file
static std::string bTaggingEfficiencies_filename;
static TFile* bTaggingEfficiencies_file;
static TString backup_b_eff_distr;
static TH2F* bTaggingEfficiencies_b_alljet   ;
static TH2F* bTaggingEfficiencies_b_tagged   ;
static TH2F* bTaggingEfficiencies_c_alljet   ;
static TH2F* bTaggingEfficiencies_c_tagged   ;
static TH2F* bTaggingEfficiencies_udsg_alljet;
static TH2F* bTaggingEfficiencies_udsg_tagged;
static struct bTaggingEfficiencyHistograms bEffs;
static BTagCalibration* btagCalib;
static BTagCalibrationReader* btagCal;

void set_bSF_calibrators()
	{
	bTaggingEfficiencies_filename = std::string(std::getenv("CMSSW_BASE")) + "/src/UserCode/ttbar-leptons-80X/jobing/configs/b-tagging-efficiencies.root";
	bTaggingEfficiencies_file = TFile::Open(bTaggingEfficiencies_filename.c_str());

	btagCalib = new BTagCalibration("CSVv2", std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");
	//BTagCalibration* btagCalib; // ("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");
	btagCal = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
	// BTagCalibrationReader btagCal(BTagEntry::OP_TIGHT,  // operating point
				     "central",             // central sys type
				     {"up", "down"});      // other sys types

	btagCal->load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_B,      // btag flavour
	//            "comb");              // they say "comb" is better precision, but mujets are independent from ttbar dilepton channels
		"mujets");                //
	btagCal->load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_C,      // btag flavour
		"mujets");               // measurement type
	btagCal->load(*btagCalib,              // calibration instance
		BTagEntry::FLAV_UDSG,   // btag flavour
		"incl");                // measurement type
	}

void set_bSF_effs_for_dtag(TString dtag)
	{
	// read-in and setup all the b-tagging crap

	//TString bTaggingEfficiencies_filename   = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/jobing/configs/b-tagging-efficiencies.root";
	//gSystem->ExpandPathName(bTaggingEfficiencies_filename);

	//cout << "b-tagging eff-s, filename: " << bTaggingEfficiencies_filename << endl;

	backup_b_eff_distr = "MC2016_Summer16_DYJetsToLL_50toInf_madgraph";

	//TString dtag("foo");
	// ( ? : ) would look too much here
	if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_b_hadronFlavour_candidates_tagged"))
		{
		bTaggingEfficiencies_b_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_b_hadronFlavour_candidates");
		bTaggingEfficiencies_b_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_b_hadronFlavour_candidates_tagged");
		}
	else
		{
		bTaggingEfficiencies_b_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates");
		bTaggingEfficiencies_b_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates_tagged");
		}
	if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_c_hadronFlavour_candidates_tagged"))
		{
		bTaggingEfficiencies_c_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_c_hadronFlavour_candidates");
		bTaggingEfficiencies_c_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_c_hadronFlavour_candidates_tagged");
		}
	else
		{
		bTaggingEfficiencies_c_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates");
		bTaggingEfficiencies_c_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates_tagged");
		}
	if (bTaggingEfficiencies_file->GetListOfKeys()->Contains(dtag + "_btag_udsg_hadronFlavour_candidates_tagged"))
		{
		bTaggingEfficiencies_udsg_alljet = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_udsg_hadronFlavour_candidates");
		bTaggingEfficiencies_udsg_tagged = (TH2F*) bTaggingEfficiencies_file->Get(dtag + "_btag_udsg_hadronFlavour_candidates_tagged");
		}
	else
		{
		bTaggingEfficiencies_udsg_alljet = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates");
		bTaggingEfficiencies_udsg_tagged = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates_tagged");
		}

	/*
	bTaggingEfficiencies_b_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates");
	bTaggingEfficiencies_b_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_b_hadronFlavour_candidates_tagged");
	bTaggingEfficiencies_c_alljet    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates");
	bTaggingEfficiencies_c_tagged    = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_c_hadronFlavour_candidates_tagged");
	bTaggingEfficiencies_udsg_alljet = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates");
	bTaggingEfficiencies_udsg_tagged = (TH2F*) bTaggingEfficiencies_file->Get(backup_b_eff_distr + "_btag_udsg_hadronFlavour_candidates_tagged");
	*/

	bEffs.b_alljet    = bTaggingEfficiencies_b_alljet   ;
	bEffs.b_tagged    = bTaggingEfficiencies_b_tagged   ;
	bEffs.c_alljet    = bTaggingEfficiencies_c_alljet   ;
	bEffs.c_tagged    = bTaggingEfficiencies_c_tagged   ;
	bEffs.udsg_alljet = bTaggingEfficiencies_udsg_alljet;
	bEffs.udsg_tagged = bTaggingEfficiencies_udsg_tagged;
	}

double b_taggin_SF (double jet_pt, double jet_eta, double jet_b_discr, int jet_hadronFlavour, double b_tag_WP)
	{
	double weight_btag_sf = 1;

	//FLoat_t b_discriminator = NT_jet_b_discr_0;
	//bool hasCSVtag(b_discriminator > b_tag_WP);
	if (jet_pt > 0) // jet exists
		weight_btag_sf *= calc_btag_sf_weight(*btagCal, bEffs, jet_b_discr > b_tag_WP, jet_hadronFlavour, jet_pt, jet_eta);
	/*
	if (jet1_p4.pt() > 0) // jet exists                                                                                    
		weight_btag_sf *= calc_btag_sf_weight(btagCal, jet1_b_discr > b_tag_WP, (int) jet1_hadronFlavour, jet1_pt, jet1_eta);
	if (jet2_p4.pt() > 0) // jet exists                                                                                    
		weight_btag_sf *= calc_btag_sf_weight(btagCal, jet2_b_discr > b_tag_WP, (int) jet2_hadronFlavour, jet2_pt, jet2_eta);
	if (jet3_p4.pt() > 0) // jet exists                                                                                    
		weight_btag_sf *= calc_btag_sf_weight(btagCal, jet3_b_discr > b_tag_WP, (int) jet3_hadronFlavour, jet3_pt, jet3_eta);
	if (jet4_p4.pt() > 0) // jet exists                                                                                    
		weight_btag_sf *= calc_btag_sf_weight(btagCal, jet4_b_discr > b_tag_WP, (int) jet4_hadronFlavour, jet4_pt, jet4_eta);
	*/

	//weight *= weight_btag_sf;
	return weight_btag_sf;
	}


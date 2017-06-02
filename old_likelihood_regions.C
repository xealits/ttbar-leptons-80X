//#include <iostream>
//
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
////#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit
//
//#include "TSystem.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TCanvas.h"
//#include "TLegend.h"
//#include "THStack.h"
//#include "TH1F.h"                                                                                                                                       
//#include "TH2F.h"
//#include "TH3F.h"
//#include "TProfile.h"
//#include "TProfile2D.h"
//#include "TEventList.h"
//#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>
#include "TMath.h"
#include "Math/GenVector/LorentzVector.h"
//
//#include "TStyle.h"
//
//#include <map>
//#include <string>
//#include <vector>
//
//#include "dtag_xsecs.h"
//
//#include "TNtuple.h"
//
//
//using namespace std;

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

double pileup_ratio[] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0};


struct bTaggingEfficiencyHistograms {
	TH2F* b_alljet   ;
	TH2F* b_tagged   ;
	TH2F* c_alljet   ;
	TH2F* c_tagged   ;
	TH2F* udsg_alljet;
	TH2F* udsg_tagged;
	};

double bTagging_b_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
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

double bTagging_c_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
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

double bTagging_udsg_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
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

struct bTaggingEfficiencyHistograms bEffs;

double calc_btag_sf_weight(BTagCalibrationReader& btagCal, bool hasCSVtag, int flavId, double pt, double eta)
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






int setup_distrs_for_mc_channel(map<TString, TH1D*>& event_yields_in_cathegories,
	map<TString, map<TString, map<TString, TH1D*>>>& distrs_in_cathegories,
	vector<TString>& cathegories,
	TString& mc_channel)
	{
	cout << "for channel " << mc_channel << " setting distrs in cathegories";
	for (int i=0; i<cathegories.size(); i++) cout << " " << cathegories[i];
	cout << endl;

	Int_t nx = cathegories.size();
	event_yields_in_cathegories[mc_channel] = new TH1D(mc_channel + "event_yields", ";;", nx, 0, nx);
	/* somehow the labels don't stack or whatever
	for (int i=1; i<=nx; i++)
		event_yields_in_cathegories[mc_channel]->GetXaxis()->SetBinLabel(i, cathegories[i-1]);
	*/
	// set up other channels dynamically in the loop in the same fashion

	// same MC channels here
	for (int i=0; i<cathegories.size(); i++)
		{
		distrs_in_cathegories[mc_channel][cathegories[i]]["lj_peak_distance_inclusive"] = new TH1D(mc_channel + cathegories[i] + "lj_peak_distance_inclusive", ";;", 100, 0, 10000);
		distrs_in_cathegories[mc_channel][cathegories[i]]["lj_peak_distance"]    = new TH1D(mc_channel + cathegories[i] + "lj_peak_distance", ";;", 100, 0, 10000);

		distrs_in_cathegories[mc_channel][cathegories[i]]["tau_energy_ratio"]    = new TH1D(mc_channel + cathegories[i] + "tau_energy_ratio", ";;", 40, 0, 1.2);
		distrs_in_cathegories[mc_channel][cathegories[i]]["met_pt"] = new TH1D(mc_channel + cathegories[i] + "met_pt", ";;", 40, 40, 200);
		distrs_in_cathegories[mc_channel][cathegories[i]]["tau_pt"]              = new TH1D(mc_channel + cathegories[i] + "tau_pt", ";;", 100, 0, 200);
		distrs_in_cathegories[mc_channel][cathegories[i]]["lep_pt"]              = new TH1D(mc_channel + cathegories[i] + "lep_pt", ";;", 100, 0, 200);
		distrs_in_cathegories[mc_channel][cathegories[i]]["jet_pt"]              = new TH1D(mc_channel + cathegories[i] + "jet_pt", ";;", 100, 0, 200);
		distrs_in_cathegories[mc_channel][cathegories[i]]["tau_eta"]             = new TH1D(mc_channel + cathegories[i] + "tau_eta", ";;", 100, -2.5, 2.5);
		distrs_in_cathegories[mc_channel][cathegories[i]]["lep_eta"]             = new TH1D(mc_channel + cathegories[i] + "lep_eta", ";;", 100, -2.5, 2.5);
		distrs_in_cathegories[mc_channel][cathegories[i]]["jet_eta"]             = new TH1D(mc_channel + cathegories[i] + "jet_eta", ";;", 100, -2.5, 2.5);

		distrs_in_cathegories[mc_channel][cathegories[i]]["nvtx"]     = new TH1D(mc_channel + cathegories[i] + "nvtx", ";;", 50, 0, 50);
		distrs_in_cathegories[mc_channel][cathegories[i]]["rho"]      = new TH1D(mc_channel + cathegories[i] + "rho", ";;", 50, 0, 50);
		distrs_in_cathegories[mc_channel][cathegories[i]]["nvtx_raw"] = new TH1D(mc_channel + cathegories[i] + "nvtx_raw", ";;", 50, 0, 50);
		distrs_in_cathegories[mc_channel][cathegories[i]]["rho_raw"]  = new TH1D(mc_channel + cathegories[i] + "rho_raw", ";;", 50, 0, 50);

		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_pu"]    = new TH1D(mc_channel + cathegories[i] + "weight_pu", ";;", 100, -1.2, 1.2);

		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_muon_sfs"]   = new TH1D(mc_channel + cathegories[i] + "weight_muon_sfs", ";;", 100, -1.2, 1.2);
		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_muon_trig"]  = new TH1D(mc_channel + cathegories[i] + "weight_muon_trig", ";;", 100, -1.2, 1.2);
		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_btag_sf"]    = new TH1D(mc_channel + cathegories[i] + "weight_btag_sf", ";;", 100, -1.2, 1.2);
		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_without_pu"] = new TH1D(mc_channel + cathegories[i] + "weight_without_pu", ";;", 100, -1.2, 1.2);
		distrs_in_cathegories[mc_channel][cathegories[i]]["weight_final"]      = new TH1D(mc_channel + cathegories[i] + "weight_final", ";;", 100, -1.2, 1.2);
		}

	return 0;
	}

#define INPUT_DTAGS_START 5
#define NTUPLE_NAME "ntuple"

int old_likelihood_regions(bool be_verbose, Long64_t event_prescale, TString dir, TString dtag)
{
//bool be_verbose = true;
//double lumi = 4631.304;
//Long64_t event_prescale = 50;
//TString dir("outdir/v11.6/merged-sets/");
//TString dtag("Data13TeV_SingleMuon2016B_03Feb2017_ver2");

char usage_string[128] = "be_verbose event_prescale dir dtags";
//if (argc < INPUT_DTAGS_START)
//	{
//	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
//	return 1;
//	}

gROOT->Reset();

/*
TString verb(argv[1]);
bool be_verbose = false;
if (verb==TString("T") || verb==TString("Y"))
	be_verbose = true;
double lumi = atof(argv[2]);
Long64_t event_prescale = atoi(argv[3]);
TString dir(argv[4]);
TString dtag1(argv[5]);
*/

TString dtag_file = dir + "/" + dtag + ".root";
cout << be_verbose  << endl;
cout << event_prescale   << endl;
cout << dtag_file << endl;


/* ------------------------------------------------------------------
 * Muon ID SFs
 */

TString muon_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/muon-effs/";
gSystem->ExpandPathName(muon_effs_dirname    );

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

/* ------------------------------------------------------------------
 * Muon trigger SF
 */

TFile* muon_effs_trg_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_RunBtoF.root").Data() );
TFile* muon_effs_trg_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_Period4.root").Data() );
TH2D* muon_effs_trg_BCDEF_histo = (TH2D*) ((TDirectoryFile*) muon_effs_trg_BCDEF_file->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins"))->Get("pt_abseta_ratio");
TH2D* muon_effs_trg_GH_histo    = (TH2D*) ((TDirectoryFile*) muon_effs_trg_GH_file->   Get("IsoMu24_OR_IsoTkMu24_PtEtaBins"))->Get("pt_abseta_ratio");
cout << "Y trigger" << endl;


/* ------------------------------------------------------------------
 * B-TAGGING SFs
 */
string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
float b_tag_WP = 0.8484; // medium

TString bTaggingEfficiencies_filename   = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/jobing/configs/b-tagging-efficiencies.root";
gSystem->ExpandPathName(bTaggingEfficiencies_filename);
TFile* bTaggingEfficiencies_file = TFile::Open(bTaggingEfficiencies_filename.Data());

cout << "b-tagging eff-s, filename: " << bTaggingEfficiencies_filename << endl;

TH2F* bTaggingEfficiencies_b_alljet   ;
TH2F* bTaggingEfficiencies_b_tagged   ;
TH2F* bTaggingEfficiencies_c_alljet   ;
TH2F* bTaggingEfficiencies_c_tagged   ;
TH2F* bTaggingEfficiencies_udsg_alljet;
TH2F* bTaggingEfficiencies_udsg_tagged;

TString backup_b_eff_distr("MC2016_Summer16_DYJetsToLL_50toInf_madgraph");

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


//struct bTaggingEfficiencyHistograms bEffs;

bEffs.b_alljet    = bTaggingEfficiencies_b_alljet   ;
bEffs.b_tagged    = bTaggingEfficiencies_b_tagged   ;
bEffs.c_alljet    = bTaggingEfficiencies_c_alljet   ;
bEffs.c_tagged    = bTaggingEfficiencies_c_tagged   ;
bEffs.udsg_alljet = bTaggingEfficiencies_udsg_alljet;
bEffs.udsg_tagged = bTaggingEfficiencies_udsg_tagged;

BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");

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

// ------------------------------------------------------------------




TFile* file = TFile::Open(dtag_file);

if (!file->GetListOfKeys()->Contains(NTUPLE_NAME))
	{
	if (be_verbose) cout << "no " << NTUPLE_NAME << endl;
	return 1;
	}

//TH1D * eventflow = (TH1D*) file->Get("eventflow");
TH1D * eventflow = (TH1D*) file->Get("weightflow_elel_NOMINAL");

// it's actually TNtuple, which descends from TTree thus should open ok
TTree* NT_output_ttree = (TTree*) file->Get(NTUPLE_NAME);
cout << "ttree pointer " << NT_output_ttree << endl;
//NT_output_ttree->Print();

if (be_verbose)
	{
	cout << "getting the cathegories " << NT_output_ttree->GetEntries() << '\t' << event_prescale << '\t' << NT_output_ttree->GetEntries() / event_prescale;
	cout << endl;
	}

//TString histo_name = dtag + TString("_cathegories");
// the histogram with yields in event categories
/*
const Int_t nx = 15;
const char *cathegories[nx] = {"1b3j", "1b4j", "1b5j", "1b6j",
	"2b3j", "2b4j", "2b5j", "2b6j",
	"lj_1b3j", "lj_1b4j", "lj_1b5j", "lj_1b6j",
	 "lj_2b4j", "lj_2b5j", "lj_2b6j"}
*/

vector<TString> cathegories;
cathegories.push_back("1b3j");
cathegories.push_back("1b4j");
cathegories.push_back("1b5j");
cathegories.push_back("2b3j");
cathegories.push_back("2b4j");
cathegories.push_back("2b5j");
cathegories.push_back("lj_peak");

	/*
	"lj_1b3j", "lj_1b4j", "lj_1b5j",
	//lj_2b3j",
	"lj_2b4j", "lj_2b5j"};
	*/

Int_t nx = cathegories.size();

// everything is per MC decay channel
map<TString, TH1D*> event_yields_in_cathegories;
map<TString, map<TString, map<TString, TH1D*>>> distrs_in_cathegories;

//cout << "frame pull cathegories" << endl;
//cout << "from ttree" << endl;
//NT_output_ttree->Print();
// the interface to NT_output_ttree
#define NTUPLE_INTERFACE_OPEN
#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_v11_3.h"

bool is_MC = dtag.Contains("MC");
bool is_W0Jets = dtag.Contains("W0Jet"); // needed for handling WNJets
bool is_aMCatNLO = dtag.Contains("amcatnlo");
bool is_TTbar = dtag.Contains("TT");

NT_output_ttree->GetEntry(0);

Long64_t n_events = (Long64_t) NT_output_ttree->GetEntries() / event_prescale;
Long64_t one100th_events = (Long64_t) NT_output_ttree->GetEntries() / 100;
Long64_t j = 0;
cout << "to the loop" << endl;
/* Loop through events, find weights, skip stuff, record if they pass to any category
 */
for (Long64_t i=0; j<n_events && i<NT_output_ttree->GetEntries(); i++)
//for (Long64_t i=0; i<10; i++)
	{
	//cout << i;
	NT_output_ttree->GetEntry(i);
	// test:
	//cout << NT_NUP_gen << '\t' << NT_aMCatNLO_weight << '\t' << NT_nvtx_gen << '\t' << (int) NT_nvtx_gen << endl;
	//cout << "\tgot entry" << endl;

	if (NT_njets==1) continue; // the prescale of large jobs run
	j++;

	TString mc_decay("");
	int event_category = -1;
	double weight_without_pu =1, weight = 1; // for MC
	double weight_muon_sfs = 1, weight_muon_trig = 1;
	double weight_btag_sf = 1;

	//cout << NT_met_corrected << endl;
	//cout << *NT_met_corrected << endl;

	// general requirement for all events:
	//if (NT_met_corrected->pt() < 40 || fabs(NT_leps_ID) != 13 || NT_lep0_id * NT_tau0_id > 0 || NT_njets < 3 || NT_nbjets < 1) continue;
	if (NT_met_corrected < 40 || NT_HLT_mu != 1 || fabs(NT_leps_ID) != 13 || NT_lep_id_0 * NT_tau_id_0 > 0 || NT_njets < 3 || NT_nbjets < 1) continue;
	// the event with taus (2 = medium, 1 = loose):
	if (NT_tau_IDlev_0 < 1. && NT_tau_IDlev_1 < 1.) continue;

	// find basic weight:
	//  * -1 aMCatNLO & NUP for WJets
	//  * pileup

	// NUP == 5 for W0Jets
	if (is_W0Jets && NT_NUP_gen != 5) continue;

	if (is_TTbar)
		{
		// set mc_decay
		int decay_id = abs(NT_gen_t_w_decay_id * NT_gen_tb_w_decay_id);
		// = id of lepton (11/13/15, but the sign means which product is lepton: minus=1, plus=2) or 1 for quarks

		if (decay_id == 15*15) mc_decay = TString("_aattuu");
		else if (decay_id == 15*11) mc_decay = TString("_aeltu");
		else if (decay_id == 15*13) mc_decay = TString("_amtuu");
		else if (decay_id == 15*1)  mc_decay = TString("_aqtu");
		else if (decay_id == 11*11) mc_decay = TString("_eell");
		else if (decay_id == 11*13) mc_decay = TString("_elmu");
		else if (decay_id == 11*1)  mc_decay = TString("_elq");
		else if (decay_id == 13*13) mc_decay = TString("_mmuu");
		else if (decay_id == 13)    mc_decay = TString("_mqu");
		else if (decay_id == 1)     mc_decay = TString("_qq");
		}

	// amcatnlo MC has these weights
	if (is_aMCatNLO)
		weight *= (NT_aMCatNLO_weight > 0? 1 : -1);

	if (is_MC)
		{
		weight *= 0.97; // tau ID SF, Medium = 0.97
		//weight *= 0.98; // average mu trig SF
		//weight *= 0.97; // average mu ID

		// muon ID SFs:
		//double weight_muon_sfs = 1;
		double bcdef_weight = 1;
		double gh_weight = 1;
		double SingleMuon_data_bcdef_fraction = 19716.274 / (19716.274 + 15931.028);
		double SingleMuon_data_gh_fraction    = 15931.028 / (19716.274 + 15931.028);

		double abs_eta = abs(NT_lep_eta_0);
		double pt = NT_lep_pt_0;

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

		//fill_1d(string("weight_muon_effs_BCDEF_trk"),  200, 0., 1.1,   bcdef_weight_trk, 1);
		//fill_1d(string("weight_muon_effs_BCDEF_id"),   200, 0., 1.1,   bcdef_weight_id,  1);
		//fill_1d(string("weight_muon_effs_BCDEF_iso"),  200, 0., 1.1,   bcdef_weight_iso, 1);
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

		//fill_1d(string("weight_muon_effs_GH_trk"),  200, 0., 1.1,   gh_weight_trk, 1);
		//fill_1d(string("weight_muon_effs_GH_id"),   200, 0., 1.1,   gh_weight_id,  1);
		//fill_1d(string("weight_muon_effs_GH_iso"),  200, 0., 1.1,   gh_weight_iso, 1);
		gh_weight *= gh_weight_trk * gh_weight_id * gh_weight_iso;
		// and merge them:
		weight_muon_sfs = SingleMuon_data_bcdef_fraction * bcdef_weight + SingleMuon_data_gh_fraction * gh_weight;

		weight *= weight_muon_sfs;


		// MUON Trigger
		{
		double no_mu_trig = 1;
		//double mu_trig_weight = 1;
		// calculate it the inverse-probbility way
		double abs_eta = abs(NT_lep_eta_0);
		double pt = NT_lep_pt_0;
		double bin_x = (pt < muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax() - 1);
		double bin_y = (abs_eta < muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() - 1);
		no_mu_trig *= SingleMuon_data_bcdef_fraction * (1 - muon_effs_trg_BCDEF_histo->GetBinContent( muon_effs_trg_BCDEF_histo->FindBin(bin_x, bin_y) )) +
				SingleMuon_data_gh_fraction * (1 - muon_effs_trg_GH_histo->GetBinContent( muon_effs_trg_GH_histo->FindBin(bin_x, bin_y) ));
		//mu_trig_weight = 1 - no_trig; // so for 1 muon it will = to the SF, for 2 there will be a mix
		//fill_1d(string("weight_trigger_no_muon"),  200, 0., 1.1,   no_mu_trig, 1);

		weight_muon_trig = 1 - no_mu_trig;
		}
		weight *= weight_muon_trig;


		/* b-tagging SF
		 * for each jet check if it is tagged
		 * apply SF accordingly
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

		//FLoat_t b_discriminator = NT_jet_b_discr_0;
		//bool hasCSVtag(b_discriminator > b_tag_WP);
		if (NT_jet_pt_0 > 0) // jet exists
			weight_btag_sf *= calc_btag_sf_weight(btagCal, NT_jet_b_discr_0 > b_tag_WP, (int) NT_jet_hadronFlavour_0, NT_jet_pt_0, NT_jet_eta_0);
		if (NT_jet_pt_1 > 0) // jet exists
			weight_btag_sf *= calc_btag_sf_weight(btagCal, NT_jet_b_discr_1 > b_tag_WP, (int) NT_jet_hadronFlavour_1, NT_jet_pt_1, NT_jet_eta_1);
		if (NT_jet_pt_2 > 0) // jet exists
			weight_btag_sf *= calc_btag_sf_weight(btagCal, NT_jet_b_discr_2 > b_tag_WP, (int) NT_jet_hadronFlavour_2, NT_jet_pt_2, NT_jet_eta_2);
		if (NT_jet_pt_3 > 0) // jet exists
			weight_btag_sf *= calc_btag_sf_weight(btagCal, NT_jet_b_discr_3 > b_tag_WP, (int) NT_jet_hadronFlavour_3, NT_jet_pt_3, NT_jet_eta_3);
		if (NT_jet_pt_4 > 0) // jet exists
			weight_btag_sf *= calc_btag_sf_weight(btagCal, NT_jet_b_discr_4 > b_tag_WP, (int) NT_jet_hadronFlavour_4, NT_jet_pt_4, NT_jet_eta_4);

		weight *= weight_btag_sf;
		}


	// and pile-up:
	weight_without_pu = weight;
	if (is_MC)
		weight *= pileup_ratio[(int)NT_nvtx_gen];

	// select event categories
	const double lj_dist = 800;
	if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets == 3)) event_category = 0;
	else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets == 4)) event_category = 1;
	else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets >= 5)) event_category = 2;
	//else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets > 5) ) event_category = 3;
	else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets == 3)) event_category = 3;
	else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets == 4)) event_category = 4;
	else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets >= 5)) event_category = 5;
	//else if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets > 5) ) event_category = 7;

	// 1 inclusive cathegory for the peak
	else if (NT_lj_peak_distance < lj_dist && NT_nbjets >= 1 && NT_njets >= 3) event_category = 6;
	/*
	else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets == 3)) event_category = 6;
	else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets == 4)) event_category = 7;
	else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets >= 5)) event_category = 8;
	//else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets > 5) ) event_category =11;
	else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets == 4)) event_category = 9;
	else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets >= 5)) event_category =10;
	//else if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets > 5) ) event_category =14;
	*/
	else continue;

	// if mc channel is not present in distrs -- set it up
	if (event_yields_in_cathegories.find(mc_decay) == event_yields_in_cathegories.end())
		setup_distrs_for_mc_channel(event_yields_in_cathegories, distrs_in_cathegories, cathegories, mc_decay);

	// the inclusive distr will go to 1 event cathegory
	distrs_in_cathegories[mc_decay][cathegories[6]]["lj_peak_distance_inclusive"]->Fill(NT_lj_peak_distance, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["lj_peak_distance"]->Fill(NT_lj_peak_distance, weight);

	distrs_in_cathegories[mc_decay][cathegories[event_category]]["tau_energy_ratio"]->Fill(NT_tau_hcalEnergyLeadChargedHadrCand / NT_tau_hcalEnergy, weight);
	//distrs_in_cathegories[mc_decay][cathegories[event_category]]["met_transverse_mass"]->Fill(NT_met_corrected->Mt(), weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["met_pt"] ->Fill(NT_met_corrected, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["tau_pt"] ->Fill(NT_tau_pt_0, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["lep_pt"] ->Fill(NT_lep_pt_0, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["jet_pt"] ->Fill(NT_jet_pt_0, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["tau_eta"] ->Fill(NT_tau_eta_0, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["lep_eta"] ->Fill(NT_lep_eta_0, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["jet_eta"] ->Fill(NT_jet_eta_0, weight);
	//distrs_in_cathegories[mc_decay][cathegories[event_category]]["bjet_pt"]->Fill(NT_tau0_p4->pt(), weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["nvtx"]->Fill(NT_nvtx, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["rho"] ->Fill(NT_fixedGridRhoFastjetAll, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["nvtx_raw"]->Fill(NT_nvtx, weight_without_pu);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["rho_raw"] ->Fill(NT_fixedGridRhoFastjetAll, weight_without_pu);

	distrs_in_cathegories[mc_decay][cathegories[event_category]]["weight_muon_sfs"]->Fill(weight_muon_sfs, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["weight_muon_trig"]->Fill(weight_muon_trig, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["weight_btag_sf"] ->Fill(weight_btag_sf, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["weight_without_pu"]->Fill(weight_without_pu, weight);
	distrs_in_cathegories[mc_decay][cathegories[event_category]]["weight_final"]     ->Fill(weight, weight);

	event_yields_in_cathegories[mc_decay]->Fill(event_category, weight);
	}

cout << "done " << j << " events from " << n_events << " max " << NT_output_ttree->GetEntries() << endl;


// loop through all MC channels in the sample
for (auto const& imap: event_yields_in_cathegories)
	{
	TString mc_decay = imap.first;
	cout << "for " << mc_decay << " output" << endl;
	// yields in event cathegories
	TString yields_file = dir + "/jobsums/yields/" + dtag + mc_decay + ".root";
	TFile* out_f = TFile::Open(yields_file, "RECREATE");

	if (be_verbose) cout << "opened output file" << endl;
	eventflow->Write();
	imap.second->SetName("event_yields");
	imap.second->Write();
	//event_yields_in_cathegories[mc_decay]->Write();

	if (be_verbose) cout << "wrote objects" << endl;

	out_f->Write("event_yields");
	out_f->Close();

	// distrs
	for (int i=0; i<cathegories.size(); i++)
		{
		TString distr_file = dir + "/jobsums/cathegories/" + cathegories[i] + "/" + dtag + mc_decay + ".root";
		out_f = TFile::Open(distr_file, "RECREATE");
		eventflow->Write();
		//for (std::map<TString, TH1D*>::iterator it = (distrs_in_cathegories[cathegories[i]]).begin(); it != (distrs_in_cathegories[cathegories[i]]).end(); ++it)
		for (auto const& distr: distrs_in_cathegories[mc_decay][cathegories[i]])
			{
			TString controlpoint_name = distr.first;
			distr.second->SetName(controlpoint_name.Data());
			distr.second->Write();
			//distrs_in_cathegories[cathegories[i]][controlpoint_name]->SetName(controlpoint_name.Data());
			//distrs_in_cathegories[cathegories[i]][controlpoint_name]->Write();
			out_f->Write(controlpoint_name.Data());
			}
		out_f->Close();
		}
	}

return 0;
}

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

#include "UserCode/ttbar-leptons-80X/interface/wrapper.h"

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

double pu_weight(Int_t nvtx) {
   return pileup_ratio[nvtx];
}

namespace MyNamespace {
Double_t XReturner( Double_t x ) { return x; }
};


// need to open histos in this environment

//TString muon_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/muon-effs/";
//gSystem->ExpandPathName(muon_effs_dirname    );
TString muon_effs_dirname = "analysis/muon-effs/";

TFile * muon_effs_tracking_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_tracking_more_BCDEF_fits.root").Data() );
TFile * muon_effs_tracking_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_tracking_more_GH_fits.root").Data() );
TGraphAsymmErrors* muon_effs_tracking_BCDEF_graph = (TGraphAsymmErrors*) muon_effs_tracking_BCDEF_file->Get("ratio_eff_aeta_dr030e030_corr");
TGraphAsymmErrors* muon_effs_tracking_GH_graph    = (TGraphAsymmErrors*) muon_effs_tracking_GH_file->Get("ratio_eff_aeta_dr030e030_corr");
//cout << "Y tracking (reconstruction)" << endl;

TFile* muon_effs_id_BCDEF_file = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_BCDEF.root").Data() );
TFile* muon_effs_id_GH_file    = TFile::Open((muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_GH.root").Data() );
TH2D* muon_effs_id_BCDEF_histo = (TH2D*) ((TDirectoryFile*) muon_effs_id_BCDEF_file->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get("pt_abseta_ratio");
TH2D* muon_effs_id_GH_histo    = (TH2D*) ((TDirectoryFile*) muon_effs_id_GH_file   ->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get("pt_abseta_ratio");
//cout << "Y id" << endl;

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
//cout << "Y trigger" << endl;


double lepton_muon_SF(Float_t abs_eta, Float_t pt, double SingleMuon_data_bcdef_fraction, double SingleMuon_data_gh_fraction)
	{
	//weight *= 0.98; // average mu trig SF
	//weight *= 0.97; // average mu ID
	double weight_muon_sfs = 1;

	// muon ID SFs:
	//double weight_muon_sfs = 1;
	double bcdef_weight = 1;
	double gh_weight = 1;
	//double SingleMuon_data_bcdef_fraction = 19716.274 / (19716.274 + 15931.028);
	//double SingleMuon_data_gh_fraction    = 15931.028 / (19716.274 + 15931.028);

	//double abs_eta = abs(lep0_eta);
	//double pt = lep0_pt;

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

	//weight *= weight_muon_sfs;
	return weight_muon_sfs;
	}

// MUON Trigger
//double lepton_muon_trigger_SF ()
double lepton_muon_trig_SF(Float_t abs_eta, Float_t pt, double SingleMuon_data_bcdef_fraction, double SingleMuon_data_gh_fraction)
	{
	double no_mu_trig = 1;
	//double mu_trig_weight = 1;
	// calculate it the inverse-probbility way
	//double abs_eta = abs(NT_lep_eta_0);
	//double pt = NT_lep_pt_0;
	double bin_x = (pt < muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_trg_BCDEF_histo->GetXaxis()->GetXmax() - 1);
	double bin_y = (abs_eta < muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_trg_BCDEF_histo->GetYaxis()->GetXmax() - 1);
	no_mu_trig *= SingleMuon_data_bcdef_fraction * (1 - muon_effs_trg_BCDEF_histo->GetBinContent( muon_effs_trg_BCDEF_histo->FindBin(bin_x, bin_y) )) +
			SingleMuon_data_gh_fraction * (1 - muon_effs_trg_GH_histo->GetBinContent( muon_effs_trg_GH_histo->FindBin(bin_x, bin_y) ));
	//mu_trig_weight = 1 - no_trig; // so for 1 muon it will = to the SF, for 2 there will be a mix
	//fill_1d(string("weight_trigger_no_muon"),  200, 0., 1.1,   no_mu_trig, 1);

	//weight_muon_trig = 1 - no_mu_trig;
	return 1 - no_mu_trig;
	}
	//weight *= weight_muon_trig;


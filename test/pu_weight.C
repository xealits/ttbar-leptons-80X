#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <Math/VectorUtil.h>

#include "UserCode/ttbar-leptons-80X/interface/wrapper.h"

// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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
double lepton_muon_trigger_SF(Float_t abs_eta, Float_t pt, double SingleMuon_data_bcdef_fraction, double SingleMuon_data_gh_fraction)
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


/* ---------------------------------------------------
 * now, electrons have
 *      track(reconstruction) efficiency, which is recommended per eta of muon now (however there should be something about N vertices too..
 *      and ID sf
 *      also trigger
 *
 * the trig eff for dilepton case is: apply negative of it for both leptons
 */
//cout << "unpacking electron eff SFs" << endl;
TString electron_effs_dirname = "analysis/electron-effs/";

TFile* electron_effs_tracking_all_file = TFile::Open((electron_effs_dirname + "/2016_Sept23_ElectronReconstructionSF_egammaEffi.txt_EGM2D.root").Data() );
TH2D* electron_effs_tracking_all_histo = (TH2D*) electron_effs_tracking_all_file->Get("EGamma_SF2D");
//cout << "Y tracking (reconstruction)" << endl;

// for the selected electrons, Tight ID
// not for Veto
TFile* electron_effs_id_all_file = TFile::Open((electron_effs_dirname + "/2016_Sept23_ElectronID_TightCutBased_egammaEffi.txt_EGM2D.root").Data() );
TH2D* electron_effs_id_all_histo = (TH2D*) electron_effs_id_all_file->Get("EGamma_SF2D");
//cout << "Y id" << endl;

//analysis/electron-effs/2016_03Feb_TriggerSF_Run2016All_v1.root
TFile* electron_effs_trg_all_file = TFile::Open((electron_effs_dirname + "/2016_03Feb_TriggerSF_Run2016All_v1.root").Data() );
TH2D* electron_effs_trg_all_histo = (TH2D*) electron_effs_trg_all_file->Get("Ele27_WPTight_Gsf");
//cout << "Y trigger" << endl;

// --- these SFs will be applied to the selected leptons independently


double lepton_electron_SF(Float_t eta, Float_t pt)
	{
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

	return weight_reco * weight_id;
	}


double lepton_electron_trigger_SF(Float_t eta, Float_t pt)
	{
	//double no_ele_trig = 1;
	// (calculate it the inverse-probbility way)
	// no, just return the SF, assume 1-lepton case
        //pat::Electron& el = selElectrons[i];
	//double eta = el.superCluster()->position().eta();
	//double pt = el.pt();
	// here X axis is pt, Y axis is eta (from -2.5 to 2.5)
	double bin_x = (pt  < electron_effs_trg_all_histo->GetXaxis()->GetXmax() ? pt  : electron_effs_trg_all_histo->GetXaxis()->GetXmax() - 1);
	double bin_y = (eta < electron_effs_trg_all_histo->GetYaxis()->GetXmax() ? eta : electron_effs_trg_all_histo->GetYaxis()->GetXmax() - 1);
	return electron_effs_trg_all_histo->GetBinContent( electron_effs_trg_all_histo->FindBin(bin_x, bin_y) );
	//el_trig_weight = 1 - no_trig; // so for 1 lepton it will = to the SF, for 2 there will be a mix
	//fill_1d(string("weight_trigger_no_electron"),  200, 0., 1.1,   no_ele_trig, 1);
	}



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
	return exp(a + b*x);
	}

double ttbar_pT_SF(double t_pt, double tbar_pt)
	{
	return sqrt(top_pT_SF(t_pt) * top_pT_SF(tbar_pt));
	}

double transverse_mass(LorentzVector v1, LorentzVector v2)
	{
	return sqrt(2*v1.pt()*v2.pt()*(1 - cos(v1.phi() - v2.phi())));
	}


/*
 * Fake Rates
 * TH-s and application
 */

/*
 * Usage:
 *
 * jetToTauFakeRate(tau_fake_rate1_jets_histo_q, tau_fake_rate1_taus_histo_q, jet.pt(), jet.eta(), jet_rad, debug)
 *
 */
double jetToTauFakeRate(TH3D * tau_fake_rate_jets_histo1, TH3D * tau_fake_rate_taus_histo1, Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
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

	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);

	/*
	if (debug)
		{
		cout << jet_pt << " " << jet_eta << " " << jet_radius << " : " << global_bin_id << " : ";
		cout << taus_rate1 << "/" << jets_rate1 << endl;
		}
	*/

	return fakerate;
	}


/*
 * dataDriven_tauFakeRates1 :          '${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/jobing/configs/jet_to_tau_fakerates1_tauMVA.root'
 * dataDriven_tauFakeRates2 :          '${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/jobing/configs/jet_to_tau_fakerates2_tauMVA.root'
 * dataDriven_tauFakeRates_dileptons : '${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/jobing/configs/jet_to_tau_fakerates_dileptons_tauMVA.root'
 */
TString fakerates_dirname = "jobing/configs/";
TFile * tau_fake_rate1_file = TFile::Open((fakerates_dirname + "jet_to_tau_fakerates1_tauMVA.root").Data());
TFile * tau_fake_rate2_file = TFile::Open((fakerates_dirname + "jet_to_tau_fakerates2_tauMVA.root").Data());
TFile * tau_fake_rate_file_dileptons = TFile::Open((fakerates_dirname + "jet_to_tau_fakerates_dileptons_tauMVA.root").Data());
//cout << "fake rate QCD file: " << tau_fake_rate1_file << endl;

TH3D * tau_fake_rate1_jets_histo_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_jets_distr_large_bins");
TH3D * tau_fake_rate1_taus_histo_q = (TH3D *) tau_fake_rate1_file->Get("HLTjet_qcd_tau_jets_distr_large_bins");

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

double jetToTauFakeRate_qcd(Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	Int_t global_bin_id = tau_fake_rate1_jets_histo_q->FindBin(jet_pt, jet_eta, jet_radius);
	Double_t jets_rate1 = tau_fake_rate1_jets_histo_q->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate1_taus_histo_q->GetBinContent(global_bin_id);
	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);
	return fakerate;
	}

double jetToTauFakeRate_wjets(Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	Int_t global_bin_id = tau_fake_rate2_jets_histo_w->FindBin(jet_pt, jet_eta, jet_radius);
	Double_t jets_rate1 = tau_fake_rate2_jets_histo_w->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate2_taus_histo_w->GetBinContent(global_bin_id);
	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);
	return fakerate;
	}

double jetToTauFakeRate_mumu(Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	Int_t global_bin_id = tau_fake_rate_jets_histo_mumu->FindBin(jet_pt, jet_eta, jet_radius);
	Double_t jets_rate1 = tau_fake_rate_jets_histo_mumu->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate_taus_histo_mumu->GetBinContent(global_bin_id);
	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);
	return fakerate;
	}

double jetToTauFakeRate_elmu(Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	Int_t global_bin_id = tau_fake_rate_jets_histo_elmu->FindBin(jet_pt, jet_eta, jet_radius);
	Double_t jets_rate1 = tau_fake_rate_jets_histo_elmu->GetBinContent(global_bin_id);
	Double_t taus_rate1 = tau_fake_rate_taus_histo_elmu->GetBinContent(global_bin_id);
	Double_t fakerate = (jets_rate1 < 1 ? 0 : taus_rate1/jets_rate1);
	return fakerate;
	}

Int_t elmu_FindBin(Double_t jet_pt, Double_t jet_eta, Double_t jet_radius, bool debug)
	{
	Int_t global_bin_id = tau_fake_rate_jets_histo_elmu->FindBin(jet_pt, jet_eta, jet_radius);
	return global_bin_id;
	}





/*
 * Likelihood distrs:
 *    njets with/without PU
 *    lj peak
 *
 */

int njets_over_pu(double jet0_pu_discr, double jet1_pu_discr, double jet2_pu_discr, double jet3_pu_discr, double jet4_pu_discr, double pu_threshold)
	{
	return (jet0_pu_discr > pu_threshold ?1:0) + (jet1_pu_discr > pu_threshold ?1:0) + (jet2_pu_discr > pu_threshold ?1:0) + (jet3_pu_discr > pu_threshold ?1:0) + (jet4_pu_discr > pu_threshold ?1:0);
	}


/*
 * bloody root: the TTree->Draw("parameterfile", "conditionfile") works in interpreter, but not in compiled script
closest_w_mass(jet0_id, jet0_b_discr, jet0_p4.X(), jet0_p4.Y(), jet0_p4.Z(), jet0_p4.T(), jet1_id, jet1_b_discr, jet1_p4.X(), jet1_p4.Y(), jet1_p4.Z(), jet1_p4.T(), jet2_id, jet2_b_discr, jet2_p4.X(), jet2_p4.Y(), jet2_p4.Z(), jet2_p4.T(), jet3_id, jet3_b_discr, jet3_p4.X(), jet3_p4.Y(), jet3_p4.Z(), jet3_p4.T(), jet4_id, jet4_b_discr, jet4_p4.X(), jet4_p4.Y(), jet4_p4.Z(), jet4_p4.T())

 */
double closest_w_mass(
	Int_t jet0_id, Float_t jet0_b_discr, Float_t jet0_X, Float_t jet0_Y, Float_t jet0_Z, Float_t jet0_T, 
	Int_t jet1_id, Float_t jet1_b_discr, Float_t jet1_X, Float_t jet1_Y, Float_t jet1_Z, Float_t jet1_T, 
	Int_t jet2_id, Float_t jet2_b_discr, Float_t jet2_X, Float_t jet2_Y, Float_t jet2_Z, Float_t jet2_T, 
	Int_t jet3_id, Float_t jet3_b_discr, Float_t jet3_X, Float_t jet3_Y, Float_t jet3_Z, Float_t jet3_T, 
	Int_t jet4_id, Float_t jet4_b_discr, Float_t jet4_X, Float_t jet4_Y, Float_t jet4_Z, Float_t jet4_T
	)
	{
	//return (jet0_pu_discr > pu_threshold ?1:0) + (jet1_pu_discr > pu_threshold ?1:0) + (jet2_pu_discr > pu_threshold ?1:0) + (jet3_pu_discr > pu_threshold ?1:0) + (jet4_pu_discr > pu_threshold ?1:0);
	std::vector<LorentzVector> light_jets, heavy_jets;

	Float_t b_threshold = 0.8484;

	if (jet0_id >=0)
		{
		if (jet0_b_discr < b_threshold) light_jets.push_back(LorentzVector(jet0_X, jet0_Y, jet0_Z, jet0_T));
		}
	if (jet1_id >=0)
		{
		if (jet1_b_discr < b_threshold) light_jets.push_back(LorentzVector(jet1_X, jet1_Y, jet1_Z, jet1_T));
		}
	if (jet2_id >=0)
		{
		if (jet2_b_discr < b_threshold) light_jets.push_back(LorentzVector(jet2_X, jet2_Y, jet2_Z, jet2_T));
		}
	if (jet3_id >=0)
		{
		if (jet3_b_discr < b_threshold) light_jets.push_back(LorentzVector(jet3_X, jet3_Y, jet3_Z, jet3_T));
		}
	if (jet4_id >=0)
		{
		if (jet4_b_discr < b_threshold) light_jets.push_back(LorentzVector(jet4_X, jet4_Y, jet4_Z, jet4_T));
		}

	// at 2 light jets
	if (!(light_jets.size()>1)) return 0;

	double closest_W_mass = 0;
	double smallest_W_mass_distance = 999999;

	LorentzVector W(0,0,0,0);
	for (unsigned int i=0; i<light_jets.size()-1; i++)
		for (unsigned int u=i+1; u<light_jets.size(); u++)
			{
			W = light_jets[i] + light_jets[u];
			double W_mass_dist = ( 80 - W.mass());
			if (W_mass_dist < smallest_W_mass_distance)
				{
				smallest_W_mass_distance = W_mass_dist;
				closest_W_mass = W.mass();
				}
			}

	if (closest_W_mass > 160) closest_W_mass = 160;
	return closest_W_mass;
	}

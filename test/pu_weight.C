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

const double b_Medium_WP = 0.8484;

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

double pileup_ratio_up[] = {0, 0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0 };

double pileup_ratio_down[] = {0, 0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0 };


double pu_weight(Int_t nvtx, int pu_shift) {
	switch (pu_shift) {
		case 0:
			return pileup_ratio[nvtx];
		case -1:
			return pileup_ratio_down[nvtx];
		case 1:
			return pileup_ratio_up[nvtx];
	}
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

/*
bool pass_jet_met_b_requirement_jes(int njets_theshold, double met_threshold, int nbjets_theshold, int sys_shift,
	int njets, int nbjets,
	// jets should pass pt 30, met -- 40, bjets -- medium WP
	Float_t jet0_b_discr, Float_t jet0_pt, Float_t jet0_phi, Float_t jet0_jes_correction_relShift,
	Float_t jet1_b_discr, Float_t jet1_pt, Float_t jet1_phi, Float_t jet1_jes_correction_relShift,
	Float_t jet2_b_discr, Float_t jet2_pt, Float_t jet2_phi, Float_t jet2_jes_correction_relShift,
	Float_t jet3_b_discr, Float_t jet3_pt, Float_t jet3_phi, Float_t jet3_jes_correction_relShift,
	Float_t jet4_b_discr, Float_t jet4_pt, Float_t jet4_phi, Float_t jet4_jes_correction_relShift,
	Float_t met_pt, Float_t met_phi)
{
TVector3 met(0,0,0);
met.SetPtEtaPhi(met_pt, 0, met_phi);
TVector3 jet0(0,0,0);
TVector3 jet1(0,0,0);
TVector3 jet2(0,0,0);
TVector3 jet3(0,0,0);
TVector3 jet4(0,0,0);
jet0.SetPtEtaPhi(jet0_pt,0,jet0_phi);
jet1.SetPtEtaPhi(jet1_pt,0,jet1_phi);
jet2.SetPtEtaPhi(jet2_pt,0,jet2_phi);
jet3.SetPtEtaPhi(jet3_pt,0,jet3_phi);
jet4.SetPtEtaPhi(jet4_pt,0,jet4_phi);

//int n_top_jets = (jet0.pt() > 30) + (jet1.pt() > 30) + (jet2.pt() > 30) + (jet3.pt() > 30) + (jet4.pt() > 30);
// they all are > 30 now

TVector3 jet_sum = jet0 + jet1 + jet2 + jet3 + jet4;
TVector3 jet_sum_new(0,0,0);

if (sys_shift != 0)
	{
	switch(sys_shift) {
		case 1: // UP
			jet0 *= (1+jet0_jes_correction_relShift);
			jet1 *= (1+jet1_jes_correction_relShift);
			jet2 *= (1+jet2_jes_correction_relShift);
			jet3 *= (1+jet3_jes_correction_relShift);
			jet4 *= (1+jet4_jes_correction_relShift);
			break;
		case 2:
			jet0 *= (1-jet0_jes_correction_relShift);
			jet1 *= (1-jet1_jes_correction_relShift);
			jet2 *= (1-jet2_jes_correction_relShift);
			jet3 *= (1-jet3_jes_correction_relShift);
			jet4 *= (1-jet4_jes_correction_relShift);
			break;
	}

	// calculate new njets
	int n_top_jets_new = (jet0.pt() > 30) + (jet1.pt() > 30) + (jet2.pt() > 30) + (jet3.pt() > 30) + (jet4.pt() > 30);

	jet_sum_new = jet0 + jet1 + jet2 + jet3 + jet4;
	// propagate to met
	met += jet_sum - jet_sum_new;
	}

return njets > njets_theshold && nbjets > nbjets_theshold && met.pt() > met_threshold;
}
*/


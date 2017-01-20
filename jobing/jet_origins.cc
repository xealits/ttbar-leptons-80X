#include <iostream>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
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
#include <vector>

using namespace std;

/*
int stacked_histo_distr (int argc, char *argv[])

where really arguments are:
(TString distr, TString dir, TString dtags)

convert char* dtags[] to TStrings,
have dtag->xsec map around,
open dir/dtag.root files,
get the distr and weightflows (for ratios, so doesn't matter which) from there,
make stack of ratio-weighted MC-s (so, separate MC and Data)
and overlay it with data histo
+ error bars everywhere

*/

double W_lep_br = 0.108;
double W_qar_br = 0.676;

double W_lep_br2 = 0.108*0.108;
double W_qar_br2 = 0.676*0.676;

std::map<TString, double> xsecs = {
{"Data13TeV_SingleElectron2016D_PromptRecoV2_", 1},
{"Data13TeV_SingleMuon2016D_PromptRecoV2_", 1},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elmubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_eltaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mumubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mutaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qelbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qmubar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qqbar"           , 831.76 * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qtaubar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauelbar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_taumubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauqbar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tautaubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elmubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_eltaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mumubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mutaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qelbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qmubar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qqbar"         , 831.76 * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qtaubar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauelbar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_taumubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauqbar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tautaubar"     , 831.76 * W_lep_br2 / 2},

{ "MC2016_reHLT_TTJets_powheg_scaleup_elelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_elmubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_elqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_eltaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_muelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_mumubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_muqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_mutaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_qelbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_qmubar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_qqbar"           , 831.76 * W_qar_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_qtaubar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_tauelbar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_taumubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_tauqbar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaleup_tautaubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_elelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_elmubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_elqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_eltaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_muelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_mumubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_muqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_mutaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_qelbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_qmubar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_qqbar"         , 831.76 * W_qar_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_qtaubar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_tauelbar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_taumubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_tauqbar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_scaledown_tautaubar"     , 831.76 * W_lep_br2 / 2},

{"MC2016_noHLT_WJets_amcatnlo",  61526.7 },
{"MC2016_noHLT_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W1Jets_madgraph", 9493},
{"MC2016_noHLT_W2Jets_madgraph", 3120},
{"MC2016_noHLT_W3Jets_madgrapg", 942.3},
{"MC2016_noHLT_W4Jets_madgraph", 524.2},
{"MC2016_reHLT_WJets_amcatnlo",  61526.7 },
{"MC2016_reHLT_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_reHLT_DYJetsToLL_50toInf_amcatnlo", 6025.2},
{"MC2016_noHLT_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo", 6025.2},
{"MC2016_noHLT_WW", 113.89},
{"MC2016_noHLT_WZ", 47.13},
{"MC2016_noHLT_ZZ", 16.52},
{"MC2016_noHLT_SingleT_tW_5FS_powheg",    35.6},
{"MC2016_noHLT_SingleTbar_tW_5FS_powheg", 35.6},
{"MC2016_noHLT_schannel_4FS_leptonicDecays_amcatnlo", 3.36},
{"MC2016_noHLT_tchannel_antitop_4f_leptonicDecays_powheg", 70.69/2},
{"MC2016_noHLT_tchannel_top_4f_leptonicDecays_powheg", 70.69/2},
{"MC2016_noHLT_QCD_HT-100-200",  27540000},
{"MC2016_noHLT_QCD_HT-200-300",  1717000},
{"MC2016_noHLT_QCD_HT-300-500",  351300},
{"MC2016_noHLT_QCD_HT-500-700",  31630},
{"MC2016_noHLT_QCD_HT-700-1000",  6802},
{"MC2016_noHLT_QCD_HT-1000-1500",  1206},
{"MC2016_noHLT_QCD_HT-1500-2000",  120.4},
{"MC2016_noHLT_QCD_HT-2000-Inf",  25.25}
};


/*
std::map<TString, EColor> dtag_colours = {
{"Data13TeV_SingleElectron2016D_PromptRecoV2_", kAzure},
{"Data13TeV_SingleMuon2016D_PromptRecoV2_",     kAzure},
{"MC2016_noHLT_DYJetsToLL_10to50_amcatnlo.root",  kGray},
{"MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo.root", kGray},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elelbar"         , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elmubar"         , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elqbar"          , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_eltaubar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muelbar"         , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mumubar"         , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muqbar"          , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mutaubar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qelbar"          , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qmubar"          , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qqbar"           , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qtaubar"         , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauelbar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_taumubar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauqbar"         , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tautaubar"       , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elelbar"       , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elmubar"       , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elqbar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_eltaubar"      , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muelbar"       , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mumubar"       , (kGreen-3)},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muqbar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mutaubar"      , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qelbar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qmubar"        , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qqbar"         , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qtaubar"       , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauelbar"      , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_taumubar"      , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauqbar"       , kGreen},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tautaubar"     , kGreen},
{"MC2016_noHLT_W0Jets_amcatnlo.root", kOrange},
{"MC2016_noHLT_W1Jets_madgraph.root", kOrange},
{"MC2016_noHLT_W2Jets_madgraph.root", kOrange},
{"MC2016_noHLT_W3Jets_madgrapg.root", kOrange},
{"MC2016_noHLT_W4Jets_madgraph.root", kOrange},
{"MC2016_noHLT_WW.root", kCyan},
{"MC2016_noHLT_WZ.root", kCyan},
{"MC2016_noHLT_ZZ.root", kCyan},
{"MC2016_noHLT_SingleT_tW_5FS_powheg.root", kAzure},
{"MC2016_noHLT_SingleTbar_tW_5FS_powheg.root", kAzure},
{"MC2016_noHLT_schannel_4FS_leptonicDecays_amcatnlo.root", kAzure},
{"MC2016_noHLT_tchannel_antitop_4f_leptonicDecays_powheg.root", kAzure},
{"MC2016_noHLT_tchannel_top_4f_leptonicDecays_powheg.root", kAzure}
};
*/

/*
// nick and colour
std::pair<TString, Color_t> dtag_nick_colour(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("dyjets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("wjets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("singletop", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qqbar")) return std::make_pair("tt_jj", kGreen+4);
		else if (dtag.Contains("elqbar") || dtag.Contains("qelbar") ||dtag.Contains("muqbar") || dtag.Contains("qmubar") || dtag.Contains("tauqbar") || dtag.Contains("qtaubar")) return std::make_pair("tt_lj", kGreen+3);
		else if (dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt_em", kGreen-9);
		else if (dtag.Contains("elelbar")) return std::make_pair("tt_ee", kAzure-9);
		else if (dtag.Contains("mumubar")) return std::make_pair("tt_mm", kYellow-7);
		else if (dtag.Contains("eltaubar") || dtag.Contains("tauelbar")) return std::make_pair("tt_et", kOrange+4);
		else if (dtag.Contains("mutaubar") || dtag.Contains("taumubar")) return std::make_pair("tt_mt", kOrange+1);
		else return std::make_pair("tt_other", kYellow+1);
		}
	else return std::make_pair("other", kBlack);

	}
*/

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 4)
	{
	std::cout << "Usage : " << argv[0] << " job_dir dtag origin_channel [origins_channels]" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString job_dir(argv[1]);
TString dtag(argv[2]);

//std::vector <TString> channels = {"HLTjet_qcd_jet_origins", "HLTmu_qcd_jet_origins", "HLTjetmu_qcd_jet_origins",
//	"HLTjet_wjets_jet_origins", "HLTmu_wjets_jet_origins", "HLTjetmu_wjets_jet_origins" };

//std::vector <TString> channels = {"HLTjet_qcd_jet_partonFlavour", "HLTmu_qcd_jet_partonFlavour", "HLTjetmu_qcd_jet_partonFlavour",
//	"HLTjet_wjets_jet_partonFlavour", "HLTmu_wjets_jet_partonFlavour", "HLTjetmu_wjets_jet_partonFlavour" };
// in ttbar dieptons:
//std::vector <TString> channels = {"dilep_passjets_jet_partonFlavour", "dilep_passjetsNbtag_jet_partonFlavour", "elel_passjets_jet_partonFlavour"};

std::vector<TString> channels;
for (int i=3; i<argc; i++)
	{
	TString a_channel(argv[i]);
	channels.push_back(a_channel);
	}

//cout << job_dir  << endl;
//cout << dtag << endl;

TFile * file = TFile::Open(job_dir + "/" + dtag + ".root");
//TH1D * weightflow = (TH1D*) file->Get("weightflow");
TH1D * weightflow;
// for ttbar dileptons
if (file->GetListOfKeys()->Contains("weightflow_elel"))
	weightflow = (TH1D*) file->Get("weightflow_elel");
else if (file->GetListOfKeys()->Contains("weightflow"))
	weightflow = (TH1D*) file->Get("weightflow");
else
	{
	cerr << "no weightflow distro" << endl;
	return 1;
	}


double xsec = xsecs[dtag];
Double_t ratio = xsec / weightflow->GetBinContent(4);

//cout << "dtag " << dtag << " ratio = " << ratio << endl;


cout << "dtag,ratio";
for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch)
	{
	/*
	if (!file->GetListOfKeys()->Contains(*ch))
		{
		cerr << "no " << *ch << endl;
		continue;
		}
	*/
	cout << "," << *ch << "_o," << *ch << "_t," << *ch << "_b," << *ch << "_g," << *ch << "_a";
	}
cout << endl;

cout << dtag << "," << ratio;

for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch)
	{
	if (!file->GetListOfKeys()->Contains(*ch))
		{
		//cerr << "no " << *ch << endl;
		cout << ",NA,NA,NA,NA,NA";
		continue;
		}

	TH1D * jet_origins = (TH1D*) file->Get(*ch);

	//Double_t gluons_to_all = jet_origins->GetBinContent(22) / jet_origins->Integral();
	cout << "," << jet_origins->GetBinContent(1) << "," << jet_origins->GetBinContent(16) << "," << jet_origins->GetBinContent(6) << "," << jet_origins->GetBinContent(22) << "," << jet_origins->Integral();
	}
cout << endl;

/*
TH1D * jet_origins_HLTjet_QCD = (TH1D*) file->Get("HLTjet_qcd_jet_origins");
// other origins:

cout << "HLTjet_qcd_jet_origins\n";
cout << "gluons\tall\n";
cout << jet_origins_HLTjet_QCD->GetBinContent(22)*ratio << "\t" << jet_origins_HLTjet_QCD->Integral()*ratio << endl;
*/

return 0;
}


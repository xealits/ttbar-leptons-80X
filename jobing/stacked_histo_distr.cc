#include <iostream>

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
{"MC2016_noHLT_DYJetsToLL_10to50_amcatnlo.root", 18610},
{"MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo.root", 6025.2},
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
{"MC2016_noHLT_W0Jets_amcatnlo.root", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W1Jets_madgraph.root", 9493},
{"MC2016_noHLT_W2Jets_madgraph.root", 3120},
{"MC2016_noHLT_W3Jets_madgrapg.root", 942.3},
{"MC2016_noHLT_W4Jets_madgraph.root", 524.2},
{"MC2016_noHLT_WW.root", 113.89},
{"MC2016_noHLT_WZ.root", 47.13},
{"MC2016_noHLT_ZZ.root", 16.52},
{"MC2016_noHLT_SingleT_tW_5FS_powheg.root",    35.6},
{"MC2016_noHLT_SingleTbar_tW_5FS_powheg.root", 35.6},
{"MC2016_noHLT_schannel_4FS_leptonicDecays_amcatnlo.root", 3.36},
{"MC2016_noHLT_tchannel_antitop_4f_leptonicDecays_powheg.root", 70.69/2},
{"MC2016_noHLT_tchannel_top_4f_leptonicDecays_powheg.root", 70.69/2}
};

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 4)
	{
	std::cout << "Usage : " << argv[0] << " distr, dir, dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString distr(argv[1]);
TString dir(argv[2]);
TString dtag1(argv[3]);
cout << distr << endl;
cout << dir   << endl;
cout << dtag1 << endl;


std::vector < TString > dtags;
std::vector < TFile * > files;
std::vector < TH1D * > histos;
//vector<int> dtags;
//dtags.reserve();

for (int i = 3; i<argc; i++)
	{
	TString dtag(argv[i]);
	dtags.push_back(dtag);
	files.push_back(TFile::Open(dir + "/" + dtag + ".root"));
	
	//dtags.push_back(5);
	}

for (int i = 0; i<dtags.size(); i++)
	{
	cout << dtags[i] << endl;
	histos.push_back((TH1D*) files[i]->Get(distr));
	histos[i]->Print();
	}

for(std::map<TString, double>::iterator it = xsecs.begin(); it != xsecs.end(); ++it)
	{
	TString dtag = it->first;
	double xsec  = it->second;
	cout << "For dtag " << dtag << " xsec " << xsec << "\n";
	}

/*
TH1D *h = (TH1D*) file->Get(distr);

Int_t size_x = h->GetNbinsX();

if (!print_header)
	{
	cout << dtag;
	for (int x=1; x<size_x; x++)
		{
		//double bin_center = h->GetXaxis()->GetBinCenter(x);
		double global_bin = h->GetBin(x);
		cout << "," << h->GetBinContent(global_bin);
		}
	cout << "\n";
	}
else
	{
	cout << "dtag";
	for (int x=1; x<size_x; x++)
		{
		double bin_center = h->GetXaxis()->GetBinCenter(x);
		cout << "," << bin_center;
		}
	cout << "\n";
	}
*/
return 0;
}


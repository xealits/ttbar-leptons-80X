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
{"MC2016_noHLT_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo", 6025.2},
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
{"MC2016_noHLT_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W1Jets_madgraph", 9493},
{"MC2016_noHLT_W2Jets_madgraph", 3120},
{"MC2016_noHLT_W3Jets_madgrapg", 942.3},
{"MC2016_noHLT_W4Jets_madgraph", 524.2},
{"MC2016_noHLT_WW", 113.89},
{"MC2016_noHLT_WZ", 47.13},
{"MC2016_noHLT_ZZ", 16.52},
{"MC2016_noHLT_SingleT_tW_5FS_powheg",    35.6},
{"MC2016_noHLT_SingleTbar_tW_5FS_powheg", 35.6},
{"MC2016_noHLT_schannel_4FS_leptonicDecays_amcatnlo", 3.36},
{"MC2016_noHLT_tchannel_antitop_4f_leptonicDecays_powheg", 70.69/2},
{"MC2016_noHLT_tchannel_top_4f_leptonicDecays_powheg", 70.69/2}
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
		else if (dtag.Contains("eltaubar") || dtag.Contains("tauelbar")) return std::make_pair("tt_et", kOrange+2);
		else if (dtag.Contains("mutaubar") || dtag.Contains("taumubar")) return std::make_pair("tt_mt", kOrange+1);
		else return std::make_pair("tt_other", kYellow+1);
		}
	else return std::make_pair("other", kBlack);

	}

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 5)
	{
	std::cout << "Usage : " << argv[0] << " lumi, distr, dir, dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

double lumi = atof(argv[1]);
TString distr(argv[2]);
TString dir(argv[3]);
TString dtag1(argv[4]);

cout << lumi  << endl;
cout << distr << endl;
cout << dir   << endl;
cout << dtag1 << endl;

for(std::map<TString, double>::iterator it = xsecs.begin(); it != xsecs.end(); ++it)
	{
	TString dtag = it->first;
	double xsec  = it->second;
	cout << "For dtag " << dtag << " xsec " << xsec << "\n";
	}


std::vector < TString > dtags;
std::vector < TFile * > files;
std::vector < TH1D * > histos;
std::vector < TH1D * > weightflows;
// nick->summed histo
std::map<TString, TH1D *> nicknamed_mc_histos;
//vector<int> dtags;
//dtags.reserve();

// make stack of MC, scaling according to ratio = lumi * xsec / weightflow4 (bin5?)
// also nickname the MC....
// per-dtag for now..

//THStack *hs = new THStack("hs","Stacked 1D histograms");
THStack *hs      = new THStack("hs", "");
TH1D    *hs_data = NULL;

TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);

//TLegend *leg = new TLegend(0.845, 0.2, 0.99, 0.99);
//leg = new TLegend(0.845, 0.2, 0.99, 0.99);
TLegend* leg = new TLegend(0.845, 0.5, 0.99, 0.99);

for (int i = 4; i<argc; i++)
	{
	TString dtag(argv[i]);
	cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	files.push_back(TFile::Open(dir + "/" + dtag + ".root"));
	//dtags.push_back(5);

	if (!files.back()->GetListOfKeys()->Contains(distr))
		{
		cout << "no " << distr << endl;
		continue;
		}

	weightflows.push_back((TH1D*) files.back()->Get("weightflow_el"));
	weightflows.back()->Print();

	histos.push_back((TH1D*) files.back()->Get(distr));
	histos.back()->Print();
	if (dtag.Contains("Data"))
		{
		cout << "summing data-stack" << endl;
		histos.back()->SetMarkerStyle(9);
		//histos.back()->SetFillColor(kRed + i);

		if (hs_data == NULL)
			{
			cout << "creating data histo" << endl;
			hs_data = (TH1D*) histos.back()->Clone();
			}
		else
			{
			cout << "add histo to data histo" << endl;
			hs_data->Add(histos.back());
			}
		}
	else
		{
		Double_t ratio = lumi * xsecs[dtag] / weightflows.back()->GetBinContent(5);
		histos.back()->Scale(ratio);
		cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		histos.back()->Print();
		//histos.back()->SetFillColor(kRed );

		std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
		TString nick = nick_colour.first;
		Color_t col = nick_colour.second;
		histos.back()->SetFillColor( col );

		histos.back()->SetMarkerStyle(20);
		histos.back()->SetLineStyle(0);
		histos.back()->SetMarkerColor(i);

		//std::map<TString, TH1D *> nicknamed_mc_histos;
		if (nicknamed_mc_histos.find(nick) == nicknamed_mc_histos.end())
			nicknamed_mc_histos[nick] = (TH1D*) histos.back()->Clone();
			//nicknamed_mc_histos[nick] = new TH1D(nick, "");
		else
			nicknamed_mc_histos[nick]->Add(histos.back());
		//hs->Add(histos.back(), "HIST");
		}
	}

//std::map<TString, TH1D *> nicknamed_mc_histos;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;
	cout << "adding to mc stack: " << nick << endl;
	distr->Print();

	hs->Add(distr, "HIST");

	//TLegendEntry *entry=leg->AddEntry("NULL","","h");
	//entry=leg->AddEntry("singlemu_altstep4rho3","Single top","F");
	leg->AddEntry(distr, nick, "F");

	//distr->SetName(controlpoint_name.c_str());
	//distr->Write();
	//out_f->Write(controlpoint_name.c_str());
	//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
	}

//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

cout << "setting title" << endl;

//hs->GetXaxis()->SetTitle(distr);
//cst->SetXaxisTile(distr);
//hs_data->GetXaxis()->SetTitle(distr);

/*
TIter next(gDirectory->GetList());
TObject *obj;
while ((obj=next())) {
	if (obj->InheritsFrom("TH1")) {
		TH1 *h = (TH1*)obj;
		h->GetXaxis()->SetTitle(distr);
		}
	}
*/


TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("x");
//cst->Update();
//cst->Modified();

cout << "drawing" << endl;

hs->Draw();
hs_data->Draw("e p same");

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


hs->GetXaxis()->SetTitle(distr);
hs_data->SetXTitle(distr);

cst->Modified();

cst->SaveAs("./teststack.png");

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


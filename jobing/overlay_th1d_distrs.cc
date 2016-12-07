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

#include "dtag_xsecs.h"

#define INPUT_DTAGS_START 6

using namespace std;



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

//histo_project3toN_ratio_distr.cc
//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 4)
	{
	//std::cout << "Usage : " << argv[0] << " lumi distr projection rebin_factor dir dtags" << std::endl;
	//std::cout << "Usage : " << argv[0] << " x|y|z proj_name distr1 distr2 rebin_factor x_axis_min_range x_axis_max_range filename" << std::endl;
	std::cout << "Usage : " << argv[0] << "outfilename name filename distrname [name filename distrname]*" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString outfilename(argv[1]);

TString name(argv[2]);
TString filename(argv[3]);
TString distrname(argv[4]);

cout << argc << endl;
cout << filename    << endl;
cout << distrname   << endl;

TFile * file = TFile::Open(filename);
TH1D * h = (TH1D*) file->Get(distrname);

TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);
TLegend* leg = new TLegend(0.845, 0.5, 0.99, 0.99);

h->Draw("ep");       // Draw the ratio plot
leg->AddEntry(h, name, "F");

for (int i=5; i<argc; i+=3)
	{
	TString name(argv[i]);
	TString filename(argv[i+1]);
	TString distrname(argv[i+2]);

	TFile * file = TFile::Open(filename);
	TH1D * h = (TH1D*) file->Get(distrname);

	leg->AddEntry(h, name, "F");

	h->SetLineColor(kOrange + i - 4);
	h->SetMarkerColor(kOrange + i - 4);
	h->Draw("epsame");       // Draw the ratio plot
	}

leg->Draw();

cst->SetLogy();

cst->Update();

cst->Modified();

cst->SaveAs(outfilename);

return 0;
}


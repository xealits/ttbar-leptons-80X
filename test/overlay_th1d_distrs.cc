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
TLegend* leg = new TLegend(0.845, 0.8, 0.99, 0.99);

h->Draw("ep");       // Draw the ratio plot
h->Print();
leg->AddEntry(h, name, "F");

for (int i=5; i<argc; i+=3)
	{
	TString name(argv[i]);
	TString filename(argv[i+1]);
	TString distrname(argv[i+2]);
	cout << filename << '\t' << distrname << endl;

	TFile * file = TFile::Open(filename);
	TH1D * h = (TH1D*) file->Get(distrname);

	h->Print();

	leg->AddEntry(h, name, "F");

	h->SetLineColor(kOrange + i - 4);
	h->SetMarkerColor(kOrange + i - 4);
	h->Draw("epsame");       // Draw the ratio plot
	}

leg->Draw();

//cst->SetLogy();

//cst->Update();

//cst->Modified();

cst->SaveAs(outfilename);

return 0;
}


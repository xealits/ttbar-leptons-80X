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
if (argc < 5)
	{
	//std::cout << "Usage : " << argv[0] << " x|y|z proj_name distr distr2 rebin_factor x_axis_min_range x_axis_max_range filename" << std::endl;
	std::cout << "Usage : " << argv[0] << " axis_name content_name distr rebin_factor x_axis_min_range x_axis_max_range filename" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString axis_name(argv[1]);
TString cont_name(argv[2]);
TString distr(argv[3]);
Int_t rebin_factor(atoi(argv[4]));
double x_axis_min_range = atof(argv[5]);
double x_axis_max_range = atof(argv[6]);
TString filename(argv[7]);

cout << distr << endl;
//cout << distr2 << endl;
cout << "x axis minimum limit = " << x_axis_min_range << endl; // apparently, it doesn't work histogram's X axis -- you cannot zoom in within the histogram
cout << "x axis maximum limit = " << x_axis_max_range << endl; // apparently, it doesn't work histogram's X axis -- you cannot zoom in within the histogram
cout << filename   << endl;

TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);
//TLegend* leg = new TLegend(0.845, 0.5, 0.99, 0.99);

TFile * file = TFile::Open(filename);

TH1D *h3 = (TH1D*) file->Get(distr);
h3->SetLineColor(kBlack);
//h3->SetMinimum(0.8);  // Define Y ..
//h3->SetMaximum(1.35); // .. range
h3->Sumw2();
h3->SetStats(0);      // No statistics on lower plot
//h3->Divide(h2);
h3->SetMarkerStyle(21);

// normalize h3 for bin width
for (Int_t i=0; i<=h3->GetSize(); i++)
	{
	//yAxis->GetBinLowEdge(3)
	double content = h3->GetBinContent(i);
	double width   = h3->GetXaxis()->GetBinUpEdge(i) - h3->GetXaxis()->GetBinLowEdge(i);
	h3->SetBinContent(i, content/width);
	}

h3->SetTitle(distr);
h3->SetXTitle(axis_name);
h3->SetYTitle(cont_name);

h3->Draw("ep");       // Draw the ratio plot

//h3->GetYaxis()->SetRange(0.0001, 1); // ranges from analysis note CMS AN-2012/489
//h3->GetYaxis()->SetRangeUser(0.0001, 1); // ranges from analysis note CMS AN-2012/489
h3->GetXaxis()->SetRange(x_axis_min_range, x_axis_max_range);
h3->GetXaxis()->SetRangeUser(x_axis_min_range, x_axis_max_range);

cout << "bin-scaled distr:" << endl;
h3->Print();

TFile f( filename.ReplaceAll(".root","") + "_" + distr + ".root", "create" );
h3->Write();
f.Write();
f.Close();


h3->Draw("ep");       // Draw the ratio plot

cst->SetLogy();
//leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");

cst->Update();

//hs->GetXaxis()->SetTitle(distr);
//hs_data->SetXTitle(distr);

cst->Modified();

cst->SaveAs( filename.ReplaceAll(".root","") + "_" + distr + ".png" );

return 0;
}


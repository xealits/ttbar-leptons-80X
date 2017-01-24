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

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = " lumi distr projection rebin_factor dir dtags";
if (argc < 5)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

double lumi = atof(argv[1]);
TString distr_selection(argv[2]);
TString projection(argv[3]);
Int_t rebin_factor(atoi(argv[4]));
TString dir(argv[5]);
TString dtag1(argv[INPUT_DTAGS_START]);

if (projection != TString("x") && projection != TString("y") && projection != TString("z"))
	{
	printf("UNKNOWN PROJECTION GIVEN\n");
	printf("supported: x, y, z\n");
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 2;
	}

cout << lumi  << endl;
cout << distr_selection << endl;
cout << dir   << endl;
cout << dtag1 << endl;


std::vector < TString > dtags;
std::vector < TFile * > files;
std::vector < TH1D * > histos;
std::vector < TH1D * > weightflows;
// nick->summed histo
//std::map<TString, TH1D *> nicknamed_mc_histos;
//vector<int> dtags;
//dtags.reserve();

/*
 * sum a data histo
 * scale and sum different jet origin histos
 * then make a stack of them
 * and plot the stack and data
 */

// make stack of MC, scaling according to ratio = lumi * xsec / weightflow4 (bin5?)
// also nickname the MC....
// per-dtag for now..

TH1D    *hs_data = NULL;

// different jet origin histos:
vector<string> mc_jet_origins = {"_jets_distr_o", "_jets_distr_t", "_jets_distr_b", "_jets_distr_q", "_jets_distr_g"};
vector<TH1D*> mc_jet_origin_ths = {NULL, NULL, NULL, NULL, NULL};
//TH1D    *hs_mc_o = NULL; // "other" objects -- not recognized by partonFlavour
//TH1D    *hs_mc_t = NULL; // taus -- matched not recognized by partonFlavour and matched to a visible part of tau in GenParticles collection
//TH1D    *hs_mc_b = NULL; // b
//TH1D    *hs_mc_q = NULL; // not b, light quarks
//TH1D    *hs_mc_g = NULL; // gluons
THStack *hs      = new THStack("hs", "");
//THStack *hs = new THStack("hs","Stacked 1D histograms");


TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);

//TLegend *leg = new TLegend(0.845, 0.2, 0.99, 0.99);
//leg = new TLegend(0.845, 0.2, 0.99, 0.99);
TLegend* leg = new TLegend(0.845, 0.5, 0.99, 0.99);

/*
 * Get each dtag file in the reduced dir,
 * open the histogram, project it, rebin the projection,
 * if the dtag is Data -> sum the projection to data-histo
 * if it is MC -> scale it according xsecs/weightflow and stack it to mc-stack
 *
 */
for (int i = INPUT_DTAGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	TString the_file = dir + "/" + dtag + ".root";
	cout << the_file << endl;
	TFile* file = TFile::Open(the_file);
	//files.push_back(file); // not needed
	//dtags.push_back(5);

	/*
	 * data containd distr_selection + _jet_distr
	 * MC: + jet_distr_g/q/b/t/o
	 */
	/*
	if (!file->GetListOfKeys()->Contains(distr_selection))
		{
		cout << "no " << distr_selection << endl;
		continue;
		}
	*/

	if (dtag.Contains("Data"))
		{
		cout << "summing data-stack" << endl;

		if (!file->GetListOfKeys()->Contains(distr_selection + string("_jets_distr")))
			{
			cout << "no " << distr_selection +string("_jets_distr") << endl;
			continue;
			}

		// get the histogram's projection
		TH1D* histo = (TH1D*) ((TH3D*) file->Get(distr_selection + string("_jets_distr")))->Project3D(projection);
		//histos.push_back();
		//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
		histo->Print();

		// normalize the histo for bin width
		for (Int_t i=0; i<=histo->GetSize(); i++)
			{
			//yAxis->GetBinLowEdge(3)
			double content = histo->GetBinContent(i);
			double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
			histo->SetBinContent(i, content/width);
			}

		histo->SetMarkerStyle(9);
		//histo->SetFillColor(kRed + i);

		if (hs_data == NULL)
			{
			cout << "creating data histo" << endl;
			hs_data = (TH1D*) histo->Clone();
			}
		else
			{
			cout << "add histo to data histo" << endl;
			hs_data->Add(histo);
			}
		}

	else
		for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
			{
			string jet_origin = mc_jet_origins[origin_n];
			if (!file->GetListOfKeys()->Contains(distr_selection + jet_origin))
				{
				cout << "no " << distr_selection + jet_origin << endl;
				continue;
				}

			// get the histogram's projection
			TH1D* histo = (TH1D*) ((TH3D*) file->Get(distr_selection + jet_origin))->Project3D(projection);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			histo->Print();

			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				}

			// different weightflow MC normalization in ttbar and jet fake rates (TODO: need to normalize it to bin10 or something everywhere)
			TH1D * weightflow;
			// actually, only 1 number will be needed:
			double normal_initial_weight = 0;
			if (file->GetListOfKeys()->Contains("weightflow_elel"))
				{
				weightflow = (TH1D*) file->Get("weightflow_elel");
				normal_initial_weight = weightflow->GetBinContent(5);
				}
			else if (file->GetListOfKeys()->Contains("weightflow"))
				{
				weightflow = (TH1D*) file->Get("weightflow");
				normal_initial_weight = weightflow->GetBinContent(4);
				}
			else
				{
				cerr << "no weightflow distro" << endl;
				return 1;
				}
			// not needed:
			//weightflows.push_back(weightflow);
			//weightflows.back()->Print();
			// --- using only the initial weight

			cout << "got weightflow init" << endl;

			// MC ratio for this dtag:
			Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
			histo->Scale(ratio);
			cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
			histo->Print();
			//histo->SetFillColor(kRed );

			/*
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
			TString nick = nick_colour.first;
			Color_t col = nick_colour.second;
			histo->SetFillColor( col );
			*/

			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor(origin_n);
			histo->SetFillColor(origin_n);

			// add the histo to an origin histo
			// or clone, if the origin hiso is not initialized already
			if (mc_jet_origin_ths[origin_n])
				mc_jet_origin_ths[origin_n]->Add(histo);
			else
				mc_jet_origin_ths[origin_n] = (TH1D*) histo->Clone();
			}

	cout << "processed dtag" << endl;
	}

for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	hs->Add(mc_jet_origin_ths[origin_n], "HIST");
	leg->AddEntry(mc_jet_origin_ths[origin_n], TString(mc_jet_origins[origin_n]), "F");
	}

//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

cout << "setting title" << endl;

//hs->GetXaxis()->SetTitle(distr);
//cst->SetXaxisTile(distr);
//hs_data->GetXaxis()->SetTitle(distr);


TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("x");
//cst->Update();
//cst->Modified();

cout << "drawing" << endl;

hs_data->Draw("e p");
hs->Draw("same");
hs_data->Draw("e p same"); // to draw it _over_ MC

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


hs->GetXaxis()->SetTitle(distr_selection);
hs_data->SetXTitle(distr_selection);

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr_selection + "_OriginStacked_" + projection +".png" );

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


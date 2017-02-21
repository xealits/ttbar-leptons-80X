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

#define INPUT_DTAGS_START 7

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = "[--verbose] [--normalize] lumi distr projection rebin_factor name_tag dir dtags";
if (argc < 6)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

unsigned int input_starts = 0;
bool be_verbose = false;
bool normalize_MC = false;
string inp1(argv[1]), inp2(argv[2]);
if (inp1 == string("--verbose"))
	{
	input_starts += 1;
	be_verbose = true;
	if (inp2 == string("--normalize"))
		{
		normalize_MC = true;
		input_starts += 1;
		}
	}
else if (inp1 == string("--normalize"))
	{
	normalize_MC = true;
	input_starts += 1;
	}

if (be_verbose) cout << "being verbose" << endl;
if (normalize_MC) cout << "will normalize MC sets to 1" << endl;
cout << "options are taken from " << input_starts << endl;

double lumi = atof(argv[input_starts + 1]);
TString distr_selection(argv[input_starts + 2]);
TString projection(argv[input_starts + 3]);
Int_t rebin_factor(atoi(argv[input_starts + 4]));
TString name_tag(argv[input_starts + 5]);
TString dir(argv[input_starts + 6]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

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
//THStack *hs      = new THStack("hs", "");
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
for (int i = input_starts + INPUT_DTAGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	if (be_verbose) cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
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

	// hack to merge HLTmu & HLTjetmu
	//string distr_selections[2] = {"HLTmu_wjets", "HLTjetmu_wjets"};
	//for (int i = 0; i<2; i++){
	//TString distr_selection(distr_selections[i]);
	if (dtag.Contains("Data"))
		{
		if (be_verbose) cout << "summing data-stack" << endl;
		TString distro_name = distr_selection + string("_jets_distr");

		if (!file->GetListOfKeys()->Contains(distro_name))
			{
			if (be_verbose) cout << "no " << distro_name << endl;
			continue;
			}

		// get the histogram's projection
		TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);
		//histos.push_back();
		//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
		if (be_verbose) histo->Print();

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
			if (be_verbose) cout << "creating data histo" << endl;
			hs_data = (TH1D*) histo->Clone();
			}
		else
			{
			if (be_verbose) cout << "add histo to data histo" << endl;
			hs_data->Add(histo);
			}
		}

	else
		for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
			{
			string jet_origin = mc_jet_origins[origin_n];
			TString distro_name = distr_selection + jet_origin;

			if (!file->GetListOfKeys()->Contains(distro_name))
				{
				if (be_verbose) cout << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			if (be_verbose) histo->Print();

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

			if (be_verbose) cout << "got weightflow init" << endl;

			// MC ratio for this dtag:
			Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
			histo->Scale(ratio);
			if (be_verbose) cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
			if (be_verbose) histo->Print();
			//histo->SetFillColor(kRed );

			/*
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
			TString nick = nick_colour.first;
			Color_t col = nick_colour.second;
			histo->SetFillColor( col );
			*/

			//histo->SetMarkerStyle(20);
			//histo->SetLineStyle(1);
			histo->SetLineWidth(3);
			histo->SetLineColor(origin_n!=0? origin_n : 6);
			//histo->SetMarkerColor(origin_n);
			//histo->SetFillColor(origin_n);

			// add the histo to an origin histo
			// or clone, if the origin hiso is not initialized already
			if (mc_jet_origin_ths[origin_n])
				mc_jet_origin_ths[origin_n]->Add(histo);
			else
				mc_jet_origin_ths[origin_n] = (TH1D*) histo->Clone();
			}
		//}

	if (be_verbose) cout << "processed dtag" << endl;
	}

double integral_data = hs_data->Integral();
double integral_MC   = 0.;
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	integral_MC += mc_jet_origin_ths[origin_n]->Integral();
	}
cout << "data  integral = " << integral_data << endl;
cout << "MC    integral = " << integral_MC   << endl;
double ratio = integral_data / integral_MC;
cout << "ratio = " << ratio << endl;

if (normalize_MC)
	{
	cout << "normalizing MCs" << endl;
	for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
		{
		double integral = mc_jet_origin_ths[origin_n]->Integral();
		double ratio = 1. / integral;
		mc_jet_origin_ths[origin_n]->Scale(ratio);
		cout << mc_jet_origins[origin_n] << " integral before = " << integral << " after = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		}
	}

// add MC to legend
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
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

//hs_data->Draw("e"); // to draw it _over_ MC

if (projection == TString("z"))
	{
	mc_jet_origin_ths[mc_jet_origins.size()-1]->GetXaxis()->SetRange(0, 0.3);
	mc_jet_origin_ths[mc_jet_origins.size()-1]->GetXaxis()->SetRangeUser(0, 0.3);
	}

mc_jet_origin_ths[mc_jet_origins.size()-1]->SetLineStyle(1);
mc_jet_origin_ths[mc_jet_origins.size()-1]->Draw("hist");
for (int origin_n=mc_jet_origins.size()-2; origin_n>=0; origin_n--)
	{
	mc_jet_origin_ths[origin_n]->SetLineStyle(1);
	mc_jet_origin_ths[origin_n]->Draw("same hist");
	}

/*
mc_jet_origin_ths[0]->GetXaxis()->SetRange(0, 0.4);
mc_jet_origin_ths[0]->GetXaxis()->SetRangeUser(0, 0.4);


mc_jet_origin_ths[0]->GetYaxis()->SetRange(0, 10000);
mc_jet_origin_ths[0]->GetYaxis()->SetRangeUser(0, 10000);
*/

bool set_logy = false;

if (projection == string("x"))
	set_logy = true;

if (set_logy)
	cst->SetLogy();

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


//hs_data->GetXaxis()->SetTitle(distr_selection);
//hs_data->SetXTitle(distr_selection);
//mc_jet_origin_ths[0]->GetXaxis()->SetTitle("Jet Radius for different origins");
mc_jet_origin_ths[0]->SetXTitle(projection);

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr_selection + "_DistrsOriginOverlayed_" + projection + "_" + name_tag + (normalize_MC ? "_normalized" : "") + (set_logy? "_logy" : "") + ".png" );

return 0;
}


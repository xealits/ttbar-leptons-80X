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

#define INPUT_DTAGS_START 10

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = " [--verbose] [--normalize] unstack largebins csvOUT lumi distr projection rebin_factor name_tag dir dtags";
if (argc < INPUT_DTAGS_START)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

// jet_distrs_origin_fractions_stacked
unsigned int input_starts = 0;
bool be_verbose = false;
bool normalize_MC_stack = false;
string inp1(argv[1]), inp2(argv[2]);
if (inp1 == string("--verbose"))
	{
	input_starts += 1;
	be_verbose = true;
	if (inp2 == string("--normalize"))
		{
		normalize_MC_stack = true;
		input_starts += 1;
		}
	}
else if (inp1 == string("--normalize"))
	{
	normalize_MC_stack = true;
	input_starts += 1;
	}

if (be_verbose) cout << "being verbose" << endl;
if (normalize_MC_stack) cout << "will normalize MC stack to Data integral" << endl;
cout << "options are taken from " << input_starts << endl;

int i = 1;

TString to_unstack(argv[input_starts + i++]);
bool unstack=false;
if (to_unstack == TString("T") || to_unstack == TString("Y"))
	{
	cout << "unstacking!" << endl;
	unstack = true;
	}

string largebins(argv[input_starts + i++]);

TString csvOUT(argv[input_starts + i++]);
bool csv_out=false;
if (csvOUT == TString("T") || csvOUT == TString("Y"))
	csv_out = true;

double lumi = atof(argv[input_starts + i++]);
TString distr_selection(argv[input_starts + i++]);
TString projection(argv[input_starts + i++]);
Int_t rebin_factor(atoi(argv[input_starts + i++]));
TString name_tag(argv[input_starts + i++]);
TString dir(argv[input_starts + i++]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

if (projection != TString("x") && projection != TString("y") && projection != TString("z"))
	{
	printf("UNKNOWN PROJECTION GIVEN\n");
	printf("supported: x, y, z\n");
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 2;
	}

cout << "largebins " << largebins  << endl;
cout << "CSV OUT " << csvOUT  << endl;
cout << "lumi " << lumi  << endl;
cout << "distr " << distr_selection << endl;
cout << "dir " << dir   << endl;
cout << "first dtag " << dtag1 << endl;


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
//vector<string> mc_jet_origins = {"_jets_distr_o", "_jets_distr_t", "_jets_distr_b", "_jets_distr_q", "_jets_distr_g"};
//vector<string> mc_jet_origins = {"_jets_distr_q", "_jets_distr_b", "_jets_distr_g", "_jets_distr_o", "_jets_distr_t"};
vector<string> mc_jet_origins = {"q", "b", "t", "g", "o"};
vector<TH1D*> mc_jet_origin_ths = {NULL, NULL, NULL, NULL, NULL};

// just for control:
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	string jet_origin = "_jets_distr" + largebins + "_" + mc_jet_origins[origin_n];
	cout << jet_origin << endl;
	}

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
			string jet_origin = "_jets_distr" + largebins + "_" + mc_jet_origins[origin_n];
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
			if (file->GetListOfKeys()->Contains("weightflow_elel_NOMINAL"))
				{
				weightflow = (TH1D*) file->Get("weightflow_elel_NOMINAL");
				normal_initial_weight = weightflow->GetBinContent(11);
				}
			else if (file->GetListOfKeys()->Contains("weightflow_elel"))
				{
				weightflow = (TH1D*) file->Get("weightflow_elel");
				normal_initial_weight = weightflow->GetBinContent(11);
				}
			else if (file->GetListOfKeys()->Contains("weightflow"))
				{
				weightflow = (TH1D*) file->Get("weightflow");
				//normal_initial_weight = weightflow->GetBinContent(11); // not yet..
				normal_initial_weight = weightflow->GetBinContent(11);
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

			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor((origin_n != 0 ? origin_n : 14));
			histo->SetFillColor((origin_n != 0 ? origin_n : 14));

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

// sum of all origins
// for getting the fractions
TH1D* mc_sum = (TH1D*) mc_jet_origin_ths[4]->Clone();
mc_sum->Reset();

double integral_data = hs_data->Integral();
double integral_MC   = 0.;
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	if (!mc_jet_origin_ths[origin_n])
		{
		cout << "no " << mc_jet_origins[origin_n] << endl;
		// let's hope the 0-s is there
		mc_jet_origin_ths[origin_n] = (TH1D*) mc_jet_origin_ths[0]->Clone();
		mc_jet_origin_ths[origin_n]->Reset();
		}
	integral_MC += mc_jet_origin_ths[origin_n]->Integral();
	mc_sum->Add(mc_jet_origin_ths[origin_n]);
	}
cout << "data  integral = " << integral_data << endl;
cout << "MC    integral = " << integral_MC   << endl;
double ratio = integral_data / integral_MC;
cout << "ratio = " << ratio << endl;

/* TODO: somehow it doesn't work
// normalize data for bin width
for (Int_t i=0; i<=hs_data->GetSize(); i++)
	{
	//yAxis->GetBinLowEdge(3)
	double content = hs_data->GetBinContent(i);
	double width   = hs_data->GetXaxis()->GetBinUpEdge(i) - hs_data->GetXaxis()->GetBinLowEdge(i);
	hs_data->SetBinContent(i, content/width);
	}

// normalize each MC origin for bin size
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	// normalize data for bin width
	for (Int_t i=0; i<=mc_jet_origin_ths[origin_n]->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = mc_jet_origin_ths[origin_n]->GetBinContent(i);
		double width   = mc_jet_origin_ths[origin_n]->GetXaxis()->GetBinUpEdge(i) - mc_jet_origin_ths[origin_n]->GetXaxis()->GetBinLowEdge(i);
		mc_jet_origin_ths[origin_n]->SetBinContent(i, content/width);
		}
	}
*/

/*
if (normalize_MC_stack) // TODO: does it interfere with the bin-wise normalization?
	{
	cout << "normalizing MC" << endl;
	for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
		{
		double integral = mc_jet_origin_ths[origin_n]->Integral();
		mc_jet_origin_ths[origin_n]->Scale(ratio);
		cout << mc_jet_origins[origin_n] << " integral before = " << integral << " after = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		}
	}
*/

// normalized bins to the sum of all jets (to get the fraction of the origin in each bin)
for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	{
	TH1D* histo = mc_jet_origin_ths[origin_n];
	for (Int_t i=0; i<=histo->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = histo->GetBinContent(i);
		double mc_sum_content = mc_sum->GetBinContent(i);
		if (mc_sum_content < 0.0001) // zero case
			histo->SetBinContent(i, origin_n/mc_jet_origins.size());
		else
			histo->SetBinContent(i, content/mc_sum_content);
		}
	}

/* not needed then
if (normalize_MC_stack) // TODO: does it interfere with the bin-wise normalization?
	{
	cout << "normalizing MC" << endl;
	for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
		{
		double integral = mc_jet_origin_ths[origin_n]->Integral();
		mc_jet_origin_ths[origin_n]->Scale(ratio);
		cout << mc_jet_origins[origin_n] << " integral before = " << integral << " after = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		}
	}
*/

// first add the light-quarks and gluon-jets
for (int origin_n=3; origin_n<mc_jet_origins.size(); origin_n++)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	hs->Add(mc_jet_origin_ths[origin_n], "HIST");
	leg->AddEntry(mc_jet_origin_ths[origin_n], TString(mc_jet_origins[origin_n]), "F");
	}
// then the rest
for (int origin_n=0; origin_n<3; origin_n++)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
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
h->SetXTitle(projection);
//cst->Update();
//cst->Modified();

cout << "drawing" << endl;

//hs_data->Draw("e p");
//hs->Draw("same");
if (!unstack)
	hs->Draw();
else
	{// loop and overay histos
	cout << 0 << " " << mc_jet_origins[0];
	if (mc_jet_origin_ths[0])
		{
		mc_jet_origin_ths[0]->GetYaxis()->SetRange(0, 1);
		mc_jet_origin_ths[0]->GetYaxis()->SetRangeUser(0, 1);
		if (projection == TString("z"))
			{
			mc_jet_origin_ths[0]->GetXaxis()->SetRange(0, 0.4);
			mc_jet_origin_ths[0]->GetXaxis()->SetRangeUser(0, 0.4);
			}
		mc_jet_origin_ths[0]->Draw("p");
		cout << " done" << endl;
		}
	else
		{
		cout << " NONE" << endl;
		}
	for (int origin_n=1; origin_n<mc_jet_origins.size(); origin_n++)
		{
		//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		//mc_jet_origin_ths[origin_n]->Print();
		cout << origin_n << " " << mc_jet_origins[origin_n];
		if (mc_jet_origin_ths[origin_n])
			{
			mc_jet_origin_ths[origin_n]->GetYaxis()->SetRange(0, 1);
			mc_jet_origin_ths[origin_n]->GetYaxis()->SetRangeUser(0, 1);
			if (projection == TString("z"))
				{
				mc_jet_origin_ths[origin_n]->GetXaxis()->SetRange(0, 0.4);
				mc_jet_origin_ths[origin_n]->GetXaxis()->SetRangeUser(0, 0.4);
				}
			mc_jet_origin_ths[origin_n]->Draw("same p");
			cout << " done" << endl;
			}
		else
			{
			cout << " NONE" << endl;
			}
		}
	}
//hs_data->Draw("e p same"); // to draw it _over_ MC

bool set_logy = false;

//if (projection == string("x") || projection == string("z"))
//	set_logy = true;

if (set_logy)
	cst->SetLogy();

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


if (!unstack)
	{
	hs->GetXaxis()->SetTitle(distr_selection);
	hs_data->SetXTitle(distr_selection);

	cst->Modified();
	}

TString output_file = dir + "/jobsums/" + distr_selection + (unstack? "_OriginFractionsUnstacked_" : "_OriginFractionsStacked_") + projection + "_" + name_tag + (normalize_MC_stack ? "_normalized" : "") + (set_logy? "_logy" : "");

cst->SaveAs( output_file + ".png" );

TFile* out_file = TFile::Open(output_file + ".root", "RECREATE");
//TFile* out_file = TFile::Open(output_file + ".root");
if (unstack)
	{// loop and overay histos
	cout << 0 << " " << mc_jet_origins[0];
	if (mc_jet_origin_ths[0])
		{
		mc_jet_origin_ths[0]->Write();
		cout << " done" << endl;
		}
	else
		{
		cout << " NONE" << endl;
		}
	for (int origin_n=1; origin_n<mc_jet_origins.size(); origin_n++)
		{
		//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		//mc_jet_origin_ths[origin_n]->Print();
		cout << origin_n << " " << mc_jet_origins[origin_n];
		if (mc_jet_origin_ths[origin_n])
			{
			mc_jet_origin_ths[origin_n]->Write();
			cout << " done" << endl;
			}
		else
			{
			cout << " NONE" << endl;
			}
		}
}
out_file->Write();
out_file->Close();

if (csv_out)
	{
	// normalized bins to the sum of all jets (to get the fraction of the origin in each bin)
	// header
	cout << "bin_x";
	for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
		{
		cout << "," + mc_jet_origins[origin_n];
		}
	cout << endl;

	TH1D* tamplate_histo = mc_jet_origin_ths[0];
	for (Int_t i=0; i<=tamplate_histo->GetSize(); i++)
		{
		cout << tamplate_histo->GetBinCenter(i);
		for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
			{
			cout << "," << mc_jet_origin_ths[origin_n]->GetBinContent(i);
			}
		cout << endl;
		}
	}

return 0;
}



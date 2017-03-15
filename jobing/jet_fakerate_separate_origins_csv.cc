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

#define INPUT_DTAGS_START 9

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = "[--verbose] [--normalize] lumi distr projection rebin_factor x_axis_min_range x_axis_max_range name_tag dir dtags";
if (argc < 8)
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
if (normalize_MC) cout << "will normalize MC stack to Data integral" << endl;
cout << "options are taken from " << input_starts << endl;

double lumi = atof(argv[input_starts + 1]);
TString distr_selection(argv[input_starts + 2]);
TString projection(argv[input_starts + 3]);
Int_t rebin_factor(atoi(argv[input_starts + 4]));

double x_axis_min_range = atof(argv[input_starts + 5]);
double x_axis_max_range = atof(argv[input_starts + 6]);
TString name_tag(argv[input_starts + 7]);
TString dir(argv[input_starts + 8]);
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
std::vector < TH3D * > histos;
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

TH1D* hs_data[2] = {NULL, NULL};

// different jet origin histos:
unsigned int n_jet_origins = 5;
vector<string> mc_jet_origins = {"_tau_jets_distr_o", "_tau_jets_distr_t", "_tau_jets_distr_b", "_tau_jets_distr_q", "_tau_jets_distr_g",
	"_jets_distr_o", "_jets_distr_t", "_jets_distr_b", "_jets_distr_q", "_jets_distr_g"};
//HLTjetmu_qcd_tau_jets_distr_q
vector<TH1D*> mc_jet_origin_ths = {NULL, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL, NULL};
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
	cout << dtag << endl;
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
		vector<TString> distro_names;
		distro_names.push_back(distr_selection + string("_tau_jets_distr"));
		distro_names.push_back(distr_selection + string("_jets_distr"));

		for (int i=0; i<distro_names.size(); i++)
			{
			TString distro_name = distro_names[i];

			if (!file->GetListOfKeys()->Contains(distro_name))
				{
				if (be_verbose) cout << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			//TH3D* histo = (TH3D*) file->Get(distro_name); //for full histo
			TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			if (be_verbose) histo->Print();

			/*
			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				}
			*/

			histo->SetMarkerStyle(9);
			//histo->SetFillColor(kRed + i);

			if (hs_data[i] == NULL)
				{
				if (be_verbose) cout << "creating data histo" << endl;
				hs_data[i] = (TH1D*) histo->Clone();
				if (be_verbose) cout <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data[i]->Integral() << endl;
				}
			else
				{
				if (be_verbose) cout << "add histo to data histo" << endl;
				hs_data[i]->Add(histo);
				if (be_verbose) cout <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data[i]->Integral() << endl;
				}
			}
		}

	else
		for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
			{
			string jet_origin = mc_jet_origins[origin_n];
			TString distro_name = distr_selection + jet_origin;

			if (!file->GetListOfKeys()->Contains(distro_name))
				{
				cout << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			//TH3D* histo = (TH3D*) file->Get(distro_name);
			TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			if (be_verbose) histo->Print();

			/*
			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				}
			*/

			// different weightflow MC normalization in ttbar and jet fake rates (TODO: need to normalize it to bin10 or something everywhere)
			TH1D * weightflow;
			// actually, only 1 number will be needed:
			double normal_initial_weight = 0;
			if (file->GetListOfKeys()->Contains("weightflow_elel"))
				{
				weightflow = (TH1D*) file->Get("weightflow_elel");
				normal_initial_weight = weightflow->GetBinContent(11);
				}
			else if (file->GetListOfKeys()->Contains("weightflow"))
				{
				weightflow = (TH1D*) file->Get("weightflow");
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
			if (be_verbose) cout << "x-sec/weight/ratio\t" << xsecs[dtag] << "\t" << normal_initial_weight << "\t" << ratio << endl;
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
			//histo->SetLineStyle(1);
			//histo->SetLineWidth(3);
			//histo->SetLineColor(1 + (origin_n<4? origin_n : origin_n+1 ));
			//histo->SetMarkerColor(1 + origin_n);
			histo->SetMarkerColor(origin_n!=0? origin_n : 6);
			histo->SetFillColor(origin_n);

			// add the histo to an origin histo
			// or clone, if the origin hiso is not initialized already
			if (mc_jet_origin_ths[origin_n])
				{
				// in principle, on addition one can get floating point saturation
				//cout << "" << hs_data[i]->Integral() << "\t" << histo->Integral() << endl;

				mc_jet_origin_ths[origin_n]->Add(histo);
				if (be_verbose) cout << "Float control: " << mc_jet_origins[origin_n] << "\t" << histo->Integral() << '\t' << mc_jet_origin_ths[origin_n]->Integral() << endl;
				}
			else
				{
				mc_jet_origin_ths[origin_n] = (TH1D*) histo->Clone();
				if (be_verbose) cout << "Float control: " << mc_jet_origins[origin_n] << "\t" << histo->Integral() << '\t' << mc_jet_origin_ths[origin_n]->Integral() << endl;
				}
			}
		//}

	if (be_verbose) cout << "processed dtag" << endl;
	}

/*
TH1D* hs_mc_all_jets = (TH1D*) mc_jet_origin_ths[0 + n_jet_origins]->Clone();
hs_mc_all_jets->Reset();
*/

double integral_data = hs_data[1]->Integral();
double integral_MC_taus  = 0.;
double integral_MC_jets  = 0.;
for (int origin_n=0; origin_n<n_jet_origins; origin_n++)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	if (!mc_jet_origin_ths[origin_n]) // let's hope the 1 (tau_jets.._o) exists
		{
		mc_jet_origin_ths[origin_n] = (TH1D*) mc_jet_origin_ths[0]->Clone();
		mc_jet_origin_ths[origin_n]->Reset();
		}
	if (!mc_jet_origin_ths[origin_n + n_jet_origins])
		{
		mc_jet_origin_ths[origin_n + n_jet_origins] = (TH1D*) mc_jet_origin_ths[0]->Clone();
		mc_jet_origin_ths[origin_n + n_jet_origins]->Reset();
		}

	integral_MC_taus += mc_jet_origin_ths[origin_n]->Integral();
	integral_MC_jets += mc_jet_origin_ths[n_jet_origins + origin_n]->Integral();

	// No mc_all_jets --- the jet origins separately
	// Float control
	//cout << "Float control: " << "mc_all_jets" << "\t" << mc_jet_origin_ths[n_jet_origins + origin_n]->Integral() << "\t" << hs_mc_all_jets->Integral() << endl;
	//hs_mc_all_jets->Add(mc_jet_origin_ths[n_jet_origins + origin_n]);
	if (be_verbose) cout << mc_jet_origins[origin_n] << "\t" << mc_jet_origins[n_jet_origins + origin_n] << endl;
	}
cout << "data    integral = " << "(" << hs_data[0]->Integral() << ")" << integral_data << endl;
cout << "MC taus integral = " << integral_MC_taus << endl;
cout << "MC jets integral = " << integral_MC_jets << endl;
//double ratio = integral_data / integral_MC;
//cout << "ratio = " << ratio << endl;

// find the fake rate in data:
hs_data[0]->Sumw2();
hs_data[0]->SetStats(0);      // No statistics on lower plot
hs_data[0]->Divide(hs_data[1]);

// divide tau_jets / jets
for (int origin_n=0; origin_n<n_jet_origins; origin_n++)
	{
	TH1D* distr = mc_jet_origin_ths[origin_n];
	// calculate fake rate
	// the origin tau jets / all jets (of all origins)
	//distr->Divide(hs_mc_all_jets);
	distr->Divide(mc_jet_origin_ths[n_jet_origins + origin_n]);
	//mc_jet_origin_ths[origin_n]->Divide(mc_jet_origin_ths[n_jet_origins + origin_n]);

	// and add the histogram to the stack and record about it to the legend
	//hs->Add(distr, "HIST");
	//leg->AddEntry(distr, TString(mc_jet_origins[origin_n]), "F");
	// rebin here, if needed
	}

cout << "setting title" << endl;


TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("x");
//cst->Update();
//cst->Modified();

cout << "drawing" << endl;

bool set_logy = true; // in case it becomes main argument

//if (projection == string("x") || projection == string("z"))
	//set_logy = true;

if (set_logy)
	cst->SetLogy();

//leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


//hs_data->GetXaxis()->SetTitle(distr_selection);
//hs_data->SetXTitle(distr_selection);
//mc_jet_origin_ths[0]->GetXaxis()->SetTitle(name_tag);
//mc_jet_origin_ths[0]->SetXTitle(projection);

//cst->Modified();

//vector<string> mc_jet_origins = {"_tau_jets_distr_o", "_tau_jets_distr_t", "_tau_jets_distr_b", "_tau_jets_distr_q", "_tau_jets_distr_g",
cout << "other,tau,b,q,g" << endl;
for (int bin=0; bin<mc_jet_origin_ths[0]->GetSize(); bin++)
	{
	double other = mc_jet_origin_ths[0]->GetBinContent(bin);
	double tau   = mc_jet_origin_ths[1]->GetBinContent(bin);
	double b     = mc_jet_origin_ths[2]->GetBinContent(bin);
	double q     = mc_jet_origin_ths[3]->GetBinContent(bin);
	double g     = mc_jet_origin_ths[4]->GetBinContent(bin);
	cout << other << ',' << tau << ',' << b << ',' << q << ',' << g << endl;
	}
cout << endl;

return 0;
}


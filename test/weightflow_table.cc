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

#define INPUT_DTAGS_START 4

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
//char usage_string[128] = "[--verbose] [--normalize] lumi distr projection rebin_factor x_axis_min_range x_axis_max_range name_tag dir dtags";
const char usage_string[128] = "[--verbose] [--weight] lumi distr  dir dtags";
if (argc < 7)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	std::cout << "        distr = {weightflow_el_NOMINAL etc}" << endl;
	return 1;
	}

gROOT->Reset();

unsigned int input_starts = 0;
bool be_verbose = false;
bool weight_MC = false;
string inp1(argv[1]), inp2(argv[2]);
if (inp1 == string("--verbose"))
	{
	input_starts += 1;
	be_verbose = true;
	if (inp2 == string("--weight"))
		{
		weight_MC = true;
		input_starts += 1;
		}
	}
else if (inp1 == string("--weight"))
	{
	weight_MC = true;
	input_starts += 1;
	}

if (be_verbose) cout << "being verbose" << endl;
if (weight_MC)  cout << "will record weighted MC" << endl;
if (be_verbose) cout << "options are taken from " << input_starts << endl;

double lumi = atof(argv[input_starts + 1]);
TString distro_name(argv[input_starts + 2]);

TString dir(argv[input_starts + 3]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

if (be_verbose)
	{
	cout << lumi  << endl;
	cout << distro_name << endl;
	cout << dir   << endl;
	cout << dtag1 << endl;
	}


/*
 * Go through each dtag, print our weightflow selection
 * print sum of data at the end
 */

TH1D* hs_data_sum = NULL;
TH1D* hs_mc_sum   = NULL;

// bins with root-shift (the overflow bin is 0, the 0 bin is 1, 1 is 2...)
int         our_selection_bins[] = {1 + 1, 2 + 1, 3 + 1, 4 + 1, 10 + 1, 11 + 1, 12 + 1, 13 + 1, 20 + 1, 30 + 1,    31 + 1, 32 + 1, 34 + 1, 38 + 1, 46 + 1, 62 + 1};
std::vector<string> our_selection_names = {"ini", "top", "gen", "PUw", "sec1", "METf", "Lumi", "Trig", "sec2", "sec3",    "s", "sj", "sjm", "sjmb", "sjmbt", "sjmbto"};

if (be_verbose)
	{
	for (int i=0; i<our_selection_names.size(); i++)
		cout << our_selection_bins[i] << '\t' << our_selection_names[i] << endl;
	}

cout << "dtag,ratio" ;
for (int i=0; i<our_selection_names.size(); i++)
	cout << "," << our_selection_names[i];
cout << endl;


/*
 * Get each dtag file in the reduced dir,
 * open the histogram, get lumi ratio (for MC, for data it = 1)
 * cout dtag,ratio,sclaed weightflow bins
 * sum data and MC
 *
 */
for (int i = input_starts + INPUT_DTAGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose)
		{
		cout << the_file << endl;
		cout << dtag << endl;
		}
	TFile* file = TFile::Open(the_file);

	if (dtag.Contains("Data"))
		{
		if (be_verbose) cout << "data dtag" << endl;

		if (!file->GetListOfKeys()->Contains(distro_name))
			{
			if (be_verbose) cout << "no " << distro_name << endl;
			continue;
			}

		TH1D* histo = (TH1D*) file->Get(distro_name);

		// print the dtag line
		cout << dtag << "," << 1 ;
		for (int i=0; i<our_selection_names.size(); i++)
			cout << "," << histo->GetBinContent(our_selection_bins[i]);
		cout << endl;

		// sum the data
		if (hs_data_sum == NULL)
			{
			if (be_verbose) cout << "creating data histo" << endl;
			hs_data_sum = (TH1D*) histo->Clone();
			if (be_verbose) cout <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data_sum->Integral() << endl;
			}
		else
			{
			if (be_verbose) cout << "add histo to data histo" << endl;
			hs_data_sum->Add(histo);
			if (be_verbose) cout <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data_sum->Integral() << endl;
			}
		}

	else
		{
		if (!file->GetListOfKeys()->Contains(distro_name))
			{
			cout << "no " << distro_name << endl;
			continue;
			}

		TH1D* histo = (TH1D*) file->Get(distro_name);

		if (be_verbose) histo->Print();

		// different weightflow MC normalization in ttbar and jet fake rates
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
			//normal_initial_weight = weightflow->GetBinContent(11); // not yet...
			normal_initial_weight = weightflow->GetBinContent(11);
			}
		else
			{
			cerr << "no weightflow distro" << endl;
			return 1;
			}

		if (be_verbose) cout << "got weightflow init" << endl;

		// MC ratio for this dtag:
		Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
		if (be_verbose) cout << "x-sec/weight/ratio\t" << xsecs[dtag] << "\t" << normal_initial_weight << "\t" << ratio << endl;
		//if (weight_MC)

		if (lumi > 0)
			{
			histo->Scale(ratio);
			if (be_verbose) cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
			if (be_verbose) histo->Print();
			}

		// print the dtag line
		cout << dtag << "," << ratio ;
		for (int i=0; i<our_selection_names.size(); i++)
			cout << "," << histo->GetBinContent(our_selection_bins[i]);
		cout << endl;

		// sum the mc
		if (hs_mc_sum)
			{
			hs_mc_sum->Add(histo);
			if (be_verbose) cout << "Float control: " << '\t' << histo->Integral() << '\t' << hs_mc_sum->Integral() << endl;
			}
		else
			{
			hs_mc_sum = (TH1D*) histo->Clone();
			if (be_verbose) cout << "Float control: " << '\t' << histo->Integral() << '\t' << hs_mc_sum->Integral() << endl;
			}
		}

	if (be_verbose) cout << "processed dtag" << endl;
	}

if (hs_data_sum)
	{
	// print the data sum line
	cout << "sum_data" << "," << 1 ;
	for (int i=0; i<our_selection_names.size(); i++)
		cout << "," << hs_data_sum->GetBinContent(our_selection_bins[i]);
	cout << endl;
	}

if (hs_mc_sum)
	{
	// print the mc sum line
	cout << "sum_mc" << "," << 1 ;
	for (int i=0; i<our_selection_names.size(); i++)
		cout << "," << hs_mc_sum->GetBinContent(our_selection_bins[i]);
	cout << endl;
	}

if (hs_data_sum && hs_mc_sum)
	{
	hs_mc_sum->Divide(hs_data_sum);
	// print the data sum line
	cout << "ratio" << "," << 1 ;
	for (int i=0; i<our_selection_names.size(); i++)
		cout << "," << hs_mc_sum->GetBinContent(our_selection_bins[i]);
	cout << endl;
	}

return 0;
}


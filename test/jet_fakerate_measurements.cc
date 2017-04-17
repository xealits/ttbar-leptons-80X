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
char usage_string[128] = "[--verbose] [--add-header] distr dir dtags";

if (argc < 3)
	{
	std::cerr << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

unsigned int input_starts = 0;
bool be_verbose = false;
bool add_header = false;
string inp1(argv[1]), inp2(argv[2]);
if (inp1 == string("--verbose"))
	{
	input_starts += 1;
	be_verbose = true;
	if (inp2 == string("--add-header"))
		{
		add_header = true;
		input_starts += 1;
		}
	}
else if (inp1 == string("--add-header"))
	{
	add_header = true;
	input_starts += 1;
	}

if (argc < input_starts + INPUT_DTAGS_START)
	{
	std::cerr << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

if (be_verbose)
	{
	cerr << "being verbose" << endl;
	if (add_header) cerr << "will add header to the std output" << endl;
	cerr << "options are taken from " << input_starts << endl;
	}

TString distr_selection(argv[input_starts + 1]);
TString dir(argv[input_starts + 2]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

cerr << distr_selection << endl;
cerr << dir   << endl;
cerr << dtag1 << endl;


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
 * then stdout distr, Data_dakes, Data_all, MC_or1_fakes, MC_or1_all, ...
 */

// make stack of MC, scaling according to ratio = lumi * xsec / weightflow4 (bin5?)
// also nickname the MC....
// per-dtag for now..

TH3D* hs_data[2] = {NULL, NULL};
vector<TString> distro_names = {distr_selection + string("_tau_jets_distr"), distr_selection + string("_jets_distr")};

// different jet origin histos:
unsigned int n_jet_origins = 5;
vector<string> mc_jet_origins = {"_tau_jets_distr_o", "_tau_jets_distr_t", "_tau_jets_distr_b", "_tau_jets_distr_q", "_tau_jets_distr_g",
	"_jets_distr_o", "_jets_distr_t", "_jets_distr_b", "_jets_distr_q", "_jets_distr_g"};
//HLTjetmu_qcd_tau_jets_distr_q
vector<TH3D*> mc_jet_origin_ths = {NULL, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL, NULL};

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
	cerr << dtag << endl;
	dtags.push_back(dtag);
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cerr << the_file << endl;
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
		if (be_verbose) cerr << "summing data-stack" << endl;

		for (int i=0; i<distro_names.size(); i++)
			{
			TString distro_name = distro_names[i];

			if (!file->GetListOfKeys()->Contains(distro_name))
				{
				if (be_verbose) cerr << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			TH3D* histo = (TH3D*) file->Get(distro_name);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			if (be_verbose) histo->Print();

			/* don't
			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				}
			*/

			//histo->SetMarkerStyle(9);
			//histo->SetFillColor(kRed + i);

			if (hs_data[i] == NULL)
				{
				if (be_verbose) cerr << "creating data histo" << endl;
				hs_data[i] = (TH3D*) histo->Clone();
				if (be_verbose) cerr <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data[i]->Integral() << endl;
				}
			else
				{
				if (be_verbose) cerr << "add histo to data histo" << endl;
				hs_data[i]->Add(histo);
				if (be_verbose) cerr <<  "Float control: data" << i << '\t' << histo->Integral() << '\t' << hs_data[i]->Integral() << endl;
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
				cerr << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			TH3D* histo = (TH3D*) file->Get(distro_name);
			//histos.push_back();
			//histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
			if (be_verbose) histo->Print();

			/* don't
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
			if (file->GetListOfKeys()->Contains("weightflow_elel_NOMINAL"))
				{
				weightflow = (TH1D*) file->Get("weightflow_elel_NOMINAL");
				normal_initial_weight = weightflow->GetBinContent(11); // normalized weightflow
				}
			else if (file->GetListOfKeys()->Contains("weightflow"))
				{
				weightflow = (TH1D*) file->Get("weightflow");
				//normal_initial_weight = weightflow->GetBinContent(11); // not yet...
				normal_initial_weight = weightflow->GetBinContent(11); // TODO: make normalized weightflow in jet study
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

			if (be_verbose) cerr << "got weightflow init" << endl;

			// MC ratio for this dtag:
			//Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
			// lumi = 1
			Double_t ratio = xsecs[dtag] / normal_initial_weight;
			if (be_verbose) cerr << "x-sec/weight/ratio\t" << xsecs[dtag] << "\t" << normal_initial_weight << "\t" << ratio << endl;
			histo->Scale(ratio);
			if (be_verbose) cerr << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
			if (be_verbose) histo->Print();
			//histo->SetFillColor(kRed );

			/*
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
			TString nick = nick_colour.first;
			Color_t col = nick_colour.second;
			histo->SetFillColor( col );
			*/

			// add the histo to an origin histo
			// or clone, if the origin hiso is not initialized already
			if (mc_jet_origin_ths[origin_n])
				{
				// in principle, on addition one can get floating point saturation
				//cerr << "" << hs_data[i]->Integral() << "\t" << histo->Integral() << endl;

				mc_jet_origin_ths[origin_n]->Add(histo);
				if (be_verbose) cerr << "[Add MC]   Float control: " << mc_jet_origins[origin_n] << "\t" << histo->Integral() << '\t' << mc_jet_origin_ths[origin_n]->Integral() << endl;
				}
			else
				{
				mc_jet_origin_ths[origin_n] = (TH3D*) histo->Clone();
				if (be_verbose) cerr << "[Clone MC] Float control: " << mc_jet_origins[origin_n] << "\t" << histo->Integral() << '\t' << mc_jet_origin_ths[origin_n]->Integral() << endl;
				}
			}
		//}

	if (be_verbose) cerr << "processed dtag" << endl;
	}


if (add_header)
	{
	cout << "selection_name, data_fakes,data_all ";
	for (int i=0; i<n_jet_origins; i++)
		{
		cout << ",mc" << mc_jet_origins[i] << ",mc" << mc_jet_origins[n_jet_origins + i];
		}
	cout << endl;
	}

/*
// NORMALIZE bins per all_jets
// In Data:
cerr << "histo size = " << hs_data[0]->GetSize() << endl;
cerr << "normalizing data" << endl;
for (int i=0; i<=hs_data[0]->GetSize(); i++)
	{
	double fakes = hs_data[0]->GetBinContent(i);
	double all   = hs_data[1]->GetBinContent(i);
	if (all <= 0.)
		{
		if (fakes > 0)
			cerr << "fakes = " << fakes << " when all jets <= 0" << endl;
		continue;
		}
	hs_data[0]->SetBinContent(i, fakes/all);
	hs_data[0]->SetBinContent(i, 1.); //all/all);
	}

cerr << "normalizing MC" << endl;
// And in MC:
for (int j=0; j<n_jet_origins; j++)
	{
	if (be_verbose)
		{
		cerr << mc_jet_origins[j] << endl;
		cerr << '\t' << mc_jet_origin_ths[j]->Integral() << '\t' << mc_jet_origin_ths[n_jet_origins + j]->Integral() << endl;
		}
	for (int i=0; i<=mc_jet_origin_ths[j]->GetSize(); i++)
		{
		double fakes = mc_jet_origin_ths[j]->GetBinContent(i);
		double all   = mc_jet_origin_ths[n_jet_origins + j]->GetBinContent(i);
		if (all <= 0.)
			{
			if (fakes > 0)
				cerr << "fakes = " << fakes << " when all jets <= 0" << endl;
			continue;
			}
		mc_jet_origin_ths[j]->SetBinContent(i, fakes/all);
		mc_jet_origin_ths[n_jet_origins + j]->SetBinContent(i, 1.);
		}
	if (be_verbose)
		cerr << '\t' << mc_jet_origin_ths[j]->Integral() << '\t' << mc_jet_origin_ths[n_jet_origins + j]->Integral() << endl;
	}
*/

// OUTPUT
// distrname, Data_fakes, Data_all, MC_or1_fakes, MC_or1_all, ...
cout << distr_selection << ',' << hs_data[0]->Integral() << ',' << hs_data[1]->Integral();
for (int i=0; i<n_jet_origins; i++)
	{
	cout << ',' << mc_jet_origin_ths[i]->Integral() << ',' << mc_jet_origin_ths[n_jet_origins + i]->Integral();
	}
cout << endl;

return 0;
}


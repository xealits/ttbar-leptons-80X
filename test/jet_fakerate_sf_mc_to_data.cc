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

#define INPUT_DTAGS_START 8

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = "[--verbose] [--normalize] lumi distr rebin_factor x_axis_min_range x_axis_max_range name_tag dir dtags";
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
//TString projection(argv[input_starts + 3]); // only xy projection (pt-eta)
Int_t rebin_factor(atoi(argv[input_starts + 3]));

double x_axis_min_range = atof(argv[input_starts + 4]);
double x_axis_max_range = atof(argv[input_starts + 5]);
TString name_tag(argv[input_starts + 6]);
TString dir(argv[input_starts + 7]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

/*
if (projection != TString("x") && projection != TString("y") && projection != TString("z"))
	{
	printf("UNKNOWN PROJECTION GIVEN\n");
	printf("supported: x, y, z\n");
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 2;
	}
*/

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

/*
struct tau_jets {
	TH1D* tau_jets;
	TH1D* jets;
};
*/ // nope

TH2D* hs_data[2] = {NULL, NULL};

TH2D* mc_all_jets = NULL;
TH2D* mc_all_taus = NULL;
std::map<TString, TH2D*> nicknamed_mc_tau_jets;

/*
// different jet origin histos:
unsigned int n_jet_origins = 5;
vector<string> mc_jet_origins = {"_tau_jets_distr_o", "_tau_jets_distr_t", "_tau_jets_distr_b", "_tau_jets_distr_q", "_tau_jets_distr_g",
	"_jets_distr_o", "_jets_distr_t", "_jets_distr_b", "_jets_distr_q", "_jets_distr_g"};
//HLTjetmu_qcd_tau_jets_distr_q
vector<TH1D*> mc_jet_origin_ths = {NULL, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL, NULL};
*/


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
			TH2D* histo = (TH2D*) ((TH3D*) file->Get(distro_name))->Project3D("xy");
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
				hs_data[i] = (TH2D*) histo->Clone();
				}
			else
				{
				if (be_verbose) cout << "add histo to data histo" << endl;
				hs_data[i]->Add(histo);
				}
			}
		}

	else // dtag doesn't contain Data = it's MC
		{
		//string jet_origin = mc_jet_origins[origin_n];
		//TString distro_name = distr_selection + jet_origin;
		//distro_names.push_back(distr_selection + string("_tau_jets_distr"));
		//distro_names.push_back(distr_selection + string("_jets_distr"));
		TString tau_jets_name = distr_selection + string("_tau_jets_distr");
		TString     jets_name = distr_selection + string("_jets_distr");

		/* for this dtag
		 * 
		 * MC xsec ratio
		 * then
		 * get jets and tau_jets histograms (if they are present)
		 * scale their bins to width and to xsec ratio
		 * add jets to MC_all_jets histogram
		 * and for tau_jets one needs nickname & colour and to the nicknamed map
		 */

		// different weightflow MC normalization in ttbar and jet fake rates (TODO: need to normalize it to bin10 or something everywhere)
		TH1D * weightflow;
		// actually, only 1 number will be needed:
		double normal_initial_weight = 0;
		if (file->GetListOfKeys()->Contains("weightflow_elel_NOMINAL"))
			{
			weightflow = (TH1D*) file->Get("weightflow_elel_NOMINAL");
			normal_initial_weight = weightflow->GetBinContent(11);
			}
		// TODO: 1 additional just "weightflow" in TTbar for MC weighting, lumi etc -- for MC ratio
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
		// not needed:
		//weightflows.push_back(weightflow);
		//weightflows.back()->Print();
		// --- using only the initial weight

		if (be_verbose) cout << "got weightflow init" << endl;

		// MC ratio for this dtag:
		Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
		//histo->Scale(ratio);
		if (be_verbose) cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		//if (be_verbose) histo->Print();

		// tau_jets -> bin width scale -> MC ratio -> color, marker -> nickames map
		if (file->GetListOfKeys()->Contains(tau_jets_name))
			{
			// get the histogram's projection
			TH2D* histo = (TH2D*) ((TH3D*) file->Get(tau_jets_name))->Project3D("xy");
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

			// Scale to MC ratio
			if (be_verbose) histo->Print();
			histo->Scale(ratio);
			if (be_verbose) histo->Print();

			// colour and nick
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
			TString nick = nick_colour.first;
			Color_t col = nick_colour.second;

			histo->SetFillColor( col );
			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor(col);

			if (mc_all_taus)
				mc_all_taus->Add(histo);
			else
				mc_all_taus = (TH2D*) histo->Clone();

			if (nicknamed_mc_tau_jets.find(nick) == nicknamed_mc_tau_jets.end())
				nicknamed_mc_tau_jets[nick] = (TH2D*) histo->Clone();
				//nicknamed_mc_histos[nick] = new TH1D(nick, "");
			else
				nicknamed_mc_tau_jets[nick]->Add(histo);
			//hs->Add(histos.back(), "HIST");
			}

		// jets -> bin width scale -> MC ratio -> add to jets histogram
		if (file->GetListOfKeys()->Contains(jets_name))
			{
			// get the histogram's projection
			TH2D* histo = (TH2D*) ((TH3D*) file->Get(jets_name))->Project3D("xy");
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

			// Scale to MC ratio
			if (be_verbose) histo->Print();
			histo->Scale(ratio);
			if (be_verbose) histo->Print();

			/*
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
			TString nick = nick_colour.first;
			Color_t col = nick_colour.second;

			histo->SetFillColor( col );
			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor(col);
			*/

			if (mc_all_jets)
				mc_all_jets->Add(histo);
				//nicknamed_mc_histos[nick] = new TH1D(nick, "");
			else
				mc_all_jets = (TH2D*) histo->Clone();
			//hs->Add(histos.back(), "HIST");
			}
		}

	if (be_verbose) cout << "processed dtag" << endl;
	}

double integral_data = hs_data[1]->Integral();
double integral_MC_taus  = 0.;
double integral_MC_jets  = 0.;
cout << "data    integral = " << "(" << hs_data[0]->Integral() << ")" << integral_data << endl;
//cout << "MC taus integral = " << integral_MC_taus   << endl;
cout << "MC jets integral = " << mc_all_jets->Integral() << endl;
cout << "MC taus integral = " << mc_all_taus->Integral() << endl;
//double ratio = integral_data / integral_MC;
//cout << "ratio = " << ratio << endl;

// find the fake rate in data:
hs_data[0]->Sumw2();
hs_data[0]->SetStats(0);      // No statistics on lower plot
hs_data[0]->Divide(hs_data[1]);

// fake rate in MC:
mc_all_taus->Sumw2();
mc_all_taus->SetStats(0);
mc_all_taus->Divide(mc_all_jets);

// SF of MC/Data:
mc_all_taus->Divide(hs_data[0]);

TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("$p_{T}$");
h->SetYTitle("$\\eta$");
cst->Modified();

mc_all_taus->Draw("COLOR");

cout << "drawing" << endl;

bool set_logy = false; // in case it becomes main argument

if (set_logy)
	cst->SetLogy();

//leg->Draw();

cst->SaveAs( dir + "/jobsums/" + distr_selection + "_MCtoDataFakeRateSF_pt-eta_" + name_tag + (normalize_MC ? "_normalized" : "") + (set_logy? "_logy" : "") + ".png" );

return 0;
}


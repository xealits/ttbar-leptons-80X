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
#include "TMath.h"
#include "TPaveText.h"

#include "TStyle.h"

#include <map>
#include <string>
#include <vector>

#include "dtag_xsecs.h"
#include "elmu_fakerates_dtag_xsecs.h"


#define INPUT_DTAGS_START 13

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[256] = " [--verbose] [--normalize] addPU addJER addJES lumi distr suffix rebin_factor x_axis_min_range x_axis_max_range name_tag distr_name dir dtags";
if (argc < INPUT_DTAGS_START)
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

unsigned int i = 1;

bool add_sys_pu  = (TString(argv[input_starts + i++]) == TString("T"));
bool add_sys_jer = (TString(argv[input_starts + i++]) == TString("T"));
bool add_sys_jes = (TString(argv[input_starts + i++]) == TString("T"));

double lumi = atof(argv[input_starts + i++]);
TString distr_selection(argv[input_starts + i++]);
string suffix(argv[input_starts + i++]);
Int_t rebin_factor(atoi(argv[input_starts + i++]));

double x_axis_min_range = atof(argv[input_starts + i++]);
double x_axis_max_range = atof(argv[input_starts + i++]);
TString name_tag(argv[input_starts + i++]);
TString distr_name(argv[input_starts + i++]);
TString dir(argv[input_starts + i++]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

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

TH1D* hs_data[2] = {NULL, NULL};

TH1D* mc_all_jets = NULL;
std::map<TString, TH1D*> nicknamed_mc_tau_jets;

struct hs_mc_sys_shifts {
	TH1D* PU_UP = NULL;
	TH1D* PU_DOWN = NULL;
	TH1D* JER_UP = NULL;
	TH1D* JER_DOWN = NULL;
	TH1D* JES_UP = NULL;
	TH1D* JES_DOWN = NULL;
};

struct hs_mc_sys_shifts mc_all_jets_sys;
struct hs_mc_sys_shifts mc_all_tau_jets_sys;
TH1D* mc_all_tau_jets = NULL;
//std::map<TString, struct hs_mc_sys_shifts> nicknamed_mc_tau_jets_sys;

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
TLegend* leg = new TLegend(0.75, 0.12, 0.9, 0.45);

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
		distro_names.push_back(distr_selection + string("_tau_jets_distr") + suffix);
		distro_names.push_back(distr_selection + string("_jets_distr") + suffix);

		for (int i=0; i<distro_names.size(); i++)
			{
			TString distro_name = distro_names[i];

			if (!file->GetListOfKeys()->Contains(distro_name))
				{
				if (be_verbose) cout << "no " << distro_name << endl;
				continue;
				}

			// get the histogram's projection
			//TH1D* histo = (TH1D*) ((TH3D*) file->Get(distro_name))->Project3D(projection);
			TH1D* histo = (TH1D*) file->Get(distro_name);
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
		TString tau_jets_name = distr_selection + string("_tau_jets_distr") + suffix;
		TString     jets_name = distr_selection + string("_jets_distr") + suffix;

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
			//TH1D* histo = (TH1D*) ((TH3D*) file->Get(tau_jets_name))->Project3D(projection);
			TH1D* histo = (TH1D*) file->Get(tau_jets_name);
			if (be_verbose) cout << "tau jets" << endl;

			/*
			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double error   = histo->GetBinError(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				histo->SetBinError(i, error/width);
				}
			*/

			// Scale to MC ratio
			if (be_verbose) histo->Print();
			histo->Scale(ratio);
			if (be_verbose) histo->Print();

			// colour and nick
			std::pair<TString, Color_t> nick_colour = dtag_nick_colour_elmufakerates(dtag);
			TString nick = nick_colour.first;
			if (nick == TString("tt_{other}")) cout << "tt_{other}     " << dtag << endl;
			Color_t col = nick_colour.second;

			histo->SetFillColor( col );
			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor(col);

			if (nicknamed_mc_tau_jets.find(nick) == nicknamed_mc_tau_jets.end())
				nicknamed_mc_tau_jets[nick] = (TH1D*) histo->Clone();
				//nicknamed_mc_histos[nick] = new TH1D(nick, "");
			else
				nicknamed_mc_tau_jets[nick]->Add(histo);
			//hs->Add(histos.back(), "HIST");

			// add to all taus distr:
			if (mc_all_tau_jets)
				mc_all_tau_jets->Add(histo);
			else
				mc_all_tau_jets = (TH1D*) histo->Clone();

			// also pull all systematic shifts for this micknamed tau_jets
			TString distr_PU_UP = tau_jets_name;
			TString distr_PU_DOWN = tau_jets_name;
			TString distr_JER_UP = tau_jets_name;
			TString distr_JER_DOWN = tau_jets_name;
			TString distr_JES_UP = tau_jets_name;
			TString distr_JES_DOWN = tau_jets_name;
			distr_PU_UP    .ReplaceAll("_vtightTaus", "_PU_UP_vtightTaus"); distr_PU_DOWN  .ReplaceAll("_vtightTaus", "_PU_DOWN_vtightTaus"); distr_JER_UP   .ReplaceAll("_vtightTaus", "_JER_UP_vtightTaus"); distr_JER_DOWN .ReplaceAll("_vtightTaus", "_JER_DOWN_vtightTaus"); distr_JES_UP   .ReplaceAll("_vtightTaus", "_JES_UP_vtightTaus"); distr_JES_DOWN .ReplaceAll("_vtightTaus", "_JES_DOWN_vtightTaus");
			distr_PU_UP    .ReplaceAll("_tightTaus", "_PU_UP_tightTaus"); distr_PU_DOWN  .ReplaceAll("_tightTaus", "_PU_DOWN_tightTaus"); distr_JER_UP   .ReplaceAll("_tightTaus", "_JER_UP_tightTaus"); distr_JER_DOWN .ReplaceAll("_tightTaus", "_JER_DOWN_tightTaus"); distr_JES_UP   .ReplaceAll("_tightTaus", "_JES_UP_tightTaus"); distr_JES_DOWN .ReplaceAll("_tightTaus", "_JES_DOWN_tightTaus");
			distr_PU_UP    .ReplaceAll("_mediumTaus", "_PU_UP_mediumTaus"); distr_PU_DOWN  .ReplaceAll("_mediumTaus", "_PU_DOWN_mediumTaus"); distr_JER_UP   .ReplaceAll("_mediumTaus", "_JER_UP_mediumTaus"); distr_JER_DOWN .ReplaceAll("_mediumTaus", "_JER_DOWN_mediumTaus"); distr_JES_UP   .ReplaceAll("_mediumTaus", "_JES_UP_mediumTaus"); distr_JES_DOWN .ReplaceAll("_mediumTaus", "_JES_DOWN_mediumTaus");
			distr_PU_UP    .ReplaceAll("_looseTaus", "_PU_UP_looseTaus"); distr_PU_DOWN  .ReplaceAll("_looseTaus", "_PU_DOWN_looseTaus"); distr_JER_UP   .ReplaceAll("_looseTaus", "_JER_UP_looseTaus"); distr_JER_DOWN .ReplaceAll("_looseTaus", "_JER_DOWN_looseTaus"); distr_JES_UP   .ReplaceAll("_looseTaus", "_JES_UP_looseTaus"); distr_JES_DOWN .ReplaceAll("_looseTaus", "_JES_DOWN_looseTaus");
			distr_PU_UP    .ReplaceAll("_vlooseTaus", "_PU_UP_vlooseTaus"); distr_PU_DOWN  .ReplaceAll("_vlooseTaus", "_PU_DOWN_vlooseTaus"); distr_JER_UP   .ReplaceAll("_vlooseTaus", "_JER_UP_vlooseTaus"); distr_JER_DOWN .ReplaceAll("_vlooseTaus", "_JER_DOWN_vlooseTaus"); distr_JES_UP   .ReplaceAll("_vlooseTaus", "_JES_UP_vlooseTaus"); distr_JES_DOWN .ReplaceAll("_vlooseTaus", "_JES_DOWN_vlooseTaus");
			vector<TString*> sys_distrs_names = { &distr_PU_UP, &distr_PU_DOWN, &distr_JER_UP, &distr_JER_DOWN, &distr_JES_UP, &distr_JES_DOWN };
			//vector<TH1D*> sys_distrs = { nicknamed_mc_tau_jets_sys[nick].PU_UP, nicknamed_mc_tau_jets_sys[nick].PU_DOWN,
			//	nicknamed_mc_tau_jets_sys[nick].JER_UP, nicknamed_mc_tau_jets_sys[nick].JER_DOWN,
			//	nicknamed_mc_tau_jets_sys[nick].JES_UP, nicknamed_mc_tau_jets_sys[nick].JES_DOWN};
			// it should be a simple list of pairs actually
			vector<TH1D**> sys_distrs = { &mc_all_tau_jets_sys.PU_UP, &mc_all_tau_jets_sys.PU_DOWN,
				&mc_all_tau_jets_sys.JER_UP, &mc_all_tau_jets_sys.JER_DOWN,
				&mc_all_tau_jets_sys.JES_UP, &mc_all_tau_jets_sys.JES_DOWN};

			for (int i=0; i< sys_distrs_names.size(); i++)
				{
				TString& distr = *sys_distrs_names[i];
				//cout << "sys:\t" << distr << '\t';
				if (!file->GetListOfKeys()->Contains(distr))
					{
					cout << "no " << distr << endl;
					continue;
					}

				TH1D* histo = (TH1D*) file->Get(distr);

				histo->Scale(ratio);
				if (rebin_factor != 1) histo->Rebin(rebin_factor); // should rebin the new histo inplace
				// and normalize per bin-width:
				/* no need for fake rate
				for (Int_t i=0; i<=histo->GetSize(); i++)
					{
					//yAxis->GetBinLowEdge(3)
					double content = histo->GetBinContent(i);
					double error   = histo->GetBinError(i);
					double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
					histo->SetBinContent(i, content/width);
					histo->SetBinError(i, error/width);
					}
				*/

				// the pointers in sys_distr are aligned to names correctly -- get them by index
				if (*sys_distrs[i])
					{
					if (be_verbose) cout << "adding" << '\t';
					//nicknamed_mc_tau_jets_sys[nick].PU_UP->Add(histo);
					(*sys_distrs[i])->Add(histo);
					if (be_verbose) cout << "added" << endl;
					}
				else
					{
					if (be_verbose) cout << "clonning" << '\t';
					//hs_mc_sys[i] = (TH1D*) histo->Clone();
					(*sys_distrs[i]) = (TH1D*) histo->Clone();
					if (be_verbose) cout << "cloned" << endl;
					}
				}
			}

		// jets -> bin width scale -> MC ratio -> add to jets histogram
		if (file->GetListOfKeys()->Contains(jets_name))
			{
			// get the histogram's projection
			//TH1D* histo = (TH1D*) ((TH3D*) file->Get(jets_name))->Project3D(projection);
			TH1D* histo = (TH1D*) file->Get(jets_name);
			if (be_verbose) cout << "jets" << endl;

			/*
			// normalize the histo for bin width
			for (Int_t i=0; i<=histo->GetSize(); i++)
				{
				//yAxis->GetBinLowEdge(3)
				double content = histo->GetBinContent(i);
				double error   = histo->GetBinError(i);
				double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
				histo->SetBinContent(i, content/width);
				histo->SetBinError(i, error/width);
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
				mc_all_jets = (TH1D*) histo->Clone();
			//hs->Add(histos.back(), "HIST");

			if (be_verbose) cout << "added to all jets" << endl;

			// and add up the systematic shifts:
			// also pull all systematic shifts for this micknamed jets
			TString distr_PU_UP = jets_name;
			TString distr_PU_DOWN = jets_name;
			TString distr_JER_UP = jets_name;
			TString distr_JER_DOWN = jets_name;
			TString distr_JES_UP = jets_name;
			TString distr_JES_DOWN = jets_name;
			distr_PU_UP    .ReplaceAll("_vtightTaus", "_PU_UP_vtightTaus"); distr_PU_DOWN  .ReplaceAll("_vtightTaus", "_PU_DOWN_vtightTaus"); distr_JER_UP   .ReplaceAll("_vtightTaus", "_JER_UP_vtightTaus"); distr_JER_DOWN .ReplaceAll("_vtightTaus", "_JER_DOWN_vtightTaus"); distr_JES_UP   .ReplaceAll("_vtightTaus", "_JES_UP_vtightTaus"); distr_JES_DOWN .ReplaceAll("_vtightTaus", "_JES_DOWN_vtightTaus");
			distr_PU_UP    .ReplaceAll("_tightTaus", "_PU_UP_tightTaus"); distr_PU_DOWN  .ReplaceAll("_tightTaus", "_PU_DOWN_tightTaus"); distr_JER_UP   .ReplaceAll("_tightTaus", "_JER_UP_tightTaus"); distr_JER_DOWN .ReplaceAll("_tightTaus", "_JER_DOWN_tightTaus"); distr_JES_UP   .ReplaceAll("_tightTaus", "_JES_UP_tightTaus"); distr_JES_DOWN .ReplaceAll("_tightTaus", "_JES_DOWN_tightTaus");
			distr_PU_UP    .ReplaceAll("_mediumTaus", "_PU_UP_mediumTaus"); distr_PU_DOWN  .ReplaceAll("_mediumTaus", "_PU_DOWN_mediumTaus"); distr_JER_UP   .ReplaceAll("_mediumTaus", "_JER_UP_mediumTaus"); distr_JER_DOWN .ReplaceAll("_mediumTaus", "_JER_DOWN_mediumTaus"); distr_JES_UP   .ReplaceAll("_mediumTaus", "_JES_UP_mediumTaus"); distr_JES_DOWN .ReplaceAll("_mediumTaus", "_JES_DOWN_mediumTaus");
			distr_PU_UP    .ReplaceAll("_looseTaus", "_PU_UP_looseTaus"); distr_PU_DOWN  .ReplaceAll("_looseTaus", "_PU_DOWN_looseTaus"); distr_JER_UP   .ReplaceAll("_looseTaus", "_JER_UP_looseTaus"); distr_JER_DOWN .ReplaceAll("_looseTaus", "_JER_DOWN_looseTaus"); distr_JES_UP   .ReplaceAll("_looseTaus", "_JES_UP_looseTaus"); distr_JES_DOWN .ReplaceAll("_looseTaus", "_JES_DOWN_looseTaus");
			distr_PU_UP    .ReplaceAll("_vlooseTaus", "_PU_UP_vlooseTaus"); distr_PU_DOWN  .ReplaceAll("_vlooseTaus", "_PU_DOWN_vlooseTaus"); distr_JER_UP   .ReplaceAll("_vlooseTaus", "_JER_UP_vlooseTaus"); distr_JER_DOWN .ReplaceAll("_vlooseTaus", "_JER_DOWN_vlooseTaus"); distr_JES_UP   .ReplaceAll("_vlooseTaus", "_JES_UP_vlooseTaus"); distr_JES_DOWN .ReplaceAll("_vlooseTaus", "_JES_DOWN_vlooseTaus");
			vector<TString*> sys_distrs_names = { &distr_PU_UP, &distr_PU_DOWN, &distr_JER_UP, &distr_JER_DOWN, &distr_JES_UP, &distr_JES_DOWN };
			vector<TH1D**> sys_distrs = { &mc_all_jets_sys.PU_UP, &mc_all_jets_sys.PU_DOWN,
				&mc_all_jets_sys.JER_UP, &mc_all_jets_sys.JER_DOWN,
				&mc_all_jets_sys.JES_UP, &mc_all_jets_sys.JES_DOWN};
			// it should be a simple list of pairs actually
			if (be_verbose) cout << "prepared systematics" << endl;

			for (int i=0; i< sys_distrs_names.size(); i++)
				{
				TString& distr = *sys_distrs_names[i];
				if (be_verbose) cout << "sys:\t" << distr << endl;
				if (!file->GetListOfKeys()->Contains(distr))
					{
					cout << "no " << distr << endl;
					continue;
					}

				TH1D* histo = (TH1D*) file->Get(distr);
				if (be_verbose) cout << "got distr" << endl;

				histo->Scale(ratio);
				if (be_verbose) cout << "scaled lumi" << endl;
				if (rebin_factor != 1)
					{
					histo->Rebin(rebin_factor); // should rebin the new histo inplace
					if (be_verbose) cout << "re-binned" << endl;
					}
				// and normalize per bin-width:
				/* no need to do it for fake rate
				for (Int_t i=0; i<=histo->GetSize(); i++)
					{
					//yAxis->GetBinLowEdge(3)
					double content = histo->GetBinContent(i);
					double error   = histo->GetBinError(i);
					double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
					histo->SetBinContent(i, content/width);
					histo->SetBinError(i, error/width);
					}
				*/

				if (be_verbose) cout << "scaled the distr, adding to struct via pointers" << endl;
				// the pointers in sys_distr are aligned to names correctly -- get them by index
				if (*sys_distrs[i])
					{
					if (be_verbose) cout << "adding" << '\t';
					//nicknamed_mc_tau_jets_sys[nick].PU_UP->Add(histo);
					(*sys_distrs[i])->Add(histo);
					if (be_verbose) cout << "added" << endl;
					}
				else
					{
					if (be_verbose) cout << "clonning" << '\t';
					//hs_mc_sys[i] = (TH1D*) histo->Clone();
					(*sys_distrs[i]) = (TH1D*) histo->Clone();
					if (be_verbose) cout << "cloned" << endl;
					}
				}
			}
		}

	if (be_verbose) cout << "processed dtag" << endl;
	}

double integral_data = hs_data[1]->Integral();
double integral_MC_taus  = 0.;
double integral_MC_jets  = 0.;
// INCLUSIVE FAKERATE
cout << "data    integral = " << hs_data[0]->Integral()   << " / " << hs_data[1]->Integral() << " = " << hs_data[0]->Integral() / hs_data[1]->Integral() << endl;
//cout << "MC taus integral = " << integral_MC_taus   << endl;
cout << "MC jets integral = " << mc_all_tau_jets->Integral() << " / " << mc_all_jets->Integral() << " = " << mc_all_tau_jets->Integral() / mc_all_jets->Integral() << endl;
//double ratio = integral_data / integral_MC;
//cout << "ratio = " << ratio << endl;

/*
if (normalize_MC)
	{
	cout << "normalizing MCs" << endl;
	for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
		{
		double integral = mc_jet_origin_ths[origin_n]->Integral();
		double ratio = integral_data / integral;
		mc_jet_origin_ths[origin_n]->Scale(ratio);
		cout << mc_jet_origins[origin_n] << " integral before = " << integral << " after = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
		}
	}
*/


// go through MC tau_jets nickname map:
//  1. divide every tau_jets distr by all_jets
// and build the sum of MC
TH1D    *hs_sum  = NULL;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_tau_jets.begin(); it != nicknamed_mc_tau_jets.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;

	if (be_verbose) cout << nick;

	// build the summ of MC
	if (be_verbose) cout << " summing mc";
	if (hs_sum == NULL)
		{
		hs_sum = (TH1D*) (distr->Clone("brand_new_hs_sum"));
		hs_sum->SetName("hs_sum");
		hs_sum->ResetBit(kCanDelete);
		}
	else
		hs_sum->Add(distr);

	if (be_verbose)
		{
		cout << " adding to mc stack" << endl;
		distr->Print();
		}
	}
// now this summ is not needed, since it's computed in the loop
// leaving it for control
cout << "compare" << endl;
hs_sum->Print();
mc_all_tau_jets->Print();

hs_sum = mc_all_tau_jets;
cout << "switched to the second one" << endl;

/* no need to normalize in case of fake rate
 * it's a function, not count of events in bins
// normalize all histos for bin width:
// hs_sum      (MC sum of fake taus)
// mc_all_jets (MC sum of jets)
// hs_data[0]  (Data fakes)
// hs_data[1]  (Data jets)
TH1D* histos_to_scale[4] = {hs_sum, mc_all_jets, hs_data[0], hs_data[1]};
for (Int_t h=0; h<4; h++)
	{
	TH1D* histo = histos_to_scale[h];
	for (Int_t i=0; i<=histo->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = histo->GetBinContent(i);
		double error   = histo->GetBinError(i);
		double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
		histo->SetBinContent(i, content/width);
		histo->SetBinError(i, error/width);
		}
	}
*/

// add systematic errors to the MC sum histos:
// taus:
for (int i=0; i<=hs_sum->GetSize(); i++)
	{
	double mc_content = hs_sum->GetBinContent(i);
	double mc_error   = hs_sum->GetBinError(i);

	double mc_PU_UP  = mc_all_tau_jets_sys.PU_UP ->GetBinContent(i);
	double mc_JER_UP = mc_all_tau_jets_sys.JER_UP->GetBinContent(i);
	double mc_JES_UP = mc_all_tau_jets_sys.JES_UP->GetBinContent(i);

	double mc_PU_DOWN  = mc_all_tau_jets_sys.PU_DOWN ->GetBinContent(i);
	double mc_JER_DOWN = mc_all_tau_jets_sys.JER_DOWN->GetBinContent(i);
	double mc_JES_DOWN = mc_all_tau_jets_sys.JES_DOWN->GetBinContent(i);

	double sys_PU   = (abs(mc_PU_UP  - mc_content) + abs(mc_PU_DOWN  - mc_content))/2;
	double sys_JER  = (abs(mc_JER_UP - mc_content) + abs(mc_JER_DOWN - mc_content))/2;
	double sys_JES  = (abs(mc_JES_UP - mc_content) + abs(mc_JES_DOWN - mc_content))/2;

	double sys_error = TMath::Sqrt(sys_PU*sys_PU + sys_JER*sys_JER + sys_JES*sys_JES + mc_error*mc_error);
	//cout << '(' << mc_PU_UP << '_' << mc_JER_UP << '_' << mc_JES_UP << '-' << mc_content << ',';
	//cout << mc_PU_DOWN << '_' << mc_JER_DOWN << '_' << mc_JES_DOWN << '-' << mc_content << '=';
	//cout << mc_error << '_' << sys_error << ") ";

	hs_sum->SetBinError(i, sys_error);	
	}
cout << "added tau jets systematics" << endl;

// jets:
for (int i=0; i<=mc_all_jets->GetSize(); i++)
	{
	double mc_content = mc_all_jets->GetBinContent(i);
	double mc_error   = mc_all_jets->GetBinError(i);

	double mc_PU_UP  = mc_all_jets_sys.PU_UP ->GetBinContent(i);
	double mc_JER_UP = mc_all_jets_sys.JER_UP->GetBinContent(i);
	double mc_JES_UP = mc_all_jets_sys.JES_UP->GetBinContent(i);

	double mc_PU_DOWN  = mc_all_jets_sys.PU_DOWN ->GetBinContent(i);
	double mc_JER_DOWN = mc_all_jets_sys.JER_DOWN->GetBinContent(i);
	double mc_JES_DOWN = mc_all_jets_sys.JES_DOWN->GetBinContent(i);

	//double sys_PU   = (abs(hs_mc_PU_UP  ->GetBinContent(i) - mc_content) + abs(hs_mc_PU_DOWN->GetBinContent(i) - mc_content ))/2;
	//double sys_JER  = (abs(hs_mc_JER_UP ->GetBinContent(i) - mc_content) + abs(hs_mc_JER_DOWN->GetBinContent(i) - mc_content))/2;
	//double sys_JES  = (abs(hs_mc_JES_UP ->GetBinContent(i) - mc_content) + abs(hs_mc_JES_DOWN->GetBinContent(i) - mc_content))/2;
	double sys_PU   = (add_sys_pu  ? (abs(mc_PU_UP  - mc_content) + abs(mc_PU_DOWN  - mc_content))/2 : 0);
	double sys_JER  = (add_sys_jer ? (abs(mc_JER_UP - mc_content) + abs(mc_JER_DOWN - mc_content))/2 : 0);
	double sys_JES  = (add_sys_jes ? (abs(mc_JES_UP - mc_content) + abs(mc_JES_DOWN - mc_content))/2 : 0);

	double sys_error = TMath::Sqrt(sys_PU*sys_PU + sys_JER*sys_JER + sys_JES*sys_JES + mc_error*mc_error);
	//cout << '(' << mc_PU_UP << '_' << mc_JER_UP << '_' << mc_JES_UP << '-' << mc_content << ',';
	//cout << mc_PU_DOWN << '_' << mc_JER_DOWN << '_' << mc_JES_DOWN << '-' << mc_content << '=';
	//cout << mc_error << '_' << sys_error << ") ";

	mc_all_jets->SetBinError(i, sys_error);	
	}
cout << "added all jets systematics" << endl;


// INCLUSIVE FAKERATE divide MC taus by mc jets
cout << "MC integrals: " << hs_sum->Integral() << " / " << mc_all_jets->Integral() << " = " <<  hs_sum->Integral() / mc_all_jets->Integral() << endl;
hs_sum->Divide(mc_all_jets);

// find the fake rate in data:
hs_data[0]->Sumw2();
hs_data[0]->SetStats(0);      // No statistics on lower plot
hs_data[0]->Divide(hs_data[1]);

// build the stack of MC
// go through MC tau_jets nickname map:
//  2. add it to MC histo stack
//  3. add it to the legend
THStack *hs      = new THStack("hs", "");
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_tau_jets.begin(); it != nicknamed_mc_tau_jets.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;

	if (be_verbose)
		{
		cout << nick;
		cout << " adding to mc stack" << endl;
		distr->Print();
		}

	/* no need to scale for binwidth anything
	// now scale the distr for bin-width
	for (Int_t i=0; i<=distr->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = distr->GetBinContent(i);
		double error   = distr->GetBinError(i);
		double width   = distr->GetXaxis()->GetBinUpEdge(i) - distr->GetXaxis()->GetBinLowEdge(i);
		distr->SetBinContent(i, content/width);
		distr->SetBinError  (i, error/width);
		}
	*/

	distr->Divide(mc_all_jets);

	hs->Add(distr, "HIST");
	leg->AddEntry(distr, nick, "F");
	}


//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

cout << "setting title" << endl;

//hs->GetXaxis()->SetTitle(distr);
//cst->SetXaxisTile(distr);
//hs_data->GetXaxis()->SetTitle(distr);


TH1D * hs_data_relative = (TH1D*) hs_data[0]->Clone();
TH1D * hs_sum_relative  = (TH1D*) hs_sum->Clone();

// just MC and Data legend
TLegend* leg2 = new TLegend(0.6, 0.65, 0.89, 0.89);
leg2->SetHeader(name_tag);
TLegendEntry *le_data = (TLegendEntry*) leg2->AddEntry(hs_data[0], "Data", "lpe");
TLegendEntry *le_mc   = (TLegendEntry*) leg2->AddEntry(hs_sum,     "Simulation", "lpe");

TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
//h->SetXTitle("x");
//cst->Update();
//cst->Modified();

gStyle->SetOptStat(0); // removes the stats legend-box from all pads

TPad *pad1 = new TPad("pad1","This is pad1", 0., 0.3, 1., 1.);
TPad *pad2 = new TPad("pad2","This is pad2", 0., 0.05,  1., 0.3);
//pad1->SetFillColor(11);
//pad2->SetFillColor(11);

// Upper and lower plot are joined
pad1->SetBottomMargin(0);
pad2->SetTopMargin(0);

bool set_logy = true;

if (set_logy)
	{
	cout << "setting logy" << endl;
	pad1->SetLogy();
	//gPad->SetLogy();
	}

pad1->Draw();
pad2->Draw();

pad1->cd();

hs_data[0]->GetXaxis()->SetLabelFont(63);
hs_data[0]->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_data[0]->GetYaxis()->SetLabelFont(63);
hs_data[0]->GetYaxis()->SetLabelSize(14); // labels will be 14 pixels

if (be_verbose) cout << "drawing the distr-s" << endl;

/*
hs_data->Draw("e p");
hs->Draw("same");
hs_data->Draw("e p same"); // to draw it _over_ MC

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");

hs->GetXaxis()->SetTitle(distr);
hs_data->SetXTitle(distr);
*/

hs_sum->    GetYaxis()->SetRange(0.000008, 1);
hs_sum->    GetYaxis()->SetRangeUser(0.000008, 1);
hs_data[0]->GetYaxis()->SetRange(0.000008, 1);
hs_data[0]->GetYaxis()->SetRangeUser(0.000008, 1);

hs_data[0]->SetMarkerSize(0.5);
//hs->Draw("same"); // the stack
// also draw the sum of the stack and its' errors:
//((TObjArray*) hs->GetStack())->Last()->Draw("e4 same"); // then, how to set styles with the hatched error bars?
hs_data[0]->GetXaxis()->SetLabelOffset(999);
hs_data[0]->GetXaxis()->SetLabelSize(0);

hs_data[0]->Draw("p"); // drawing data-MC-data to have them both in the range of the plot

if (hs_sum != NULL)
	{
	hs_sum->Print();
	hs_sum->GetXaxis()->SetLabelOffset(999);
	hs_sum->GetXaxis()->SetLabelSize(0);
	hs_sum->SetFillStyle(3004);
	hs_sum->SetFillColor(1);
	//hs_sum->SetMarkerColorAlpha(0, 0.1);
	hs_sum->SetMarkerStyle(25);
	hs_sum->SetMarkerColor(kRed);
	hs_sum->SetLineColor(kRed);
	hs_sum->Draw("e same"); // the errors on the stack
	}
else
	cout << "NO HS_SUM!!!!!!!!!!!!!!" << endl;

hs_data[0]->Draw("e p same"); // the data with the errors

// histo->SetTitle("boxtitle;x axis title [unit];y axis title [unit]")
cout << "setting the stack title" << endl;

TString distr = distr_selection + TString((const char*) suffix.c_str());

//hs->GetXaxis()->SetTitle(distr);
cout << "done setting the stack title" << endl;
hs_data[0]->GetYaxis()->SetTitleFont(63);
hs_data[0]->GetYaxis()->SetTitleSize(24);
hs_data[0]->GetYaxis()->SetTitleOffset(1);
hs_data[0]->SetYTitle("Fake rate");
hs_data[0]->SetXTitle("");
hs_sum->SetXTitle("");

// label texts
TPaveText* left_title = new TPaveText(0.1, 0.9, 0.4, 0.95, "brNDC");
left_title->AddText("CMS preliminary at 13 TeV");
left_title->SetTextFont(1);
left_title->SetFillColor(0);
cout << "drawing left title" << endl;
left_title->Draw("same");

TPaveText* right_title = new TPaveText(0.75, 0.9, 0.9, 0.95, "brNDC");
TString s_title(""); s_title.Form("L = %.1f fb^{-1}", lumi/1000);
right_title->AddText(s_title);
right_title->SetTextFont(132);
right_title->SetFillColor(0);
cout << "drawing right title" << endl;
right_title->Draw("same");

leg2->SetBorderSize(0);
leg2->Draw();

// THE RATIO PLOT
pad2->cd();
//cst->cd(2);

hs_data_relative->SetMarkerSize(0.5);
hs_sum_relative->SetMarkerStyle(1);
hs_sum_relative->SetFillColor(kRed);
hs_sum_relative->SetFillStyle(3004);

hs_data_relative->Print();
hs_sum_relative->Print();

hs_data_relative->SetXTitle("");
hs_sum_relative->SetXTitle("");

hs_data_relative->SetMarkerColor(1);

// to have the same font size on both pads:
hs_data_relative->GetXaxis()->SetLabelFont(63);
hs_data_relative->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_data_relative->GetYaxis()->SetLabelFont(63);
hs_data_relative->GetYaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_sum_relative->GetXaxis()->SetLabelFont(63);
hs_sum_relative->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_sum_relative->GetYaxis()->SetLabelFont(63);
hs_sum_relative->GetYaxis()->SetLabelSize(14); // labels will be 14 pixels

cout << "distr_name !!!!!!!!!!!!!!!!! " << distr_name << endl;
hs_data_relative->SetXTitle(distr_name);
hs_sum_relative->SetXTitle(distr_name);
hs_data_relative->SetYTitle("obs/exp");
hs_sum_relative->SetYTitle("obs/exp");

hs_data_relative->GetXaxis()->SetTitleFont(63);
hs_data_relative->GetXaxis()->SetTitleSize(24);
hs_sum_relative ->GetXaxis()->SetTitleFont(63);
hs_sum_relative ->GetXaxis()->SetTitleSize(24);
hs_data_relative->GetYaxis()->SetTitleFont(63);
hs_data_relative->GetYaxis()->SetTitleSize(24);
hs_sum_relative ->GetYaxis()->SetTitleFont(63);
hs_sum_relative ->GetYaxis()->SetTitleSize(24);

hs_data_relative->GetYaxis()->SetTitleOffset(1);
hs_sum_relative ->GetYaxis()->SetTitleOffset(1);
hs_data_relative->GetXaxis()->SetTitleOffset(3.5);
hs_sum_relative ->GetXaxis()->SetTitleOffset(3.5);


for (int i=0; i<=hs_sum_relative->GetSize(); i++)
	{
	double mc_content = hs_sum_relative->GetBinContent(i);
	double mc_error   = hs_sum_relative->GetBinError(i);
	//double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);

	double data_content = hs_data_relative->GetBinContent(i);
	double data_error   = hs_data_relative->GetBinError(i);
	//double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
	if (mc_content > 0)
		{
		hs_sum_relative->SetBinContent(i, 1);
		hs_sum_relative->SetBinError(i, mc_error/mc_content);	
		hs_data_relative->SetBinContent(i, data_content/mc_content);
		hs_data_relative->SetBinError(i, data_error/mc_content);
		}
	else
		{
		hs_sum_relative->SetBinContent(i, 1);
		hs_sum_relative->SetBinError(i, mc_error);
		hs_data_relative->SetBinContent(i, 1);
		hs_data_relative->SetBinError(i, data_error);
		}
	}

hs_sum_relative->SetStats(false);
hs_data_relative->SetStats(false);
hs_sum_relative->GetYaxis()->SetRange(0, 3.4);
hs_sum_relative->GetYaxis()->SetRangeUser(0, 3.4);
hs_data_relative->GetYaxis()->SetRange(0, 3.4);
hs_data_relative->GetYaxis()->SetRangeUser(0, 3.4);

/*
if (xlims_set)
	{
	hs_sum_relative->GetXaxis()->SetRange(xmin, xmax);
	hs_sum_relative->GetXaxis()->SetRangeUser(xmin, xmax);
	//hs_data_relative->GetXaxis()->SetRange(xmin, xmax);
	//hs_data_relative->GetXaxis()->SetRangeUser(xmin, xmax);
	}
*/

hs_sum_relative->Draw("e2");
hs_data_relative->Draw("e p same");


/*
cout << "drawing" << endl;

//hs_data->Draw("e"); // to draw it _over_ MC

hs_data[0]->GetYaxis()->SetRange(0.0001, 1);
hs_data[0]->GetYaxis()->SetRangeUser(0.0001, 1);
hs_data[0]->GetXaxis()->SetRange(x_axis_min_range, x_axis_max_range);
hs_data[0]->GetXaxis()->SetRangeUser(x_axis_min_range, x_axis_max_range);

hs->GetYaxis()->SetRange(0.0001, 1);
hs->GetYaxis()->SetRangeUser(0.0001, 1);
hs->GetXaxis()->SetRange(x_axis_min_range, x_axis_max_range);
hs->GetXaxis()->SetRangeUser(x_axis_min_range, x_axis_max_range);

//mc_jet_origin_ths[0]->GetXaxis()->SetRange(x_axis_min_range, x_axis_max_range);
//mc_jet_origin_ths[0]->GetXaxis()->SetRangeUser(x_axis_min_range, x_axis_max_range);
//h3->GetXaxis()->SetRange(x_axis_min_range, x_axis_max_range);
//h3->GetXaxis()->SetRangeUser(x_axis_min_range, x_axis_max_range);

//mc_jet_origin_ths[0]->GetYaxis()->SetRange(0.0001, 1);
//mc_jet_origin_ths[0]->GetYaxis()->SetRangeUser(0.0001, 1);
//h3->GetYaxis()->SetRange(0.0001, 1); // ranges from analysis note CMS AN-2012/489
//h3->GetYaxis()->SetRangeUser(0.0001, 1); // ranges from analysis note CMS AN-2012/489

bool set_logy = true; // in case it becomes main argument

if (set_logy)
	cst->SetLogy();

hs_data[0]->Draw("e p"); // to pass the title/axes settings to the canvas/plot/pad/etc root-shmuz
hs->Draw("same"); // draw the stack
hs_data[0]->Draw("e same p"); // to overlay MC

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");
*/


//cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr_selection + "_MCFakeRates_SYS_" + suffix + "_" + name_tag + (add_sys_pu ? "_wPU" : "") + (add_sys_jer ? "_wJER" : "") + (add_sys_jes ? "_wJES" : "") + (normalize_MC ? "_normalized" : "") + (set_logy? "_logy" : "") + ".png" );

return 0;
}


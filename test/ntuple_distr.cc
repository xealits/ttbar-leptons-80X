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

#include "TStyle.h"

#include <map>
#include <string>
#include <vector>

#include "dtag_xsecs.h"

#define INPUT_DTAGS_START 10
#define NT_OUTPUT_TTREE_NAME "reduced_ttree"

double pileup_ratio[] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0};

double pu_weight(Int_t nvtx) {
   return pileup_ratio[nvtx];
}

using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = "[--verbose] [--normalize] unstack lumi distr distr_cond range out_name distr_name dir dtags";
if (argc < INPUT_DTAGS_START)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();
gROOT->ProcessLine(".L pu_weight.C+");

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

bool set_logy = false;
TString logy_inp(argv[input_starts + 1]);
if (logy_inp == TString("T") || logy_inp == TString("Y"))
	{
	set_logy = true;
	}

bool unstack = false;
TString unstack_inp(argv[input_starts + 2]);
if (unstack_inp == TString("T") || unstack_inp == TString("Y"))
	{
	unstack = true;
	}

double lumi = atof(argv[input_starts + 3]);
TString distr(argv[input_starts + 4]);
TString distr_condition(argv[input_starts + 5]);
TString range(argv[input_starts + 6]);
string out_name(argv[input_starts + 7]);

TString distr_name(argv[input_starts + 8]);
TString dir(argv[input_starts + 9]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

cout << lumi  << endl;
cout << distr << endl;
cout << range << endl;
cout << dir   << endl;
cout << dtag1 << endl;


std::vector < TString > dtags;
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

std::map<TString, TH1D*> nicknamed_mc_distrs;
TH1D* mc_sum = NULL;

//THStack *hs = new THStack("hs","Stacked 1D histograms");
THStack *hs      = new THStack("hs", "");
TH1D    *hs_data = NULL;

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
TLegend* leg = new TLegend(0.7, 0.7, 0.89, 0.89);

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
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	TFile* file = TFile::Open(the_file);
	//dtags.push_back(5);

	bool isData = dtag.Contains("Data");

	/*
	 * data containd distr + _jet_distr
	 * MC: + jet_distr_g/q/b/t/o
	 */
	/*
	if (!file->GetListOfKeys()->Contains(distr))
		{
		cout << "no " << distr << endl;
		continue;
		}
	*/

	TTree* output_ttree = (TTree*) file->Get(NT_OUTPUT_TTREE_NAME); // hardcoded output name
	if (be_verbose) cout << "got ttree, drawing:" << endl;

	//TH1D* histo = new TH1D("myhist"+dtag, "title", Nbins, Xlow, Xup);
	//if (be_verbose) cout << "drawing " << distr + ">>myhist" << endl;
	//output_ttree->Draw(distr + ">>myhist"+dtag, distr_condition);
	if (isData)
		output_ttree->Draw(distr + ">>h" + range, distr_condition);
	else
		// add weights for MC
		output_ttree->Draw(distr + ">>h" + range, "(pu_weight(nvtx_gen))*(aMCatNLO_weight > 0? 1 : 0)* (" + distr_condition + ")");
	TH1D* histo = (TH1D*) output_ttree->GetHistogram();

	// hack to merge HLTmu & HLTjetmu
	if (be_verbose) cout << "got histogram, scaling to lumi, adding to stack" << endl;

	Double_t ratio = 1;
	if (!isData)
		{
		TH1D* weightflow = NULL;
		if (file->GetListOfKeys()->Contains("eventflow"))
			{
			weightflow = (TH1D*) file->Get("eventflow");
			}
		else if (file->GetListOfKeys()->Contains("weightflow_el_NOMINAL"))
			{
			weightflow = (TH1D*) file->Get("weightflow_el_NOMINAL");
			}
		else {
			cout << "NO WEIGHTFLOW: " << dtag << endl;
			}
		ratio = lumi * xsecs[dtag] / weightflow->GetBinContent(11);
		histo->Scale(ratio);
		cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		if (be_verbose) histo->Print();
		}
	//if (rebin_factor != 1) histo->Rebin(rebin_factor); // should rebin the new histo inplace
	// and normalize per bin-width:
	for (Int_t i=0; i<=histo->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = histo->GetBinContent(i);
		double error   = histo->GetBinError(i);
		double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
		histo->SetBinContent(i, content/width);
		histo->SetBinError(i, error/width);
		}

	if (be_verbose) histo->Print();
	if (isData)
		{
		if (be_verbose) cout << "summing data-stack" << endl;
		histo->SetMarkerStyle(9);

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
		{
		if (be_verbose) cout << "summing mc-stack" << endl;
		std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
		TString nick = nick_colour.first;
		Color_t col = nick_colour.second;

		if (!unstack)
			{
			histo->SetFillColor( col );
			histo->SetMarkerStyle(20);
			histo->SetLineStyle(0);
			histo->SetMarkerColor(i);
			}
		else
			{
			histo->SetLineColor(col);
			histo->SetLineStyle(1);
			}

		//std::map<TString, TH1D *> nicknamed_mc_distrs;
		if (nicknamed_mc_distrs.find(nick) == nicknamed_mc_distrs.end())
			nicknamed_mc_distrs[nick] = (TH1D*) histo->Clone();
			//nicknamed_mc_distrs[nick] = new TH1D(nick, "");
		else
			nicknamed_mc_distrs[nick]->Add(histo);

		if (be_verbose) cout << "summing mc,";
		// and add to the sum:
		if (mc_sum == NULL)
			{
			if (be_verbose) cout << " clone,";
			mc_sum = (TH1D*) histo->Clone();
			}
		else
			{
			if (be_verbose) cout << " sum,";
			mc_sum->Add(histo);
			}
		if (be_verbose) cout << " done." << endl;
		}

	if (be_verbose) cout << "processed dtag" << endl;
	}

double integral_data = hs_data->Integral();
double integral_MC_taus  = 0.;
double integral_MC_jets  = 0.;
// INCLUSIVE event yield
cout << "------------------ inclusive event yield ------------------" << endl;
cout << "data    integral: " << integral_data << endl;
//cout << "MC taus integral = " << integral_MC_taus   << endl;
cout << "MC jets integral = " << mc_sum->Integral() << endl;
//double ratio = integral_data / integral_MC;
//cout << "ratio = " << ratio << endl;


// go through MC tau_jets nickname map:
//  1. divide every tau_jets distr by all_jets
// and build the sum of MC
TH1D    *hs_sum  = NULL;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
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
mc_sum->Print();

hs_sum = mc_sum;
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

hs_data->Sumw2();
hs_data->SetStats(0);      // No statistics on lower plot

// build the stack of MC
// go through MC tau_jets nickname map:
//  2. add it to MC histo stack
//  3. add it to the legend
//THStack *hs      = new THStack("hs", "");
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
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

if (unstack)
	{
	hs_data->Scale(1/hs_data->Integral());
	hs_sum-> Scale(1/hs_sum->Integral());
	for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
		{
		TString nick = it->first;
		TH1D * distr = it->second;

		distr->Scale(1/distr->Integral());
		}
	}

TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("x");
//cst->Update();
//cst->Modified();

gStyle->SetOptStat(0); // removes the stats legend-box from all pads

TPad *pad1 = new TPad("pad1","This is pad1", 0., 0.2, 1., 1.);
TPad *pad2 = new TPad("pad2","This is pad2", 0., 0.,  1., 0.2);
//pad1->SetFillColor(11);
//pad2->SetFillColor(11);


if (set_logy)
	{
	cout << "setting logy" << endl;
	pad1->SetLogy();
	//gPad->SetLogy();
	}

pad1->Draw();
pad2->Draw();

pad1->cd();

hs_data->GetXaxis()->SetLabelFont(63);
hs_data->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_data->GetYaxis()->SetLabelFont(63);
hs_data->GetYaxis()->SetLabelSize(14); // labels will be 14 pixels

if (be_verbose) cout << "drawing the distr-s" << endl;

hs_data->Draw("p"); // drawing data-MC-data to have them both in the range of the plot

if (!unstack)
	{
	hs->Draw("same"); // the stack
	if (hs_sum != NULL)
		{
		hs_sum->Print();
		hs_sum->SetFillStyle(3004);
		hs_sum->SetFillColor(1);
		hs_sum->SetMarkerColorAlpha(0, 0.1);
		hs_sum->SetMarkerStyle(1);
		hs_sum->SetMarkerColor(0);
		hs_sum->Draw("e2 same"); // the errors on the stack
		}
	else
		cout << "NO HS_SUM!!!!!!!!!!!!!!" << endl;
	}
else
	{
	for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
		{
		TString nick = it->first;
		TH1D * distr = it->second;
		distr->Draw("Lhist same");
		}
	}

hs_data->Draw("e p same"); // the data with the errors

// histo->SetTitle("boxtitle;x axis title [unit];y axis title [unit]")
cout << "setting the stack title" << endl;

//TString distr_out = distr + TString((const char*) suffix.c_str());

//hs->GetXaxis()->SetTitle(distr);
cout << "done setting the stack title" << endl;
hs_data->SetXTitle(distr_name);
hs_sum->SetXTitle(distr_name);

leg->SetBorderSize(0);
leg->Draw();

// THE RATIO PLOT
pad2->cd();
//cst->cd(2);

TH1D * hs_data_relative = (TH1D*) hs_data->Clone();
TH1D * hs_sum_relative  = (TH1D*) hs_sum->Clone();

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
hs_sum_relative->GetYaxis()->SetRange(0.75, 1.25);
hs_sum_relative->GetYaxis()->SetRangeUser(0.75, 1.25);
hs_data_relative->GetYaxis()->SetRange(0.75, 1.25);
hs_data_relative->GetYaxis()->SetRangeUser(0.75, 1.25);

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


//cst->Modified();

cst->SaveAs( dir + "/jobsums/" + out_name + "_QuickNtupleDistr_" + (normalize_MC ? "_normalized" : "") + (set_logy? "_logy" : "") + (unstack? "_unstacked" : "") + ".png" );

return 0;
}


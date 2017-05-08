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
char usage_string[128] = "[--verbose] [--normalize] lumi distr projection rebin_factor dir dtags";
if (argc < 5)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

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

double lumi = atof(argv[input_starts + 1]);
TString distr(argv[input_starts + 2]);
TString projection(argv[input_starts + 3]);
Int_t rebin_factor(atoi(argv[input_starts + 4]));
TString dir(argv[input_starts + 5]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

if (projection != TString("x") && projection != TString("y") && projection != TString("z"))
	{
	printf("UNKNOWN PROJECTION GIVEN\n");
	printf("supported: x, y, z\n");
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 2;
	}

if (be_verbose)
	{
	cout << lumi  << endl;
	cout << distr << endl;
	cout << dir   << endl;
	cout << dtag1 << endl;
	}

/*
for(std::map<TString, double>::iterator it = xsecs.begin(); it != xsecs.end(); ++it)
	{
	TString dtag = it->first;
	double xsec  = it->second;
	cout << "For dtag " << dtag << " xsec " << xsec << "\n";
	}
*/


std::vector < TString > dtags;
std::vector < TFile * > files;
std::vector < TH1D * > histos;
std::vector < TH1D * > weightflows;
// nick->summed histo
std::map<TString, TH1D *> nicknamed_mc_histos;
//vector<int> dtags;
//dtags.reserve();

// make stack of MC, scaling according to ratio = lumi * xsec / weightflow4 (bin5?)
// also nickname the MC....
// per-dtag for now..

//THStack *hs = new THStack("hs","Stacked 1D histograms");
THStack *hs      = new THStack("hs", "");
TH1D    *hs_data = NULL;

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
	cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	TFile* file = TFile::Open(the_file);
	//files.push_back(file); // not needed
	//dtags.push_back(5);

	if (!file->GetListOfKeys()->Contains(distr))
		{
		if (be_verbose) cout << "no " << distr << endl;
		continue;
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

	if (be_verbose) cout << "got weightflow init " << normal_initial_weight << endl;

	if (be_verbose) cout << "getting " << distr  << " projection to " << projection << endl;
	// get the histogram's projection
	histos.push_back((TH1D*) ((TH3D*) file->Get(distr))->Project3D(projection));
	histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
	if (be_verbose) histos.back()->Print();

	// normalize the histo for bin width
	for (Int_t i=0; i<=histos.back()->GetSize(); i++)
		{
		//yAxis->GetBinLowEdge(3)
		double content = histos.back()->GetBinContent(i);
		double width   = histos.back()->GetXaxis()->GetBinUpEdge(i) - histos.back()->GetXaxis()->GetBinLowEdge(i);
		histos.back()->SetBinContent(i, content/width);
		}

	if (dtag.Contains("Data"))
		{
		if (be_verbose) cout << "summing data-stack" << endl;
		histos.back()->SetMarkerStyle(9);
		//histos.back()->SetFillColor(kRed + i);

		if (hs_data == NULL)
			{
			if (be_verbose) cout << "creating data histo" << endl;
			hs_data = (TH1D*) histos.back()->Clone();
			}
		else
			{
			if (be_verbose) cout << "add histo to data histo" << endl;
			hs_data->Add(histos.back());
			}
		}

	else
		{
		Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
		histos.back()->Scale(ratio);
		if (be_verbose)
			{
			cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
			histos.back()->Print();
			}
		//histos.back()->SetFillColor(kRed );

		std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
		TString nick = nick_colour.first;
		Color_t col = nick_colour.second;

		histos.back()->SetFillColor( col );
		histos.back()->SetMarkerStyle(20);
		histos.back()->SetLineStyle(0);
		histos.back()->SetMarkerColor(i);

		//std::map<TString, TH1D *> nicknamed_mc_histos;
		if (nicknamed_mc_histos.find(nick) == nicknamed_mc_histos.end())
			nicknamed_mc_histos[nick] = (TH1D*) histos.back()->Clone();
			//nicknamed_mc_histos[nick] = new TH1D(nick, "");
		else
			nicknamed_mc_histos[nick]->Add(histos.back());
		//hs->Add(histos.back(), "HIST");
		}

	if (be_verbose) cout << "processed histogram" << endl;
	}


double integral_data = hs_data->Integral();
double integral_MC   = 0.;
//for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
for (std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	//cout << mc_jet_origins[origin_n] << "  integral = " << mc_jet_origin_ths[origin_n]->Integral() << endl;
	//mc_jet_origin_ths[origin_n]->Print();
	TH1D * distr = it->second;
	integral_MC += distr->Integral();
	}
cout << "data  integral = " << integral_data << endl;
cout << "MC    integral = " << integral_MC   << endl;
double ratio = integral_data / integral_MC;
cout << "ratio = " << ratio << endl;

if (normalize_MC_stack)
	{
	cout << "normalizing MC" << endl;
	//for (int origin_n=0; origin_n<mc_jet_origins.size(); origin_n++)
	for (std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
		{
		TString nick = it->first;
		TH1D * distr = it->second;
		double integral = distr->Integral();
		distr->Scale(ratio);
		if (be_verbose) cout << nick << " integral before = " << integral << " after = " << distr->Integral() << endl;
		}
	}

//std::map<TString, TH1D *> nicknamed_mc_histos;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;
	if (be_verbose)
		{
		cout << "adding to mc stack: " << nick << endl;
		distr->Print();
		}

	hs->Add(distr, "HIST");

	//TLegendEntry *entry=leg->AddEntry("NULL","","h");
	//entry=leg->AddEntry("singlemu_altstep4rho3","Single top","F");
	leg->AddEntry(distr, nick, "F");

	//distr->SetName(controlpoint_name.c_str());
	//distr->Write();
	//out_f->Write(controlpoint_name.c_str());
	//cout << "For channel " << channel << " writing " << controlpoint_name << "\n";
	}

//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

if (be_verbose) cout << "setting title" << endl;

//hs->GetXaxis()->SetTitle(distr);
//cst->SetXaxisTile(distr);
//hs_data->GetXaxis()->SetTitle(distr);

/*
TIter next(gDirectory->GetList());
TObject *obj;
while ((obj=next())) {
	if (obj->InheritsFrom("TH1")) {
		TH1 *h = (TH1*)obj;
		h->GetXaxis()->SetTitle(distr);
		}
	}
*/


TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
h->SetXTitle("x");
//cst->Update();
//cst->Modified();

if (be_verbose) cout << "drawing" << endl;

hs_data->Draw("e p");
hs->Draw("same");
hs_data->Draw("e p same"); // to draw it _over_ MC

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");

bool set_logy = false;

if (projection == string("x") || projection == string("z"))
	set_logy = true;

if (set_logy)
	cst->SetLogy();

hs->GetXaxis()->SetTitle(distr);
hs_data->SetXTitle(distr);

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr + "_MCstacked_" + projection + (normalize_MC_stack? "_normalized" : "") + (set_logy? "_logy" : "") + ".png" );

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


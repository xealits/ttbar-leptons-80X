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

#include "TStyle.h"

#include <map>
#include <string>
#include <vector>

#include "dtag_xsecs.h"

#define DTAG_ARGS_START 10
using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 10)
	{
	std::cout << "Usage : " << argv[0] << " normalize logy lumi distr distr_name rebin_factor xmin xmax dir dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString normalize_s(argv[1]);
bool normalize_to_data = false;
if (normalize_s == TString("T") || normalize_s == TString("Y"))
	normalize_to_data = true;

TString logy(argv[2]);
bool set_logy = false;
if (logy == TString("T") || logy == TString("Y"))
	set_logy = true;

double lumi = atof(argv[3]);
TString distr(argv[4]);
TString distr_name(argv[5]);
Int_t rebin_factor(atoi(argv[6]));

double xmin = atof(argv[7]);
double xmax = atof(argv[8]);
bool xlims_set = true;
if (xmin < 0 || xmax < 0)
	xlims_set = false;

TString dir(argv[9]);
TString dtag1(argv[10]);

bool eltau = false, mutau = false;
if (distr.Contains("singleel"))
	eltau = true;
if (distr.Contains("singlemu"))
	mutau = true;

cout << "channel: " << eltau << ' ' << mutau << endl;
cout << "logy: " << set_logy << endl;
cout << lumi  << endl;
cout << distr << endl;
cout << rebin_factor << endl;
cout << xmin << ' ' << xmax << endl;
cout << dir   << endl;
cout << dtag1 << endl;


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
TLegend* leg = new TLegend(0.7, 0.7, 0.89, 0.89);

for (int i = DTAG_ARGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	files.push_back(TFile::Open(dir + "/" + dtag + ".root"));
	//dtags.push_back(5);

	if (!files.back()->GetListOfKeys()->Contains(distr))
		{
		cout << "no " << distr << endl;
		continue;
		}

	weightflows.push_back((TH1D*) files.back()->Get("weightflow_el_NOMINAL"));
	weightflows.back()->Print();

	histos.push_back((TH1D*) files.back()->Get(distr));
	histos.back()->Rebin(rebin_factor); // should rebin the new histo inplace
	histos.back()->Print();
	if (dtag.Contains("Data"))
		{
		cout << "summing data-stack" << endl;
		histos.back()->SetMarkerStyle(9);
		//histos.back()->SetFillColor(kRed + i);

		if (hs_data == NULL)
			{
			cout << "creating data histo" << endl;
			hs_data = (TH1D*) histos.back()->Clone();
			}
		else
			{
			cout << "add histo to data histo" << endl;
			hs_data->Add(histos.back());
			}
		}
	else
		{
		Double_t ratio = lumi * xsecs[dtag] / weightflows.back()->GetBinContent(11);
		histos.back()->Scale(ratio);
		cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		histos.back()->Print();
		//histos.back()->SetFillColor(kRed );

		std::pair<TString, Color_t> nick_colour = dtag_nick_colour(dtag);
		TString nick = nick_colour.first;
		Color_t col = nick_colour.second;

		//cout << dtag.Contains("aeltu") << endl;

		/*
		if (dtag.Contains("aeltu") && eltau)
			{
			cout << "setting signal" << endl;
			col = kOrange;
			nick = TString("tt_eltau");
			}
		if (dtag.Contains("amtuu") && mutau)
			{
			cout << "setting signal" << endl;
			col = kGreen - 3;
			nick = TString("tt_mutau");
			}
		*/

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
	}

// build the sum of MC
TH1D    *hs_sum  = NULL;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;
	if (hs_sum == NULL)
		{
		hs_sum = (TH1D*) (distr->Clone("brand_new_hs_sum"));
		hs_sum->SetName("hs_sum");
		hs_sum->ResetBit(kCanDelete);
		}
	else
		hs_sum->Add(distr);
	}

// scale MC to Data
double scale = 1.0;
// normalize MC integral if requested
if (normalize_to_data)
	{
	scale = hs_sum->Integral() / hs_data->Integral();
	}

// build the stack of MC histos
//std::map<TString, TH1D *> nicknamed_mc_histos;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;

	// normalize MC sum to Data
	if (normalize_to_data)
		distr->Scale(scale);

	cout << "adding to mc stack: " << nick << endl;
	distr->Print();

	hs->Add(distr, "HIST");

	//TLegendEntry *entry=leg->AddEntry("NULL","","h");
	//entry=leg->AddEntry("singlemu_altstep4rho3","Single top","F");
	leg->AddEntry(distr, nick, "F");
	}

//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

TStyle* m_gStyle = new  TStyle();
m_gStyle->SetOptStat(111100);

cout << "setting title" << endl;

//hs->GetXaxis()->SetTitle(distr);
//cst->SetXaxisTile(distr);
//hs_data->GetXaxis()->SetTitle(distr);


//TH1F * h = cst->DrawFrame(0.,0.,1.,1.);
//h->SetXTitle("x");
//cst->Update();
//cst->Modified();

if (xlims_set)
	{
	cout << "setting X limits " << xmin << ' ' << xmax << endl;
	hs_data->GetXaxis()->SetRange(xmin, xmax);
	hs_data->GetXaxis()->SetRangeUser(xmin, xmax);
	//hs->GetXaxis()->SetRange(xmin, xmax);
	//hs->GetXaxis()->SetRangeUser(xmin, xmax);
	}

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

//cst->Divide(1, 2);

/*
TPad *pad1 = (TPad *)(cst->cd(1)); 
TPad *pad2 = (TPad *)(cst->cd(2)); 
*/

//cst->cd(1);
pad1->cd();

hs_data->GetXaxis()->SetLabelFont(63);
hs_data->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
hs_data->GetYaxis()->SetLabelFont(63);
hs_data->GetYaxis()->SetLabelSize(14); // labels will be 14 pixels

cout << "drawing" << endl;

hs_data->Draw("p"); // drawing data-MC-data to have them both in the range of the plot
hs->Draw("same"); // the stack
// also draw the sum of the stack and its' errors:
//((TObjArray*) hs->GetStack())->Last()->Draw("e4 same"); // then, how to set styles with the hatched error bars?

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

hs_data->Draw("e p same"); // the data with the errors

// histo->SetTitle("boxtitle;x axis title [unit];y axis title [unit]")
cout << "setting the stack title" << endl;
hs->GetXaxis()->SetTitle(distr_name);
cout << "done setting the stack title" << endl;
hs_data->SetXTitle(distr_name);
hs_sum->SetXTitle(distr_name);

leg->SetBorderSize(0);
leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");



//cst->SetLogy();
//cst->Modified();

pad2->cd();
//cst->cd(2);

TH1D * hs_data_relative = (TH1D*) hs_data->Clone();
TH1D * hs_sum_relative  = (TH1D*) hs_sum->Clone();

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

if (xlims_set)
	{
	hs_sum_relative->GetXaxis()->SetRange(xmin, xmax);
	hs_sum_relative->GetXaxis()->SetRangeUser(xmin, xmax);
	//hs_data_relative->GetXaxis()->SetRange(xmin, xmax);
	//hs_data_relative->GetXaxis()->SetRangeUser(xmin, xmax);
	}

hs_sum_relative->Draw("e2");
hs_data_relative->Draw("e p same");

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr + (set_logy? "_logy" : "") + (normalize_to_data? "_normalizedToData.png" : ".png") );

return 0;
}


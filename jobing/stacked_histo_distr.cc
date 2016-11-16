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

#define DTAG_ARGS_START 5
using namespace std;

// nick and colour
std::pair<TString, Color_t> dtag_nick_colour(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("dyjets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("wjets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("singletop", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qqbar")) return std::make_pair("tt_jj", kGreen+4);
		else if (dtag.Contains("elqbar") || dtag.Contains("qelbar") ||dtag.Contains("muqbar") || dtag.Contains("qmubar") || dtag.Contains("tauqbar") || dtag.Contains("qtaubar")) return std::make_pair("tt_lj", kGreen+3);
		else if (dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt_em", kGreen-9);
		else if (dtag.Contains("elelbar")) return std::make_pair("tt_ee", kAzure-9);
		else if (dtag.Contains("mumubar")) return std::make_pair("tt_mm", kYellow-7);
		else if (dtag.Contains("eltaubar") || dtag.Contains("tauelbar")) return std::make_pair("tt_et", kOrange+4);
		else if (dtag.Contains("mutaubar") || dtag.Contains("taumubar")) return std::make_pair("tt_mt", kOrange+1);
		else return std::make_pair("tt_other", kYellow+1);
		}
	else return std::make_pair("other", kBlack);

	}

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 5)
	{
	std::cout << "Usage : " << argv[0] << " lumi distr rebin_factor dir dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

double lumi = atof(argv[1]);
TString distr(argv[2]);
Int_t rebin_factor(atoi(argv[3]));
TString dir(argv[4]);
TString dtag1(argv[5]);

cout << lumi  << endl;
cout << distr << endl;
cout << dir   << endl;
cout << dtag1 << endl;

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

	weightflows.push_back((TH1D*) files.back()->Get("weightflow_el"));
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
		Double_t ratio = lumi * xsecs[dtag] / weightflows.back()->GetBinContent(5);
		histos.back()->Scale(ratio);
		cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		histos.back()->Print();
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
	}

//std::map<TString, TH1D *> nicknamed_mc_histos;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;
	cout << "adding to mc stack: " << nick << endl;
	distr->Print();

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

cout << "setting title" << endl;

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

cout << "drawing" << endl;

hs->Draw();
hs_data->Draw("e p same");

leg->Draw();

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");


hs->GetXaxis()->SetTitle(distr);
hs_data->SetXTitle(distr);

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr + ".png" );

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


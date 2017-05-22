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

#include "TNtuple.h"

#define INPUT_DTAGS_START 4
#define NTUPLE_NAME "ntuple"

using namespace std;

TH1D* pull_likelihood_regions(TNtuple* ntuple, TString& histo_name)
	{
	TH1D* h = new TH1D(histo_name, ";;", 15, 0, 15)
	int N;
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets == 1) && (njets == 3) ")
	h->SetBinContent(0,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets == 1) && (njets == 4) ")
	h->SetBinContent(1,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets == 1) && (njets == 5) ")
	h->SetBinContent(2,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets == 1) && (njets > 5) ")
	h->SetBinContent(3,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets > 1) && (njets == 3) ")
	h->SetBinContent(4,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets > 1) && (njets == 4) ")
	h->SetBinContent(5,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets > 1) && (njets == 5) ")
	h->SetBinContent(6,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (nbjets > 1) && (njets > 5) ")
	h->SetBinContent(7,N)                                                                                                                                                           
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets == 1) && (njets == 3) ")
	h->SetBinContent(8,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets == 1) && (njets == 4) ")
	h->SetBinContent(9,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets == 1) && (njets == 5) ")
	h->SetBinContent(10,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets == 1) && (njets > 5) ")
	h->SetBinContent(11,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets > 1) && (njets == 4) ")
	h->SetBinContent(12,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets > 1) && (njets == 5) ")
	h->SetBinContent(13,N)
	N = ntuple->Draw("njets", "(met_corrected > 40) && (fabs(leps_ID) == 13) && (lj_peak_distance < 500) && (nbjets > 1) && (njets > 5) ")
	h->SetBinContent(14,N)
	//h->Draw()

	return h;
	}

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
char usage_string[128] = "be_verbose lumi dir dtags";
if (argc < 4)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

TString verb(argv[1]);
bool be_verbose = false;
if (verb==TString("T") || verb==TString("Y"))
	be_verbose = true;
double lumi = atof(argv[2]);
TString dir(argv[3]);
TString dtag1(argv[4]);

cout << be_verbose  << endl;
cout << lumi  << endl;
cout << dir   << endl;
cout << dtag1 << endl;

TH1D* data_histo = NULL;
std::map<TString, TH1D*> nicknamed_mc_histos;


/*
 * Get each dtag file in the reduced dir,
 * open the ntuple, puul the N-events per event region distribution
 * if the dtag is Data -> sum the projection to data-histo
 * if it is MC -> scale it according xsecs/weightflow (and stack it to mc-stack)

 * save them all to the jobsums/likelihood_regions...root file
 */
for (int i = INPUT_DTAGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	if (be_verbose) cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	TFile* file = TFile::Open(the_file);

	if (dtag.Contains("Data"))
		{
		if (be_verbose) cout << "summing data-stack" << endl;
		//pull_likelihood_regions if the ntuple present

		if (!file->GetListOfKeys()->Contains(NTUPLE_NAME))
			{
			if (be_verbose) cout << "no " << NTUPLE_NAME << endl;
			continue;
			}

		TNtuple* ntuple = (TNtuple*) file->Get(NTUPLE_NAME);

		TH1D* histo = pull_likelihood_regions(ntuple, dtag + "_regions");

		histo->SetMarkerStyle(9);

		if (data_histo == NULL)
			{
			if (be_verbose) cout << "creating data histo" << endl;
			data_histo = (TH1D*) histo->Clone();
			data_histo->SetName("data_event_regions");
			}
		else
			{
			if (be_verbose) cout << "add histo to data histo" << endl;
			data_histo->Add(histo);
			}
		}
	else
		{
		/* in MC do the same
		 * but scale the histo by weightflow
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

		if (be_verbose) cout << "got weightflow init" << endl;

		// MC ratio for this dtag:
		Double_t ratio = lumi * xsecs[dtag] / normal_initial_weight;
		if (be_verbose) cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;


		if (!file->GetListOfKeys()->Contains(NTUPLE_NAME))
			{
			if (be_verbose) cout << "no " << NTUPLE_NAME << endl;
			continue;
			}

		TNtuple* ntuple = (TNtuple*) file->Get(NTUPLE_NAME);

		TH1D* histo = pull_likelihood_regions(ntuple, dtag + "_regions");

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

		// save the histogram according to its' nick
		if (nicknamed_mc_histos.find(nick) == nicknamed_mc_histos.end())
			{
			nicknamed_mc_histos[nick] = (TH1D*) histo->Clone(); // TODO: should I really change the name here?
			nicknamed_mc_histos[nick].SetName(nick + "_regions");
			}
		else
			nicknamed_mc_histos[nick]->Add(histo);
		}
	}

// let's also build the stack of MC histos
THStack *hs      = new THStack("hs", "");

for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;

	if (be_verbose) cout << nick;

	hs->Add(distr, "HIST");
	}


TFile* out_f = TFile::Open (dir + TString("/jobsums/Likelihood_Regions.root"), "CREATE");
data_histo->Write();

for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos->begin(); it != nicknamed_mc_histos->end(); ++it)
	{
	TString nickname = it->first;
	TH1D * distr = it->second;
	//distr->SetName(controlpoint_name.c_str());
	distr->Write();
	//out_f->Write(controlpoint_name.c_str()); // TODO: check if it's needed here
	}

hs->Write();

out_f->Write();
out_f->Close();
}


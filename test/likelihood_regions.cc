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

TH1D* pull_likelihood_regions(TTree& NT_output_ttree, TString& histo_name, TString& dtag)
	{
	// the interface to NT_output_ttree
	#define NTUPLE_INTERFACE_OPEN
	//#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_v11_3.h"
	#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput.h"

	bool is_MC = dtag.Contains("MC");
	bool is_W0Jets = dtag.Contains("W0Jet"); // needed for handling WNJets
	bool is_aMCatNLO = dtag.Contains("amcatnlo");

	// the histogram with yields in event categories
	TH1D* h = new TH1D(histo_name, ";;", 15, 0, 15);

	/* Loop through events, find weights, skip stuff, record if they pass to any category
	 */
	for (Long64_t i=0; i<NT_output_ttree.GetEntries(); i++)
	//for (Long64_t i=0; i<10; i++)
		{
		NT_output_ttree.GetEntry(i);
		// test:
		//cout << NT_NUP_gen << '\t' << NT_aMCatNLO_weight << '\t' << NT_nvtx_gen << '\t' << (int) NT_nvtx_gen << endl;

		int event_category = -1;
		double weight = 1; // for MC

		// general requirement for all events:
		if (NT_met_corrected.pt() < 40 || fabs(NT_leps_ID) != 13) continue;
		// the event with taus (2 = medium):
		if (NT_tau_IDlev_0 < 2. && NT_tau_IDlev_1 < 2.) continue;

		// find basic weight:
		//  * -1 aMCatNLO & NUP for WJets
		//  * pileup

		// NUP == 5 for W0Jets
		if (is_W0Jets && NT_NUP_gen != 5) continue;

		// amcatnlo MC has these weights
		if (is_aMCatNLO)
			weight *= (NT_aMCatNLO_weight > 0? 1 : -1);

		// and pile-up:
		if (is_MC)
			weight *= pileup_ratio[(int)NT_nvtx_gen];

		// select event categories
		const double lj_dist = 500;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets == 3)) event_category = 0;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets == 4)) event_category = 1;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets == 5)) event_category = 2;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets == 1) && (NT_njets > 5) ) event_category = 3;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets == 3)) event_category = 4;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets == 4)) event_category = 5;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets == 5)) event_category = 6;
		if ((NT_lj_peak_distance > lj_dist) && (NT_nbjets > 1)  && (NT_njets > 5) ) event_category = 7;

		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets == 3)) event_category = 8;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets == 4)) event_category = 9;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets == 5)) event_category =10;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets == 1) && (NT_njets > 5) ) event_category =11;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets == 4)) event_category =12;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets == 5)) event_category =13;
		if ((NT_lj_peak_distance < lj_dist) && (NT_nbjets > 1)  && (NT_njets > 5) ) event_category =14;

		h->Fill(event_category, weight);
		}

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
	TString histo_name = dtag + TString("_regions");
	if (be_verbose) cout << "processing " << dtag << endl;
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	TFile* file = TFile::Open(the_file);

	if (!file->GetListOfKeys()->Contains(NTUPLE_NAME))
		{
		if (be_verbose) cout << "no " << NTUPLE_NAME << endl;
		continue;
		}

	// it's actually TNtuple, which descends from TTree thus should open ok
	TTree* ntuple = (TTree*) file->Get(NTUPLE_NAME);

	TH1D* histo = pull_likelihood_regions(*ntuple, histo_name, dtag);


	if (dtag.Contains("Data"))
		{
		if (be_verbose) cout << "summing data-stack" << endl;
		//pull_likelihood_regions if the ntuple present

		if (!file->GetListOfKeys()->Contains(NTUPLE_NAME))
			{
			if (be_verbose) cout << "no " << NTUPLE_NAME << endl;
			continue;
			}

		/*
		// it's actually TNtuple, which descends from TTree thus should open ok
		TTree* ntuple = (TTree*) file->Get(NTUPLE_NAME);
		TH1D* histo = pull_likelihood_regions(*ntuple, histo_name);
		*/

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

		// removing QCD from mu-tau selection (suddenly couple events show up in 2 jets or whatever)
		if (dtag.Contains("QCD")) continue;

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
			TString mc_histoname = nick + TString("_regions");
			nicknamed_mc_histos[nick]->SetName(mc_histoname);
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

if (be_verbose) cout << "built MC stack" << endl;

// and write everything out:
// data histo
// separate MC histos (separated by their nicknames)
// and the MC stack
TFile* out_f = TFile::Open (dir + TString("/jobsums/Likelihood_Regions_2_with_taus_without_qcd.root"), "CREATE");
if (be_verbose) cout << "opened output file" << endl;
data_histo->Write();

for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nickname = it->first;
	TH1D * distr = it->second;
	//distr->SetName(controlpoint_name.c_str());
	distr->Write();
	//out_f->Write(controlpoint_name.c_str()); // TODO: check if it's needed here
	}

hs->Write();

if (be_verbose) cout << "wrote objects" << endl;

out_f->Write();
out_f->Close();
}


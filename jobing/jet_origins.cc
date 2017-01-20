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

using namespace std;

/*
int stacked_histo_distr (int argc, char *argv[])

where really arguments are:
(TString distr, TString dir, TString dtags)

convert char* dtags[] to TStrings,
have dtag->xsec map around,
open dir/dtag.root files,
get the distr and weightflows (for ratios, so doesn't matter which) from there,
make stack of ratio-weighted MC-s (so, separate MC and Data)
and overlay it with data histo
+ error bars everywhere

*/

/*
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
*/

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 4)
	{
	std::cout << "Usage : " << argv[0] << " job_dir dtag origin_channel [origins_channels]" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString job_dir(argv[1]);
TString dtag(argv[2]);

//std::vector <TString> channels = {"HLTjet_qcd_jet_origins", "HLTmu_qcd_jet_origins", "HLTjetmu_qcd_jet_origins",
//	"HLTjet_wjets_jet_origins", "HLTmu_wjets_jet_origins", "HLTjetmu_wjets_jet_origins" };

//std::vector <TString> channels = {"HLTjet_qcd_jet_partonFlavour", "HLTmu_qcd_jet_partonFlavour", "HLTjetmu_qcd_jet_partonFlavour",
//	"HLTjet_wjets_jet_partonFlavour", "HLTmu_wjets_jet_partonFlavour", "HLTjetmu_wjets_jet_partonFlavour" };
// in ttbar dieptons:
//std::vector <TString> channels = {"dilep_passjets_jet_partonFlavour", "dilep_passjetsNbtag_jet_partonFlavour", "elel_passjets_jet_partonFlavour"};

std::vector<TString> channels;
for (int i=3; i<argc; i++)
	{
	TString a_channel(argv[i]);
	channels.push_back(a_channel);
	}

//cout << job_dir  << endl;
//cout << dtag << endl;

// Find the MC ratio
double xsec = xsecs[dtag];
Double_t ratio = 0; // xsec / weightflow->GetBinContent(4);

TFile * file = TFile::Open(job_dir + "/" + dtag + ".root");
//TH1D * weightflow = (TH1D*) file->Get("weightflow");
TH1D * weightflow;

// for ttbar dileptons
if (file->GetListOfKeys()->Contains("weightflow_elel"))
	{
	weightflow = (TH1D*) file->Get("weightflow_elel");
	ratio = xsec / weightflow->GetBinContent(5);
	}
else if (file->GetListOfKeys()->Contains("weightflow"))
	{
	weightflow = (TH1D*) file->Get("weightflow");
	ratio = xsec / weightflow->GetBinContent(4);
	}
else
	{
	cerr << "no weightflow distro" << endl;
	return 1;
	}


//cout << "dtag " << dtag << " ratio = " << ratio << endl;


cout << "dtag,ratio";
for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch)
	{
	/*
	if (!file->GetListOfKeys()->Contains(*ch))
		{
		cerr << "no " << *ch << endl;
		continue;
		}
	*/
	cout << "," << *ch << "_o," << *ch << "_t," << *ch << "_b," << *ch << "_g," << *ch << "_a";
	}
cout << endl;

cout << dtag << "," << ratio;

for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch)
	{
	if (!file->GetListOfKeys()->Contains(*ch))
		{
		//cerr << "no " << *ch << endl;
		cout << ",NA,NA,NA,NA,NA";
		continue;
		}

	TH1D * jet_origins = (TH1D*) file->Get(*ch);

	//Double_t gluons_to_all = jet_origins->GetBinContent(22) / jet_origins->Integral();
	cout << "," << jet_origins->GetBinContent(1) << "," << jet_origins->GetBinContent(16) << "," << jet_origins->GetBinContent(6) << "," << jet_origins->GetBinContent(22) << "," << jet_origins->Integral();
	}
cout << endl;

/*
TH1D * jet_origins_HLTjet_QCD = (TH1D*) file->Get("HLTjet_qcd_jet_origins");
// other origins:

cout << "HLTjet_qcd_jet_origins\n";
cout << "gluons\tall\n";
cout << jet_origins_HLTjet_QCD->GetBinContent(22)*ratio << "\t" << jet_origins_HLTjet_QCD->Integral()*ratio << endl;
*/

return 0;
}


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
#include "TPaveText.h"
#include "TMath.h"

#include "TStyle.h"

#include <map>
#include <string>
#include <vector>

#include "dtag_xsecs.h"

#define DTAG_ARGS_START 13
using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < 13)
	{
	std::cout << "Usage : " << argv[0] << " normalize logy lumi distr mc_distr distr_name yname rebin_factor xmin xmax ratio_range dir dtags" << std::endl;
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
TString distr_data(argv[4]);
TString distr_mc(argv[5]);
if (distr_mc == TString("f") || distr_mc == TString("same"))
	distr_mc = distr_data;
TString distr_name(argv[6]);
TString yname(argv[7]);

Int_t rebin_factor(atoi(argv[8]));

double xmin = atof(argv[9]);
double xmax = atof(argv[10]);
double ratio_range = atof(argv[11]);
bool xlims_set = true;
if (xmin < 0 || xmax < 0)
	xlims_set = false;

TString dir(argv[12]);
TString dtag1(argv[13]);

bool eltau = false, mutau = false;
if (distr_data.Contains("singleel"))
	eltau = true;
if (distr_data.Contains("singlemu"))
	mutau = true;

cout << "channel: " << eltau << ' ' << mutau << endl;
cout << "logy: " << set_logy << endl;
cout << lumi  << endl;
cout << distr_data << endl;
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

vector<TH1D*> hs_mc_sys = {NULL, NULL, NULL, NULL, NULL, NULL};

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
	bool isData = dtag.Contains("Data");
	//TString distr = (isData? distr_data : distr_mc);
	TString distr = distr_data;

	if (!files.back()->GetListOfKeys()->Contains(distr))
		{
		cout << "no " << distr << endl;
		continue;
		}

	weightflows.push_back((TH1D*) files.back()->Get("weightflow_el_NOMINAL"));
	weightflows.back()->Print();

	TH1D* histo = (TH1D*) files.back()->Get(distr);

	Double_t ratio = 1;
	if (!isData)
		{
		ratio = lumi * xsecs[dtag] / weightflows.back()->GetBinContent(11);
		histo->Scale(ratio);
		cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		histo->Print();
		//histos.back()->SetFillColor(kRed );
		}
	if (rebin_factor != 1) histo->Rebin(rebin_factor); // should rebin the new histo inplace
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

	histos.push_back(histo);

	histos.back()->Print();
	if (isData)
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

		// for MC also add up the systematic histograms
		// elmu_passjets_mediumTaus_jets_distr_pt  elmu_passjets_PU_UP_mediumTaus_jets_distr_pt
		TString distr_PU_UP = distr;
		TString distr_PU_DOWN = distr;
		TString distr_JER_UP = distr;
		TString distr_JER_DOWN = distr;
		TString distr_JES_UP = distr;
		TString distr_JES_DOWN = distr;
		distr_PU_UP    .ReplaceAll("_vtightTaus", "_PU_UP_vtightTaus"); distr_PU_DOWN  .ReplaceAll("_vtightTaus", "_PU_DOWN_vtightTaus"); distr_JER_UP   .ReplaceAll("_vtightTaus", "_JER_UP_vtightTaus"); distr_JER_DOWN .ReplaceAll("_vtightTaus", "_JER_DOWN_vtightTaus"); distr_JES_UP   .ReplaceAll("_vtightTaus", "_JES_UP_vtightTaus"); distr_JES_DOWN .ReplaceAll("_vtightTaus", "_JES_DOWN_vtightTaus");
		distr_PU_UP    .ReplaceAll("_tightTaus", "_PU_UP_tightTaus"); distr_PU_DOWN  .ReplaceAll("_tightTaus", "_PU_DOWN_tightTaus"); distr_JER_UP   .ReplaceAll("_tightTaus", "_JER_UP_tightTaus"); distr_JER_DOWN .ReplaceAll("_tightTaus", "_JER_DOWN_tightTaus"); distr_JES_UP   .ReplaceAll("_tightTaus", "_JES_UP_tightTaus"); distr_JES_DOWN .ReplaceAll("_tightTaus", "_JES_DOWN_tightTaus");
		distr_PU_UP    .ReplaceAll("_mediumTaus", "_PU_UP_mediumTaus"); distr_PU_DOWN  .ReplaceAll("_mediumTaus", "_PU_DOWN_mediumTaus"); distr_JER_UP   .ReplaceAll("_mediumTaus", "_JER_UP_mediumTaus"); distr_JER_DOWN .ReplaceAll("_mediumTaus", "_JER_DOWN_mediumTaus"); distr_JES_UP   .ReplaceAll("_mediumTaus", "_JES_UP_mediumTaus"); distr_JES_DOWN .ReplaceAll("_mediumTaus", "_JES_DOWN_mediumTaus");
		distr_PU_UP    .ReplaceAll("_looseTaus", "_PU_UP_looseTaus"); distr_PU_DOWN  .ReplaceAll("_looseTaus", "_PU_DOWN_looseTaus"); distr_JER_UP   .ReplaceAll("_looseTaus", "_JER_UP_looseTaus"); distr_JER_DOWN .ReplaceAll("_looseTaus", "_JER_DOWN_looseTaus"); distr_JES_UP   .ReplaceAll("_looseTaus", "_JES_UP_looseTaus"); distr_JES_DOWN .ReplaceAll("_looseTaus", "_JES_DOWN_looseTaus");
		distr_PU_UP    .ReplaceAll("_vlooseTaus", "_PU_UP_vlooseTaus"); distr_PU_DOWN  .ReplaceAll("_vlooseTaus", "_PU_DOWN_vlooseTaus"); distr_JER_UP   .ReplaceAll("_vlooseTaus", "_JER_UP_vlooseTaus"); distr_JER_DOWN .ReplaceAll("_vlooseTaus", "_JER_DOWN_vlooseTaus"); distr_JES_UP   .ReplaceAll("_vlooseTaus", "_JES_UP_vlooseTaus"); distr_JES_DOWN .ReplaceAll("_vlooseTaus", "_JES_DOWN_vlooseTaus");
		vector<TString*> sys_distrs = { &distr_PU_UP, &distr_PU_DOWN, &distr_JER_UP, &distr_JER_DOWN, &distr_JES_UP, &distr_JES_DOWN };

		for (int i=0; i< sys_distrs.size(); i++)
			{
			TString& distr = *sys_distrs[i];
			cout << "sys:\t" << distr << '\t';
		        if (!files.back()->GetListOfKeys()->Contains(distr))
		                {
		                cout << "no " << distr << endl;
		                continue;
				}

			TH1D* histo = (TH1D*) files.back()->Get(distr);

			histo->Scale(ratio);
			if (rebin_factor != 1) histo->Rebin(rebin_factor); // should rebin the new histo inplace
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

			if (hs_mc_sys[i])
				{
				cout << "adding" << '\t';
				hs_mc_sys[i]->Add(histo);
				cout << "added" << endl;
				}
			else
				{
				cout << "clonning" << '\t';
				hs_mc_sys[i] = (TH1D*) histo->Clone();
				cout << "cloned" << endl;
				}
			}

		}
	}
// build the sum of MC
// it has to be done before building the stack because ROOT is crap
// yep, you cannot scale a stack after it's been built
TH1D    *hs_sum  = NULL;
for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_histos.begin(); it != nicknamed_mc_histos.end(); ++it)
	{
	TString nick = it->first;
	TH1D * distr = it->second;
	cout << nick;

	// build the summ of MC
	cout << " summing mc" << endl;
	if (hs_sum == NULL)
		{
		hs_sum = (TH1D*) (distr->Clone("brand_new_hs_sum"));
		hs_sum->SetName("hs_sum");
		hs_sum->ResetBit(kCanDelete);
		}
	else
		hs_sum->Add(distr);
	}

// now hs_sum VS systematics can be done
TH1D* hs_mc_PU_UP    = hs_mc_sys[0];
TH1D* hs_mc_PU_DOWN  = hs_mc_sys[1];
TH1D* hs_mc_JER_UP   = hs_mc_sys[2];
TH1D* hs_mc_JER_DOWN = hs_mc_sys[3];
TH1D* hs_mc_JES_UP   = hs_mc_sys[4];
TH1D* hs_mc_JES_DOWN = hs_mc_sys[5];
cout << "sys subtraction " << endl;
hs_sum->Print();
hs_mc_PU_UP   ->Print();
hs_mc_PU_DOWN ->Print();
hs_mc_JES_UP   ->Print();
hs_mc_JES_DOWN ->Print();
hs_mc_JER_UP   ->Print();
hs_mc_JER_DOWN ->Print();
//hs_mc_PU_UP   ->Add(hs_sum, -1);
//hs_mc_PU_DOWN ->Add(hs_sum, -1);
//hs_mc_JER_UP  ->Add(hs_sum, -1);
//hs_mc_JER_DOWN->Add(hs_sum, -1);
//hs_mc_JES_UP  ->Add(hs_sum, -1);
//hs_mc_JES_DOWN->Add(hs_sum, -1);
TFile* test_f = TFile::Open( dir + "/jobsums/TEST_" + distr_mc + (set_logy? "_logy" : "") + (normalize_to_data? "_normalizedToData.root" : ".root"), "RECREATE" );

hs_mc_PU_UP   ->Write();
hs_mc_PU_DOWN ->Write();
hs_mc_JER_UP  ->Write();
hs_mc_JER_DOWN->Write();
hs_mc_JES_UP  ->Write();
hs_mc_JES_DOWN->Write();

test_f->Write();
test_f->Close();

// now set this errors in mc sum
for (int i=0; i<=hs_sum->GetSize(); i++)
	{
	double mc_content = hs_sum->GetBinContent(i);
	double mc_error   = hs_sum->GetBinError(i);

	double mc_PU_UP  = hs_mc_PU_UP ->GetBinContent(i);
	double mc_JER_UP = hs_mc_JER_UP->GetBinContent(i);
	double mc_JES_UP = hs_mc_JES_UP->GetBinContent(i);

	double mc_PU_DOWN  = hs_mc_PU_DOWN ->GetBinContent(i);
	double mc_JER_DOWN = hs_mc_JER_DOWN->GetBinContent(i);
	double mc_JES_DOWN = hs_mc_JES_DOWN->GetBinContent(i);

	double sys_PU   = (abs(hs_mc_PU_UP  ->GetBinContent(i) - mc_content) + abs(hs_mc_PU_DOWN->GetBinContent(i) - mc_content ))/2;
	double sys_JER  = (abs(hs_mc_JER_UP ->GetBinContent(i) - mc_content) + abs(hs_mc_JER_DOWN->GetBinContent(i) - mc_content))/2;
	double sys_JES  = (abs(hs_mc_JES_UP ->GetBinContent(i) - mc_content) + abs(hs_mc_JES_DOWN->GetBinContent(i) - mc_content))/2;

	double sys_error = TMath::Sqrt(sys_PU*sys_PU + sys_JER*sys_JER + sys_JES*sys_JES + mc_error*mc_error);
	//cout << '(' << mc_PU_UP << '_' << mc_JER_UP << '_' << mc_JES_UP << '-' << mc_content << ',';
	//cout << mc_PU_DOWN << '_' << mc_JER_DOWN << '_' << mc_JES_DOWN << '-' << mc_content << '=';
	//cout << mc_error << '_' << sys_error << ") ";

	hs_sum->SetBinError(i, sys_error);	
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
	cout << nick;

	// each histogram has to be normalized before going into Stack
	// because ROOT is crap
	if (normalize_to_data)
		{
		distr->Scale(scale);
		}

	cout << " adding to mc stack" << endl;
	distr->Print();

	hs->Add(distr, "HIST");

	//TLegendEntry *entry=leg->AddEntry("NULL","","h");
	//entry=leg->AddEntry("singlemu_altstep4rho3","Single top","F");
	leg->AddEntry(distr, nick, "F");
	}

/* scale the stack and the sum
*/
// normalize MC sum to Data
if (normalize_to_data)
	{
	cout << "normalizeing MC" << endl;
	hs_sum->Scale(scale);
	//hs->Scale(scale); // ROOT is crap, THStack doesn't have Scale method
	}

cout << "Data sum = " << hs_data->Integral() << endl;
cout << "MC   sum = " << hs_sum->Integral()  << endl;

//cst->SetFillColor(41);
//cst->Divide(1,1);
// in top left pad, draw the stack with defaults
//cst->cd(1);

TStyle* m_gStyle = new  TStyle();
m_gStyle->SetOptStat(111100);

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
cout << "setting the stack titles " << distr_name << " " << yname << endl;
hs->GetXaxis()->SetTitle(distr_name);
hs->GetYaxis()->SetTitle(yname);
cout << "done setting the stack title" << endl;
hs_data->SetXTitle(distr_name);
hs_sum->SetXTitle(distr_name);
hs_data->SetYTitle(yname);
hs_sum-> SetYTitle(yname);

hs->     GetYaxis()->SetTitleOffset(1.4);
hs_data->GetYaxis()->SetTitleOffset(1.4);
hs_sum-> GetYaxis()->SetTitleOffset(1.4);

cout << "setting title" << endl;
// title for the plot
TPaveText* left_title = new TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC");
left_title->AddText("CMS preliminary at 13 TeV");
left_title->SetTextFont(1);
left_title->SetFillColor(0);
cout << "drawing left title" << endl;
left_title->Draw("same");

TPaveText* right_title = new TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC");
TString s_title(""); s_title.Form("L = %.1f fb^{-1}", lumi/1000);
right_title->AddText(s_title);
right_title->SetTextFont(132);
right_title->SetFillColor(0);
cout << "drawing right title" << endl;
right_title->Draw("same");

leg->SetBorderSize(0);
leg->Draw("same");

//MC_stack_124->GetXaxis()->SetTitle("#rho");
//mcrelunc929->GetYaxis()->SetTitle("Data/#Sigma MC");



//cst->SetLogy();
//cst->Modified();

// THE RATIO PLOT
pad2->cd();
//cst->cd(2);

TH1D * hs_data_relative = (TH1D*) hs_data->Clone();
TH1D * hs_sum_relative  = (TH1D*) hs_sum->Clone();

hs_data_relative->SetXTitle("");
hs_sum_relative-> SetXTitle("");
hs_data_relative->SetYTitle("");
hs_sum_relative-> SetYTitle("");

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
hs_sum_relative->GetYaxis()->SetRange(1 - ratio_range, 1 + ratio_range);
hs_sum_relative->GetYaxis()->SetRangeUser(1 - ratio_range, 1 + ratio_range);
hs_data_relative->GetYaxis()->SetRange(1 - ratio_range, 1 + ratio_range);
hs_data_relative->GetYaxis()->SetRangeUser(1 - ratio_range, 1 + ratio_range);

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

cst->SaveAs( dir + "/jobsums/" + distr_mc + (set_logy? "_logy" : "") + (normalize_to_data? "_normalizedToData.png" : ".png") );

return 0;
}


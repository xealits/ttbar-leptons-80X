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

#include "TStyle.h"

#include <map>
#include <string>
#include <vector>

#include "dtag_xsecs.h"

#define DTAG_ARGS_START 18
using namespace std;

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < DTAG_ARGS_START)
	{
	std::cout << "Usage : " << argv[0] << "scale_width inclusive_xsec normalize logy lumi distr mc_distr name_suffix distr_name yname rebin_factor xmin xmax ymin ymax ratio_range dir dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString scaling_width(argv[1]);
bool scale_width = false;
if (scaling_width == TString("T") || scaling_width == TString("Y"))
	scale_width = true;

TString inclusive(argv[2]); // for ntuple output
bool inclusive_xsec = false;
if (inclusive == TString("T") || inclusive == TString("Y"))
	inclusive_xsec = true;

TString normalize_s(argv[3]);
bool normalize_to_data = false;
if (normalize_s == TString("T") || normalize_s == TString("Y"))
	normalize_to_data = true;

TString logy(argv[4]);
bool set_logy = false;
if (logy == TString("T") || logy == TString("Y"))
	set_logy = true;

double lumi = atof(argv[5]);
TString distr_data(argv[6]);
TString distr_mc(argv[7]);
if (distr_mc == TString("f") || distr_mc == TString("same"))
	distr_mc = distr_data;
TString suffix(argv[8]);
TString distr_name(argv[9]);
TString yname(argv[10]);

Int_t rebin_factor(atoi(argv[11]));

double xmin = atof(argv[12]);
double xmax = atof(argv[13]);
double y_min = atof(argv[14]);
double y_max = atof(argv[15]);
double ratio_range = atof(argv[16]);
bool xlims_set = true;
if (xmin < 0 || xmax < 0)
	xlims_set = false;

TString dir(argv[17]);
TString dtag1(argv[18]);

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

TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);

//TLegend *leg = new TLegend(0.845, 0.2, 0.99, 0.99);
//leg = new TLegend(0.845, 0.2, 0.99, 0.99);
TLegend* leg = new TLegend(0.7, 0.7, 0.89, 0.89);

for (int i = DTAG_ARGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	//cout << "processing " << dtag << endl;
	dtags.push_back(dtag);
	files.push_back(TFile::Open(dir + "/" + dtag + ".root"));
	//dtags.push_back(5);
	bool isData = dtag.Contains("Data");
	TString distr = (isData? distr_data : distr_mc);

	if (!files.back()->GetListOfKeys()->Contains(distr))
		{
		cout << "no " << distr << endl;
		continue;
		}

	if (files.back()->GetListOfKeys()->Contains("weightflow_el_NOMINAL"))
		weightflows.push_back((TH1D*) files.back()->Get("weightflow_el_NOMINAL"));
	else if (files.back()->GetListOfKeys()->Contains("eventflow"))
		weightflows.push_back((TH1D*) files.back()->Get("eventflow"));
	else if (files.back()->GetListOfKeys()->Contains("weightflow_elel_NOMINAL"))
		weightflows.push_back((TH1D*) files.back()->Get("weightflow_elel_NOMINAL"));
	else
		{
		cout << "no weightflow distr" << endl;
		return 2;
		}
	//weightflows.back()->Print();

	TH1D* histo = (TH1D*) files.back()->Get(distr);

	if (!isData)
		{
		Double_t ratio = 1;
		if (inclusive_xsec && dtag.Contains("TTJ"))
			ratio = lumi * xsecs_inclusive_tt[dtag] / weightflows.back()->GetBinContent(11);
		else
			ratio = lumi * xsecs[dtag] / weightflows.back()->GetBinContent(11);
		histo->Scale(ratio);
		//cout << "scaling and adding a stack histo " << dtag << " ratio = " << ratio << endl;
		cout << dtag << "\t" << ratio << "\t" << histo->Integral() << endl;
		//histo->Print();
		//histos.back()->SetFillColor(kRed );
		}
	if (rebin_factor != 1) histo->Rebin(rebin_factor); // should rebin the new histo inplace
	//cout << "rebined" << endl;

	// and normalize per bin-width:
	if (scale_width)
		for (Int_t i=0; i<=histo->GetSize(); i++)
			{
			//yAxis->GetBinLowEdge(3)
			double content = histo->GetBinContent(i);
			double error   = histo->GetBinError(i);
			double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
			histo->SetBinContent(i, content/width);
			histo->SetBinError(i, error/width);
			}

	//cout << "scaled bin widths" << endl;

	histos.push_back(histo);

	//histos.back()->Print();
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

		//hs->Add(histos.back(), "HIST");
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

	cout << " adding to mc stack " << distr->GetBinContent(1) << "\t" << distr->GetBinContent(2) << "\t" << distr->GetFillColor() << ":" << kAzure << endl;
	distr->Print();

	hs->Add(distr, "HIST");

	//TLegendEntry *entry=leg->AddEntry("NULL","","h");
	//entry=leg->AddEntry("singlemu_altstep4rho3","Single top","F");
	leg->AddEntry(distr, nick, "F");
	}

cout << "mc sum " << hs_sum->GetBinContent(1) << "\t" << hs_sum->GetBinContent(2) << "\t" << endl;

/* scale the stack and the sum
*/
// normalize MC sum to Data
if (normalize_to_data)
	{
	cout << "normalizing MC" << endl;
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

/* I don't get how root chooses default values
 * it constantly screwws everything up
 */
if (y_min > 0 && y_max > 0)
	{
	hs_data->GetYaxis()->    SetRange(y_min, y_max);
	hs_data->GetYaxis()->SetRangeUser(y_min, y_max);
	//hs->GetYaxis()->    SetRange(0, 5500);
	//hs->GetYaxis()->SetRangeUser(0, 5500);
	hs_sum->GetYaxis()->    SetRange(y_min, y_max);
	hs_sum->GetYaxis()->SetRangeUser(y_min, y_max);
	}

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
cout << "relative distr" << endl;

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

cout << "set labels etc" << endl;

Int_t HS_size = hs_sum_relative->GetSize();
cout << "histo size " << HS_size << endl;
for (int i=0; i<=HS_size; i++)
	{
	double mc_content = hs_sum_relative->GetBinContent(i);
	double mc_error   = hs_sum_relative->GetBinError(i);
	//double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);

	double data_content = hs_data_relative->GetBinContent(i);
	double data_error   = hs_data_relative->GetBinError(i);
	//double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
	//cout << "got contents" << endl;
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
	//cout << "reset contents, now size is " << hs_sum_relative->GetSize() << " i = " << i << endl;
	}
cout << "scaled the distrs" << endl;

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

cout << "drawn" << endl;

cst->Modified();

cst->SaveAs( dir + "/jobsums/" + distr_mc + suffix + (set_logy? "_logy" : "") + (normalize_to_data? "_normalizedToData.png" : ".png") );

return 0;
}


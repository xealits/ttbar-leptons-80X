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
#include "TPaveText.h"

#include <map>
#include <string>
#include <vector>

#include <fstream>

#include "dtag_xsecs.h"

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

// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

/*
double pu_weight(Int_t nvtx) {
   return pileup_ratio[nvtx];
}
*/

double closest_t_mass(
	Int_t jet0_id, Float_t jet0_b_discr, LorentzVector* jet0_p4, 
	Int_t jet1_id, Float_t jet1_b_discr, LorentzVector* jet1_p4, 
	Int_t jet2_id, Float_t jet2_b_discr, LorentzVector* jet2_p4, 
	Int_t jet3_id, Float_t jet3_b_discr, LorentzVector* jet3_p4, 
	Int_t jet4_id, Float_t jet4_b_discr, LorentzVector* jet4_p4
	)
	{
	//return (jet0_pu_discr > pu_threshold ?1:0) + (jet1_pu_discr > pu_threshold ?1:0) + (jet2_pu_discr > pu_threshold ?1:0) + (jet3_pu_discr > pu_threshold ?1:0) + (jet4_pu_discr > pu_threshold ?1:0);
	std::vector<LorentzVector> light_jets, heavy_jets;

	Float_t b_threshold = 0.8484;

	if (jet0_id >=0)
		{
		if (jet0_b_discr < b_threshold)
			light_jets.push_back(LorentzVector(jet0_p4->X(), jet0_p4->Y(), jet0_p4->Z(), jet0_p4->T()));
		else
			heavy_jets.push_back(LorentzVector(jet0_p4->X(), jet0_p4->Y(), jet0_p4->Z(), jet0_p4->T()));
		}
	if (jet1_id >=0)
		{
		if (jet1_b_discr < b_threshold)
			light_jets.push_back(LorentzVector(jet1_p4->X(), jet1_p4->Y(), jet1_p4->Z(), jet1_p4->T()));
		else
			heavy_jets.push_back(LorentzVector(jet1_p4->X(), jet1_p4->Y(), jet1_p4->Z(), jet1_p4->T()));
		}
	if (jet2_id >=0)
		{
		if (jet2_b_discr < b_threshold)
			light_jets.push_back(LorentzVector(jet2_p4->X(), jet2_p4->Y(), jet2_p4->Z(), jet2_p4->T()));
		else
			heavy_jets.push_back(LorentzVector(jet2_p4->X(), jet2_p4->Y(), jet2_p4->Z(), jet2_p4->T()));
		}
	if (jet3_id >=0)
		{
		if (jet3_b_discr < b_threshold)
			light_jets.push_back(LorentzVector(jet3_p4->X(), jet3_p4->Y(), jet3_p4->Z(), jet3_p4->T()));
		else
			heavy_jets.push_back(LorentzVector(jet3_p4->X(), jet3_p4->Y(), jet3_p4->Z(), jet3_p4->T()));
		}
	if (jet4_id >=0)
		{
		if (jet4_b_discr < b_threshold)
			light_jets.push_back(LorentzVector(jet4_p4->X(), jet4_p4->Y(), jet4_p4->Z(), jet4_p4->T()));
		else
			heavy_jets.push_back(LorentzVector(jet4_p4->X(), jet4_p4->Y(), jet4_p4->Z(), jet4_p4->T()));
		}

	// at least 1 heavy jet and 2 light
	if (!(heavy_jets.size()>0 && light_jets.size()>1)) return 0;

	double closest_W_mass = 0;
	double smallest_W_mass_distance = 999999;
	double closest_T_mass = 0;
	double smallest_T_mass_distance = 999999;

	LorentzVector W(0,0,0,0);
	for (unsigned int i=0; i<light_jets.size()-1; i++)
		for (unsigned int u=i+1; u<light_jets.size(); u++)
			{
			W = light_jets[i] + light_jets[u];
			double W_mass_dist = ( 80 - W.mass());
			if (W_mass_dist < smallest_W_mass_distance)
				{
				smallest_W_mass_distance = W_mass_dist;
				closest_W_mass = W.mass();
				}
			}

	for (unsigned int j=0; j<heavy_jets.size()-1; j++)
		{
		LorentzVector T = W + heavy_jets[j];
		double T_mass_dist = (170 - T.mass());
		if (T_mass_dist < smallest_T_mass_distance)
			{
			smallest_T_mass_distance = T_mass_dist;
			closest_T_mass = T.mass();
			}
		}

	if (closest_T_mass > 1000) closest_T_mass = 1000;
	return closest_T_mass;
	}





double closest_w_mass(
	Int_t jet0_id, Float_t jet0_b_discr, LorentzVector* jet0_p4, 
	Int_t jet1_id, Float_t jet1_b_discr, LorentzVector* jet1_p4, 
	Int_t jet2_id, Float_t jet2_b_discr, LorentzVector* jet2_p4, 
	Int_t jet3_id, Float_t jet3_b_discr, LorentzVector* jet3_p4, 
	Int_t jet4_id, Float_t jet4_b_discr, LorentzVector* jet4_p4
	)
	{
	//return (jet0_pu_discr > pu_threshold ?1:0) + (jet1_pu_discr > pu_threshold ?1:0) + (jet2_pu_discr > pu_threshold ?1:0) + (jet3_pu_discr > pu_threshold ?1:0) + (jet4_pu_discr > pu_threshold ?1:0);
	std::vector<LorentzVector*> light_jets, heavy_jets;

	Float_t b_threshold = 0.8484;

	if (jet0_id >=0)
		{
		if (jet0_b_discr < b_threshold) light_jets.push_back(jet0_p4);
		}
	if (jet1_id >=0)
		{
		if (jet1_b_discr < b_threshold) light_jets.push_back(jet1_p4);
		}
	if (jet2_id >=0)
		{
		if (jet2_b_discr < b_threshold) light_jets.push_back(jet2_p4);
		}
	if (jet3_id >=0)
		{
		if (jet3_b_discr < b_threshold) light_jets.push_back(jet3_p4);
		}
	if (jet4_id >=0)
		{
		if (jet4_b_discr < b_threshold) light_jets.push_back(jet4_p4);
		}

	// at 2 light jets
	if (!(light_jets.size()>1)) return 0;

	double closest_W_mass = 0;
	double smallest_W_mass_distance = 999999;

	LorentzVector W(0,0,0,0);
	for (unsigned int i=0; i<light_jets.size()-1; i++)
		for (unsigned int u=i+1; u<light_jets.size(); u++)
			{
			W = *light_jets[i] + *light_jets[u];
			double W_mass_dist = (80 - W.mass());
			if (W_mass_dist < smallest_W_mass_distance)
				{
				smallest_W_mass_distance = W_mass_dist;
				closest_W_mass = W.mass();
				}
			}

	if (closest_W_mass > 1000) closest_W_mass = 1000;
	return closest_W_mass;
	}

void deep_distribution(TH2D* histo, TTree* NT_output_ttree, bool isMC, double mu_channel, int distr, bool tt_sig, int sig_id)
	{
	#define NTUPLE_INTERFACE_OPEN
	#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput.h"

	cout << "to the loop! " << NT_output_ttree->GetEntries() << endl;
	for (long i = 0; i<NT_output_ttree->GetEntries(); i++)
		{
		NT_output_ttree->GetEntry(i);
		bool pass = true;
		if (mu_channel) pass &= NT_HLT_mu && abs(NT_leps_ID) == 13 ;
		else            pass &= NT_HLT_el && abs(NT_leps_ID) == 11 ;
		pass &= abs(NT_lep0_p4->eta()) < 2.4 && NT_lep0_dxy < 0.01 && NT_lep0_dz < 0.02 && NT_njets > 2 && NT_met_corrected->pt() > 40 && NT_nbjets > 0;
		if (sig_id > 0) // than it is ttbar (weird, yep)
			{
			if (tt_sig) // than select signal
				pass &= (abs(NT_gen_t_w_decay_id * NT_gen_tb_w_decay_id) == sig_id);
			else // select the rest
				pass &= (abs(NT_gen_t_w_decay_id * NT_gen_tb_w_decay_id) != sig_id);
			}
		//else // not ttbar
		//	{
		//	}
		if (!pass) continue;

		double weight = 1;
		if (isMC)
			{
			weight *= pileup_ratio[NT_nvtx];
			}

		double var1 = closest_w_mass(NT_jet0_id, NT_jet0_b_discr, NT_jet0_p4, NT_jet1_id, NT_jet1_b_discr, NT_jet1_p4,
			NT_jet2_id, NT_jet2_b_discr, NT_jet2_p4, NT_jet3_id, NT_jet3_b_discr, NT_jet3_p4, NT_jet4_id, NT_jet4_b_discr, NT_jet4_p4);
		double var2 = closest_t_mass(NT_jet0_id, NT_jet0_b_discr, NT_jet0_p4, NT_jet1_id, NT_jet1_b_discr, NT_jet1_p4,
			NT_jet2_id, NT_jet2_b_discr, NT_jet2_p4, NT_jet3_id, NT_jet3_b_discr, NT_jet3_p4, NT_jet4_id, NT_jet4_b_discr, NT_jet4_p4);
		histo->Fill(var1, var2, weight);
		}
	}

/*
bool cutfilename_1el3j40met1b1t()
	{
	return HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1;
	}

bool cutfilename_1mu3j40met1b1t()
	{
	return HLT_mu && abs(leps_ID) == 13 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0 && tau0_IDlev > 1;
	}
*/

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}


int add_nicknamed_mc_histo(std::map<TString, TH2D *>& nicknamed_mc_distrs, TH2D** mc_sum, TString& nick, TH2D* histo, bool be_verbose) {
	//std::map<TString, TH1D *> nicknamed_mc_distrs;
	if (nicknamed_mc_distrs.find(nick) == nicknamed_mc_distrs.end())
		nicknamed_mc_distrs[nick] = (TH2D*) histo->Clone();
		//nicknamed_mc_distrs[nick] = new TH2D(nick, "");
	else
		nicknamed_mc_distrs[nick]->Add(histo);

	if (be_verbose) cout << "summing mc,";
	// and add to the sum:
	if (*mc_sum == NULL)
		{
		if (be_verbose) cout << " clone,";
		//(*mc_sum) = (TH2D*) histo->Clone();
		(*mc_sum) = new TH2D("mc_sum", ";;", 50, 1, 300, 50, 1, 500);
		//mc_sum->SetName("mc_sum");
		}
	else
		{
		if (be_verbose) cout << " sum,";
		(*mc_sum)->Add(histo);
		}
	if (be_verbose) cout << " done." << endl;

	return 0;
}

using namespace std;


#define INPUT_DTAGS_START 21
const char usage_string[256] = " [--verbose] [--normalize] distr_number sig global_scale tau_ID_SF with_b_SF with_top_pt with_lep_SF with_lep_trig_SF with_PU_weight set_logy unstack lumi_bcdef lumi_gh distr distr_cond range out_name distr_name Y_name dir dtags";

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < INPUT_DTAGS_START)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();
gROOT->ProcessLine(".L pu_weight.C+");
// load shared library with the wrapper for b-tag SF
// tmp/slc6_amd64_gcc530/src/UserCode/ttbar-leptons-80X/src/UserCodettbar-leptons-80X/libUserCodettbar-leptons-80X.so
gSystem->Load("libUserCodettbar-leptons-80X.so");
//gROOT->ProcessLine("#include <vector>");

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

int i = 1;

int distr_number = atoi(argv[input_starts + i++]);

cout << i << ' ';
TString sig(argv[input_starts + i++]);
if (sig != TString("el") && sig != TString("mu"))
	{
	cout << "not supported tt signal channel: " << sig << endl;
	}

float global_scale = atof(argv[input_starts + i++]);

float tau_ID_SF = atof(argv[input_starts + i++]);

bool with_b_SF = false;
TString set_b_SF(argv[input_starts + i++]);
cout << i << ' ';
if (set_b_SF == TString("T") || set_b_SF == TString("Y"))
	{
	with_b_SF = true;
	gROOT->ProcessLine("set_bSF_calibrators()");
	cout << "with b SF weight" << endl;
	}
else
	cout << "NO b SF weight!" << endl;

bool with_Top_PT = false;
TString set_top_pt(argv[input_starts + i++]);
cout << i << ' ';
if (set_top_pt == TString("T") || set_top_pt == TString("Y"))
	{
	with_Top_PT = true;
	cout << "with Top pT SF weight in TTbar MC" << endl;
	}
else
	cout << "NO Top pT SF weight!" << endl;

bool with_lep_SF = false;
TString set_lep_SF(argv[input_starts + i++]);
cout << i << ' ';
if (set_lep_SF == TString("T") || set_lep_SF == TString("Y"))
	{
	with_lep_SF = true;
	cout << "with lep SF weight" << endl;
	}
else
	cout << "NO lep SF weight!" << endl;

bool with_lep_trig_SF = false;
TString set_lep_trig_SF(argv[input_starts + i++]);
cout << i << ' ';
if (set_lep_trig_SF == TString("T") || set_lep_trig_SF == TString("Y"))
	{
	with_lep_trig_SF = true;
	cout << "with lep trig SF weight" << endl;
	}
else
	cout << "NO lep trig SF weight!" << endl;

bool with_PU_weight = false;
TString set_PU(argv[input_starts + i++]);
cout << i << ' ';
if (set_PU == TString("T") || set_PU == TString("Y"))
	{
	with_PU_weight = true;
	cout << "with PU weight" << endl;
	}
else
	cout << "NO PU weight!" << endl;

bool set_logy = false;
TString logy_inp(argv[input_starts + i++]);
if (logy_inp == TString("T") || logy_inp == TString("Y"))
	{
	set_logy = true;
	}

bool unstack = false;
TString unstack_inp(argv[input_starts + i++]);
if (unstack_inp == TString("T") || unstack_inp == TString("Y"))
	{
	unstack = true;
	}

double lumi_bcdef = atof(argv[input_starts + i++]);
double lumi_gh    = atof(argv[input_starts + i++]);
double lumi = lumi_bcdef + lumi_gh;
TString distr(argv[input_starts + i++]);
TString distr_condition_init(argv[input_starts + i++]);
TString range(argv[input_starts + i++]);
string out_name(argv[input_starts + i++]);

TString distr_name(argv[input_starts + i++]);
TString Y_name(argv[input_starts + i++]);
TString dir(argv[input_starts + i++]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

cout << lumi_bcdef << " + " << lumi_gh << " = " << lumi  << endl;
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

std::map<TString, TH2D*> nicknamed_mc_distrs;
TH2D* mc_sum = NULL;

//THStack *hs = new THStack("hs","Stacked 1D histograms");
THStack *hs      = new THStack("hs", "");
TH2D    *hs_data = NULL;

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

bool inclusive_TT = true ;
bool in_proxy   = true;
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
	if (be_verbose) cout << "processing " << dtag << ' ';
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	if (!is_file_exist(the_file.Data()))
		{
		cout << "file " << the_file << "doesn't exist" << endl;
		continue;
		}

	TFile* file = TFile::Open(the_file);

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
		{
		TString draw_command = distr + (in_proxy? "" : ">>h");
		TString condition_command = distr_condition_init;
		cout << "output_ttree->Draw(" << draw_command << ", " << condition_command << ");" << endl;
		TH2D* histo = new TH2D("h_d", ";;", 50, 1, 300, 50, 1, 500);
		//bool deep_distribution(TH1D* histo, TTree* NT_output_ttree, bool isMC, double mu_channel)
		deep_distribution(histo, output_ttree, false, (sig == TString("mu")), distr_number, false, -1);

		//output_ttree->Draw(draw_command, condition_command);
		if (be_verbose) cout << "done drawing, getting the histogram" << endl;
		//if (in_proxy) histo = (TH1D*) gROOT->FindObject("htemp");
		//else          histo = (TH1D*) output_ttree->GetHistogram();
		if (be_verbose) cout << histo->Integral();

		/*
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
		*/

		if (be_verbose) cout << " summing data-stack " << histo->Integral();
		histo->SetMarkerStyle(9);

		if (hs_data == NULL)
			{
			if (be_verbose) cout << " creating data histo" << endl;
			hs_data = (TH2D*) histo->Clone();
			}
		else
			{
			if (be_verbose) cout << " add histo to data histo" << endl;
			hs_data->Add(histo);
			}
		}
	else
		{
		// MC lumi ratio
		Double_t ratio = 1, xsec = 1;
		TString distr_condition = distr_condition_init;

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
		xsec = xsecs[dtag];
		ratio = lumi * xsec / weightflow->GetBinContent(11);

		/*
		// add weights for MC
		TString weight_cond("");
		if (with_PU_weight)
			{
			weight_cond = "(pu_weight(nvtx_gen))*";
			// PU integral:
			weight_cond += "0.986*";
			}

		if (dtag.Contains("amcatnlo"))
			weight_cond += "(aMCatNLO_weight > 0? 1 : -1)*";

		// lepton SF
		if (with_lep_SF)
			{
			TString lepton_SF_call("");
			if (sig == TString("mu"))
				{
				lepton_SF_call.Form("(lepton_muon_SF(abs(lep0_p4.eta()), lep0_p4.pt(), %f, %f))*", lumi_bcdef/lumi, lumi_gh/lumi);
				}
			else // assume single-electron as only alternative
				{
				lepton_SF_call.Form("(lepton_electron_SF(lep0_p4.eta(), lep0_p4.pt()))*");
				}
			weight_cond += lepton_SF_call;
			}

		if (with_lep_trig_SF)
			{
			TString lepton_trig_SF_call("");
			if (sig == TString("mu"))
				{
				lepton_trig_SF_call.Form("(lepton_muon_trigger_SF(abs(lep0_p4.eta()), lep0_p4.pt(), %f, %f))*", lumi_bcdef/lumi, lumi_gh/lumi);
				}
			else
				{
				lepton_trig_SF_call.Form("(lepton_electron_trigger_SF(lep0_p4.eta(), lep0_p4.pt()))*");
				}
			weight_cond += lepton_trig_SF_call;
			}

		// btag SF
		//b_taggin_SF (double jet_pt, double jet_eta, double jet_b_discr, int jet_hadronFlavour, double b_tag_WP);
		if (with_b_SF)
			{
			// first initialize b SF distrs for this dtag:
			TString init_bSF_call = "set_bSF_effs_for_dtag(\"" + dtag + "\");";
			cout << "init b SFs with: " << init_bSF_call << endl;
			gROOT->ProcessLine(init_bSF_call);
			TString b_SF_call("(b_taggin_SF(jet0_p4.pt(), jet0_p4.eta(), jet0_b_discr, jet0_hadronFlavour, 0.8484))*");
			b_SF_call += "(b_taggin_SF(jet1_p4.pt(), jet1_p4.eta(), jet1_b_discr, jet1_hadronFlavour, 0.8484))*";
			b_SF_call += "(b_taggin_SF(jet2_p4.pt(), jet2_p4.eta(), jet2_b_discr, jet2_hadronFlavour, 0.8484))*";
			b_SF_call += "(b_taggin_SF(jet3_p4.pt(), jet3_p4.eta(), jet3_b_discr, jet3_hadronFlavour, 0.8484))*";
			b_SF_call += "(b_taggin_SF(jet4_p4.pt(), jet4_p4.eta(), jet4_b_discr, jet4_hadronFlavour, 0.8484))*";
			weight_cond += b_SF_call;
			}

		// Tau ID SF
		if (tau_ID_SF > 0)
			{
			TString tau_ID_SF_call("");
			tau_ID_SF_call.Form("%f*", tau_ID_SF);
			weight_cond += tau_ID_SF_call;
			}
			
		// NUP cut for W0Jets
		if (dtag.Contains("W0Jet"))
			distr_condition += "&& NUP_gen < 6";
		*/

		// draw distribution
		// in case of inclusive samples (like TTbar) several distributions are drawn -- for each sub-channel
		// thus the loop
		map<TString, TH2D*> histos;
		//map<TString, TString> histo_conditions;


		TString nick("");// = dtag_nick(dtag);
		if (! dtag.Contains("TT"))
			{
			nick = dtag_nick(dtag);
			TString distr_name_for_drawing = "h_" + nick;
			TH2D* histo = new TH2D(distr_name_for_drawing, ";;", 50, 1, 300, 50, 1, 500);
			deep_distribution(histo, output_ttree, true, (sig == TString("mu")), distr_number, false, -1);
			histos[nick] = histo;
			}
		else
			{
			nick = "tt_{other}";
			TString distr_name_for_drawing = "h_" + nick;
			TH2D* histo = new TH2D(distr_name_for_drawing, ";;", 50, 1, 300, 50, 1, 500);
			deep_distribution(histo, output_ttree, true, (sig == TString("mu")), distr_number, false, (sig == TString("mu") ? 13*15 : 11*15));
			histos[nick] = histo;

			nick = (sig==TString("mu") ? "tt_{\\mu\\tau}" : "tt_{e\\tau}");
			distr_name_for_drawing = "h_" + nick;
			histo = new TH2D(distr_name_for_drawing, ";;", 50, 1, 300, 50, 1, 500);
			deep_distribution(histo, output_ttree, true, (sig == TString("mu")), distr_number, true, (sig == TString("mu") ? 13*15 : 11*15));
			histos[nick] = histo;
			}

		/*
		if (! dtag.Contains("TT"))
			{
			TString nick = dtag_nick(dtag);
			histo_conditions[nick] = weight_cond + " (" + distr_condition + ")";
			//output_ttree->Draw(distr + ">>h" + range, weight_cond + " (" + distr_condition + ")");
			//TH1D* histo = (TH1D*) output_ttree->GetHistogram();
			//histos[nick] = histo;
			}
		else
			{
			if (with_Top_PT)
				{
				weight_cond += "(ttbar_pT_SF(gen_t_pt, gen_tb_pt))*";
				}

			TString nick;

			if (!inclusive_TT)
			{
			if (sig == TString("mu"))
				{
				cout << "TT channels" << endl;
				nick = "tt_{other}";
				//histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 13*11 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 13*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 11*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) > 1*15" + ")";
				histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 13*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) > 1*15" + ")";

				//nick = "tt_em";
				//histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*11" + ")";

				nick = "tt_{\\mu\\tau}";
				histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*15" + ")";
				//output_ttree->Draw(distr + ">>h" + range, weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*15" + ")");
				//TH1D* histo = (TH1D*) output_ttree->GetHistogram();
				//histos[nick] = histo;
				}
			if (sig == TString("el"))
				{
				cout << "TT channels" << endl;
				nick = "tt_{other}";
				//histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 13*11 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 13*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 11*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) > 1*15" + ")";
				histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) != 11*15 && abs(gen_t_w_decay_id * gen_tb_w_decay_id) > 1*15" + ")";

				nick = "tt_{e\\tau}";
				histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*15" + ")";
				}

			nick = "tt_lj";
			histo_conditions[nick] = weight_cond + " (" + distr_condition + " && abs(gen_t_w_decay_id * gen_tb_w_decay_id) <= 1*15" + ")";
			}
			else
				{
				nick = "tt";
				histo_conditions[nick] = weight_cond + " (" + distr_condition + ")";
				}
			}
		*/

		//if (be_verbose) cout << distr << '\t' << weight_cond << '\t' << distr_condition << endl;
		//TString distr_name_for_drawing = ""; // FUCK ROOOT THESE FUNCTIONS ARE AWESOMELY USELESS
		//distr_name_for_drawing.Form("h_%s", nick);

		for (std::map<TString, TH2D*>::iterator it = histos.begin(); it != histos.end(); ++it)
			{
			TString nick = it->first;
			TH2D* histo = it->second;
			if (be_verbose) cout << dtag << ' ' << nick << '\t' << histo->Integral() << '\t' << weightflow->GetBinContent(11) << endl;
			if (be_verbose) histo->Print();

			cout << "scaling and adding a stack histo " << dtag << " xsec = " << xsec << " ratio = " << ratio << " norm = " << histo->Integral() / weightflow->GetBinContent(11) << endl;
			histo->Scale(ratio);
			if (be_verbose) histo->Print();

			/*
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
			*/

			if (be_verbose) cout << "summing mc-stack " << histo->Integral() << endl;
			Color_t col  = nick_colour(nick);
			//TString nick = nick_colour.first;
			//Color_t col = nick_colour.second;

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

			//int add_nicknamed_mc_histo(std::map<TString, TH1D *>& nicknamed_mc_distrs, TH1D* mc_sum, TString& nick, TH1D* histo, be_verbose)
			add_nicknamed_mc_histo(nicknamed_mc_distrs, &mc_sum, nick, histo, be_verbose);
			}
		}

	//if (be_verbose) cout << "processed dtag" << endl;
	}

double integral_data = hs_data->Integral();
double integral_MC_taus  = 0.;
double integral_MC_jets  = 0.;
// INCLUSIVE event yield
cout << "------------------ inclusive event yield ------------------" << endl;
cout << "data integral: " << integral_data << endl;
//cout << "MC taus integral = " << integral_MC_taus   << endl;
cout << "MC   integral = " << mc_sum->Integral() << endl;
//double ratio = integral_data / integral_MC;
//cout << "ratio = " << ratio << endl;


/*
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
*/
// no stacks, no drawing, just save them to file
cout << "compare" << endl;
mc_sum->Print();

cout << "switched to the second one" << endl;

/*
 * Do normalize for comparison with cut-n-count
if (global_scale > 0)
	{
	Double_t scale = global_scale*integral_data/hs_sum->Integral();
	hs_sum->Scale(scale);
	for(std::map<TString, TH1D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
		{
		TH1D * distr = it->second;
		distr->Scale(scale);
		}
	}
 */

// ---- Write the separate distributions to file
//TString ntuple_output_filename = outdir + TString(string("/") + dtag_s + string("_") + job_num + string(".root"));
TFile* f_out = TFile::Open(dir + "/jobsums/" + "QuickNtupleDistr2D_" + out_name + ".root", "RECREATE");
//cst->SaveAs( dir + "/jobsums/" + out_name + "_QuickNtupleDistr_" + (normalize_MC ? "_normalized" : "") + (set_logy? "_logy" : "") + (unstack? "_unstacked" : "") + ".png" );

hs_data->Write();
hs_data->SetName(TString("data"));
f_out->Write(TString("data"));

mc_sum->Write();
mc_sum->SetName(TString("mc_sum"));
f_out->Write(TString("mc_sum"));
for(std::map<TString, TH2D*>::iterator it = nicknamed_mc_distrs.begin(); it != nicknamed_mc_distrs.end(); ++it)
	{
	TString nick = it->first;
	TH2D * distr = it->second;
	distr->SetName(nick);
	distr->Write();
	f_out->Write(nick);
	}

f_out->Close();

return 0;
}


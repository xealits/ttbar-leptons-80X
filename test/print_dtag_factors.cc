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

/*
double pu_weight(Int_t nvtx) {
   return pileup_ratio[nvtx];
}
*/

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}


int add_nicknamed_mc_histo(std::map<TString, TH1D *>& nicknamed_mc_distrs, TH1D** mc_sum, TString& nick, TH1D* histo, bool be_verbose) {
	//std::map<TString, TH1D *> nicknamed_mc_distrs;
	if (nicknamed_mc_distrs.find(nick) == nicknamed_mc_distrs.end())
		nicknamed_mc_distrs[nick] = (TH1D*) histo->Clone();
		//nicknamed_mc_distrs[nick] = new TH1D(nick, "");
	else
		nicknamed_mc_distrs[nick]->Add(histo);

	if (be_verbose) cout << "summing mc,";
	// and add to the sum:
	if (*mc_sum == NULL)
		{
		if (be_verbose) cout << " clone,";
		(*mc_sum) = (TH1D*) histo->Clone();
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


#define INPUT_DTAGS_START 2
const char usage_string[256] = " [--verbose] [--normalize] dir dtags";

//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
if (argc < INPUT_DTAGS_START)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

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

cout << i << ' ';

TString dir(argv[input_starts + i++]);
TString dtag1(argv[input_starts + i++]);

cout << dir   << endl;
cout << dtag1 << endl;


/*
 * Get each dtag file in the reduced dir,
 * calculate the factor and print it
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

	if (!isData)
		{
		// MC lumi ratio
		Double_t ratio = 1, xsec = 1, mc_factor = 1;

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
		//ratio = lumi * xsec / weightflow->GetBinContent(11);
		// no lumi, just normalization:
		mc_factor = xsec / weightflow->GetBinContent(11);

		cout << dtag << ',' << xsec << ',' << mc_factor << endl;
		}

	//if (be_verbose) cout << "processed dtag" << endl;
	}

return 0;
}


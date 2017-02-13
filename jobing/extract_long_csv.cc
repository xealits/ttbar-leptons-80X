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

//#include "dtag_xsecs.h" // TODO: maybe pass xsec right here??

#define INPUT_DTAGS_START 4

using namespace std;

int main (int argc, char *argv[])
{
char usage_string[128] = "[--verbose] distr_name distr_nick dir dtags";
if (argc < 5)
	{
	std::cout << "Usage : " << argv[0] << usage_string << std::endl;
	return 1;
	}

gROOT->Reset();

unsigned int input_starts = 0;
bool be_verbose = false;
string inp1(argv[1]), inp2(argv[2]);
if (inp1 == string("--verbose"))
	{
	input_starts += 1;
	be_verbose = true;
	}

if (be_verbose)
	{
	cout << "v:being verbose" << endl;
	cout << "v:options are taken from " << input_starts << endl;
	}

TString distr_name(argv[input_starts + 1]);
TString distr_nick(argv[input_starts + 2]);
TString dir(argv[input_starts + 3]);
TString dtag1(argv[input_starts + INPUT_DTAGS_START]);

if (be_verbose)
	{
	cout << "v:" << distr_name   << endl;
	cout << "v:" << distr_nick   << endl;
	cout << "v:" << dir   << endl;
	cout << "v:" << dtag1 << endl;
	}

/*
 * Get each dtag file in the reduced dir,
 * open the histogram, if it exists
 * loop through bins printing long-format lines
 * dtag,nbin,bin_content
 */
cout << "dtag,nbin," << distr_nick << endl;

for (int i = input_starts + INPUT_DTAGS_START; i<argc; i++)
	{
	TString dtag(argv[i]);
	if (be_verbose) cout << dtag << endl;
	TString the_file = dir + "/" + dtag + ".root";
	if (be_verbose) cout << the_file << endl;
	TFile* file = TFile::Open(the_file);

	if (!file->GetListOfKeys()->Contains(distr_name))
		{
		cout << "v:no " << distr_name << endl;
		continue;
		}
	TH1D* histo = (TH1D*) file->Get(distr_name);

	// print each bin on separate line
	for (Int_t nbin=0; nbin<=histo->GetSize(); nbin++)
		{
		double content = histo->GetBinContent(nbin);
		cout << dtag << ',' << nbin << ',' << content << endl;
		}

	if (be_verbose) cout << "v:processed dtag" << endl;
	}

return 0;
}


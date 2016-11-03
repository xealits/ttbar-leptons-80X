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

int stacked_histo_distr (int argc, char *argv[])
{
if (argc < 4)
	{
	std::cout << "Usage : " << argv[0] << " distr, dir, dtags" << std::endl;
	exit (0);
	}

gROOT->Reset();

TString distr(argv[1]);

cout << distr << endl;

/*
TFile *file = TFile::Open(dir + "/" + dtag + ".root");
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


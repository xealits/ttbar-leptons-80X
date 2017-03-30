int csv_distr_TH1D(TString distr, TString dir, TString dtag, bool print_header=false)
{
gROOT->Reset();

TFile *file = TFile::Open(dir + "/" + dtag + ".root");
TH1D *h = (TH1D*) file->Get(distr);

Int_t size_x = h->GetNbinsX();

if (print_header)
	{
	cout << "dtag,selection";
	for (int x=1; x<size_x; x++)
		{
		double bin_center = h->GetXaxis()->GetBinCenter(x);
		cout << "," << bin_center;
		}
	cout << "\n";
	}
else
	{ // print the contents of the distr
	cout << dtag << "," << distr;
	for (int x=1; x<size_x; x++)
		{
		//double bin_center = h->GetXaxis()->GetBinCenter(x);
		double global_bin = h->GetBin(x);
		cout << "," << h->GetBinContent(global_bin);
		}
	cout << "\n";
	}

return 0;
}

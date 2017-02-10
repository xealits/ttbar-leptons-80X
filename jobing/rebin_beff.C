static Float_t beff_bins_pt[51] = { 0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
	55, 56, 57, 58, 59, 60, 65, 70, 75, 80, 100, 150, 200, 300, 500 }; // 48 bins 49 edges                                                                                                   
static int beff_n_bins_pt = 50;

static Float_t beff_bins_eta[51] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2,
	-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 }; // 50 bins 51 edges
static int beff_n_bins_eta = 50;

int rebin_beff(TString filename, TString output_filename)
{
TFile* file = TFile::Open(filename);
vector<TString> histo_names = {"btag_b_hadronFlavour_candidates", "btag_b_hadronFlavour_candidates_tagged",
	"btag_c_hadronFlavour_candidates", "btag_c_hadronFlavour_candidates_tagged",
	"btag_udsg_hadronFlavour_candidates", "btag_udsg_hadronFlavour_candidates_tagged"};

TFile* file_out = TFile::Open(output_filename, "RECREATE");

for (int i = 0; i<histo_names.size(); i++)
	{
	TString histo_name = histo_names[i];
	TH2 *old = (TH2D*) file->Get(histo_name); //the original histogram
	old->SetName(histo_name + "_old"); // for root name management
	//create a new TH2 with your bin arrays spec

	TH2F *h = new TH2F(histo_name, old->GetTitle(), beff_n_bins_pt, beff_bins_pt, beff_n_bins_eta, beff_bins_eta);
	TAxis *xaxis = old->GetXaxis();
	TAxis *yaxis = old->GetYaxis();
	for (int j=1; j<=yaxis->GetNbins();j++) {
		for (int i=1; i<=xaxis->GetNbins();i++) {
			h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j));
			}
		}

	h->Write();
	file_out->Write();
	}
file_out->Close();

return 0;
}


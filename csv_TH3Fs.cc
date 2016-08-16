
int csv_TH3Fs(TString filename1, TString filename2, TString out_filename) {
	// read histo_name from the two files
	// add the histos
	// save to out_filename

	gROOT->Reset();


	TFile *file1 = TFile::Open(filename1);
	TFile *file2 = TFile::Open(filename2);

	FILE * out_f;
	out_f = fopen(out_filename.Data(), "w");

	TH3F * h1 = (TH3F*) file1->Get("jets_distr");
	cout << "got jets_distr\n";
	TH3F * h2 = (TH3F*) file2->Get("tau_jets_distr");
	cout << "got tau_jets_distr\n";

	h1->Print();
	h2->Print();

	fprintf(out_f, "i_bin,x,y,z,tau_jets_bin,jets_bin");

	for (int i=0; h1->GetSize(); ++i)
		{
		double xcenter1 = h1->GetZaxis()->GetBinCenter(i);
		double xcenter2 = h2->GetZaxis()->GetBinCenter(i);
		double ycenter1 = h1->GetZaxis()->GetBinCenter(i);
		double ycenter2 = h2->GetZaxis()->GetBinCenter(i);
		double zcenter1 = h1->GetZaxis()->GetBinCenter(i);
		double zcenter2 = h2->GetZaxis()->GetBinCenter(i);
		if (xcenter1 != xcenter2) return 1;

		double n_jets = h1->GetBinContent(i);
		double n_tau_jets = h2->GetBinContent(i);
		//fprintf(out_f, "%d,%g,%g,%g,%g,%g\n", i, xcenter1, ycenter1, zcenter1, n_jets, n_tau_jets);
		fprintf(out_f, "%d,%g,%g\n", i, n_jets, n_tau_jets);
		}

	file1->Close();
	file2->Close();
	fclose(out_f);


	return 0;
}


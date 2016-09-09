
int csv_TH3Fs(TString filename, TString out_filename) {
	// read jets_distr and tau_jets_distr from the filename
	// loop through bins, printing the contents in long csv format
	// save to out_filename

	gROOT->Reset();


	TFile *file = TFile::Open(filename);

	FILE * out_f;
	out_f = fopen(out_filename.Data(), "w");

	TH3F * h1 = (TH3F*) file->Get("jets_distr");
	cout << "got jets_distr\n";
	TH3F * h2 = (TH3F*) file->Get("tau_jets_distr");
	cout << "got tau_jets_distr\n";

	h1->Print();
	h2->Print();

	fprintf(out_f, "i_bin,x,y,z,tau_jets_bin,jets_bin\n");

	Int_t size = h1->GetSize();
	bool passed_midwork = 0;
	cout << "the histos size = " << size << " and " << h2->GetSize() << "\n";

	Int_t size_x = h1->GetNbinsX();
	Int_t size_y = h1->GetNbinsY();
	Int_t size_z = h1->GetNbinsZ();
	cout << "the histos size = " << size_x << " , " << size_y << " , " << size_z << "\n";
	for (int iz=0; iz<size_z; ++iz)
	{
	for (int ix=0; ix<size_x; ++ix)
	{
	for (int iy=0; iy<size_y; ++iy)
		{
		/*
		if ((i / size) > 0.5 && ! passed_midwork)
			{
			passed_midwork = 1;
			cout << "Passed midwork!\n";
			}
		*/
		double xcenter1 = (double) h1->GetXaxis()->GetBinCenter(ix);
		double xcenter2 = (double) h2->GetXaxis()->GetBinCenter(ix);
		double ycenter1 = (double) h1->GetYaxis()->GetBinCenter(iy);
		double ycenter2 = (double) h2->GetYaxis()->GetBinCenter(iy);
		double zcenter1 = (double) h1->GetZaxis()->GetBinCenter(iz);
		double zcenter2 = (double) h2->GetZaxis()->GetBinCenter(iz);
		if (xcenter1 != xcenter2) return 1;

		// some brilliance from ROOT histograms (f*cking 3d array!)
		// Int_t bin = h->GetBin(binx,biny,binz);
		Int_t global_bin = h1->GetBin(ix,iy,iz);
		Int_t global_bin2 = h2->GetBin(ix,iy,iz);
		if (global_bin != global_bin2) return 2;

		double n_jets     = h1->GetBinContent( global_bin );
		double n_tau_jets = h2->GetBinContent( global_bin );
		fprintf(out_f, "%d,%g,%g,%g,%g,%g\n", global_bin, xcenter1, ycenter1, zcenter1, n_jets, n_tau_jets);
		//fprintf(out_f, "%d,%g,%g\n", i, n_jets, n_tau_jets);
		}
	}
	}

	file->Close();
	fclose(out_f);


	return 0;
}


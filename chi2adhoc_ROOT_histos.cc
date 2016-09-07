
int chi2adhoc_ROOT_histos(TString filename1, TString filename2 ) {
	// read histo_name from the two files
	// add the histos
	// save to out_filename

	gROOT->Reset();


	TFile *file1 = TFile::Open(filename1);
	TFile *file2 = TFile::Open(filename2);

	TH3F * jets_h1 = (TH3F*) file1->Get("jets_distr");
	TH3F * jets_h2 = (TH3F*) file2->Get("jets_distr");
	jets_h1->Print();
	jets_h2->Print();

	TH3F * taus_h1 = (TH3F*) file1->Get("tau_jets_distr");
	TH3F * taus_h2 = (TH3F*) file2->Get("tau_jets_distr");
	taus_h1->Print();
	taus_h2->Print();

	// normalize:
	//

	jets_h1->Scale( 1/jets_h1->Integral() );
	jets_h2->Scale( 1/jets_h2->Integral() );

	taus_h1->Scale( 1/taus_h1->Integral() );
	taus_h2->Scale( 1/taus_h2->Integral() );

	jets_h1->Print();
	jets_h2->Print();
	taus_h1->Print();
	taus_h2->Print();


	Int_t size_x = jets_h1->GetNbinsX();
	Int_t size_y = jets_h1->GetNbinsY();
	Int_t size_z = jets_h1->GetNbinsZ();

	Double_t jets_chi2 = 0;
	Double_t taus_chi2 = 0;

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
		/*
		double xcenter1 = (double) h1->GetXaxis()->GetBinCenter(ix);
		double xcenter2 = (double) h2->GetXaxis()->GetBinCenter(ix);
		double ycenter1 = (double) h1->GetYaxis()->GetBinCenter(iy);
		double ycenter2 = (double) h2->GetYaxis()->GetBinCenter(iy);
		double zcenter1 = (double) h1->GetZaxis()->GetBinCenter(iz);
		double zcenter2 = (double) h2->GetZaxis()->GetBinCenter(iz);
		if (xcenter1 != xcenter2) return 1;
		*/

		// some brilliance from ROOT histograms (f*cking 3d array!)
		// Int_t bin = h->GetBin(binx,biny,binz);
		Int_t global_bin  = jets_h1->GetBin(ix,iy,iz);
		Int_t global_bin2 = jets_h2->GetBin(ix,iy,iz);
		if (global_bin != global_bin2) return 2;

		Double_t n_jets1     = jets_h1->GetBinContent( global_bin );
		Double_t n_jets2     = jets_h2->GetBinContent( global_bin );
		Double_t n_taus1 = taus_h1->GetBinContent( global_bin );
		Double_t n_taus2 = taus_h2->GetBinContent( global_bin );
		//fprintf(out_f, "%d,%g,%g,%g,%g,%g\n", global_bin, xcenter1, ycenter1, zcenter1, n_jets, n_tau_jets);
		//fprintf(out_f, "%d,%g,%g\n", i, n_jets, n_tau_jets);

		if ((n_jets1 + n_jets2) != 0) jets_chi2 += 2 * (n_jets1 - n_jets2)*(n_jets1 - n_jets2)/(n_jets1 + n_jets2);
		if ((n_taus1 + n_taus2) != 0) taus_chi2 += 2 * (n_taus1 - n_taus2)*(n_taus1 - n_taus2)/(n_taus1 + n_taus2);

		}
	}
	}

	cout << "jets chi2 = " << jets_chi2 << "     " << "taus chi2 = " << taus_chi2 << "\n";

	file1->Close();
	file2->Close();


	return 0;
}


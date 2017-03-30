
void chi2_ROOT_histos( TString filename1, TString filename2 ) {
	// read jet & tau histos from the two files
	// normalize histos
	// run Chi2Test on them, print results

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

	Double_t chi1 = jets_h1->Chi2Test(jets_h2,"WW");
	Double_t chi2 = taus_h1->Chi2Test(taus_h2,"WW");
	cout << "jets chi2 = " << chi1 << "   taus chi2 = " << chi2 << "\n";
	printf("%g %g\n", chi1, chi2);

	file1->Close();
	file2->Close();
}


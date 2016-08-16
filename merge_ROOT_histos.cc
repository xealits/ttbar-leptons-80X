
void merge_ROOT_histos( TString filename1, TString filename2, TString out_filename ) {
	// read histo_name from the two files
	// add the histos
	// save to out_filename

	gROOT->Reset();


	TFile *file1 = TFile::Open(filename1);
	TFile *file2 = TFile::Open(filename2);

	TH3F * jets_h1 = (TH3F*) file1->Get("jets_distr");
	TH3F * jets_h2 = (TH3F*) file2->Get("jets_distr");

	TH3F * taus_h1 = (TH3F*) file1->Get("tau_jets_distr");
	TH3F * taus_h2 = (TH3F*) file2->Get("tau_jets_distr");

	jets_h1->Print();
	jets_h2->Print();
	jets_h1->Add(jets_h2, 1);
	jets_h1->Print();

	taus_h1->Print();
	taus_h2->Print();
	taus_h1->Add(taus_h2, 1);
	taus_h1->Print();

	TFile *out_f = TFile::Open(out_filename, "RECREATE");
	jets_h1->Write();
	taus_h1->Write();
	out_f->Write();
	out_f->Close();

	file1->Close();
	file2->Close();
}


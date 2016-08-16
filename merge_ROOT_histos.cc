
void merge_ROOT_histos(TString histo_name, TString filename1, TString filename2, TString out_filename) {
	// read histo_name from the two files
	// add the histos
	// save to out_filename

	gROOT->Reset();


	TFile *file1 = TFile::Open(filename1);
	TFile *file2 = TFile::Open(filename2);

	TH3F * h1 = (TH3F*) file1->Get(histo_name);
	TH3F * h2 = (TH3F*) file2->Get(histo_name);

	h1->Print();
	h2->Print();

	h1->Add(h2, 1);

	h1->Print();

	TFile *out_f = TFile::Open(out_filename, "RECREATE");
	h1->Write();
	out_f->Write();
	out_f->Close();

	file1->Close();
	file2->Close();
}


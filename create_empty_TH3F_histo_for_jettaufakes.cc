
void create_empty_TH3F_histo_for_jettaufakes(TString histo_name, TString filename, TString out_filename)
	{
	// clone empty copy of TH3F histo_name from filename to out_file

	gROOT->Reset();


	TFile *file = TFile::Open(filename);

	TH3F * h = (TH3F*) file->Get(histo_name);
	TH3F * copy_h = (TH3F*) h->Clone(); // copy
	copy_h->Reset(); // make it empty

	h->Print();
	copy_h->Print();

	TFile *out_f = TFile::Open(out_filename, "CREATE");
	copy_h->Write();
	out_f->Write();
	out_f->Close();

	file->Close();
	}


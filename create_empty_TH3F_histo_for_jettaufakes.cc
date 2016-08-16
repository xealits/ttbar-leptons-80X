
void create_empty_TH3F_histo_for_jettaufakes( TString filename, TString out_filename)
	{
	// clone empty copy of TH3F histo_name from filename to out_file

	gROOT->Reset();


	TFile *file = TFile::Open(filename);

	TH3F * jets_h = (TH3F*) file->Get("jets_distr");
	TH3F * copy_jets_h = (TH3F*) jets_h->Clone(); // copy
	copy_jets_h->Reset(); // make it empty

	jets_h->Print();
	copy_jets_h->Print();

	TH3F * tau_jets_h = (TH3F*) file->Get("tau_jets_distr");
	TH3F * copy_tau_jets_h = (TH3F*) tau_jets_h->Clone(); // copy
	copy_tau_jets_h->Reset(); // make it empty

	tau_jets_h->Print();
	copy_tau_jets_h->Print();

	TFile *out_f = TFile::Open(out_filename, "CREATE");
	copy_jets_h->Write();
	copy_tau_jets_h->Write();
	out_f->Write();
	out_f->Close();

	file->Close();
	}



/*
 *   KEY: TH3F	qcd_jets_distr;1	
 *   KEY: TH3F	qcd_tau_jets_distr;1	
 *   KEY: TH1D	qcd_taujet_distance;1	
 *   KEY: TH1I	qcd_taujet_origin;1	
 *   KEY: TH1I	qcd_jet_origin;1	
 *   KEY: TH3F	wjets_jets_distr;1	
 *   KEY: TH3F	wjets_tau_jets_distr;1	
 *   KEY: TH1D	wjets_taujet_distance;1	
 *   KEY: TH1I	wjets_taujet_origin;1	
 *   KEY: TH1I	wjets_jet_origin;1	
 *
 */
void create_empty_TH3F_histo_for_jettaufakes( TString filename, TString out_filename)
	{
	// clone empty copy of TH3F histo_name from filename to out_file

	gROOT->Reset();


	TFile *file = TFile::Open(filename);

	//TH3F * jets_h = (TH3F*) file->Get("jets_distr");
	TH3F * qcd_jets_h   = (TH3F*) file->Get("qcd_jets_distr");
	TH3F * wjets_jets_h = (TH3F*) file->Get("wjets_jets_distr");
	TH3F * copy_q_jets_h = (TH3F*)   qcd_jets_h->Clone(); // copy
	TH3F * copy_w_jets_h = (TH3F*) wjets_jets_h->Clone(); // copy
	copy_q_jets_h->Reset(); // make it empty
	copy_w_jets_h->Reset(); // make it empty

	qcd_jets_h->Print();
	copy_q_jets_h->Print();

	//TH3F * tau_jets_h = (TH3F*) file->Get("tau_jets_distr");
	TH3F * q_tau_jets_h = (TH3F*) file->Get("qcd_tau_jets_distr");
	TH3F * w_tau_jets_h = (TH3F*) file->Get("wjets_tau_jets_distr");
	TH3F * copy_q_tau_jets_h = (TH3F*) q_tau_jets_h->Clone(); // copy
	TH3F * copy_w_tau_jets_h = (TH3F*) w_tau_jets_h->Clone(); // copy
	copy_q_tau_jets_h->Reset(); // make it empty
	copy_w_tau_jets_h->Reset(); // make it empty

	q_tau_jets_h->Print();
	copy_q_tau_jets_h->Print();

	//fake jet distance
	TH3F * q_dist_h = (TH3F*) file->Get("qcd_taujet_distance");
	TH3F * w_dist_h = (TH3F*) file->Get("wjets_taujet_distance");
	TH3F * copy_q_dist_h = (TH3F*) q_dist_h->Clone(); // copy
	TH3F * copy_w_dist_h = (TH3F*) w_dist_h->Clone(); // copy
	copy_q_dist_h->Reset(); // make it empty
	copy_w_dist_h->Reset(); // make it empty

	//jet origins
	//*   KEY: TH1I	qcd_taujet_origin;1	
	//*   KEY: TH1I	qcd_jet_origin;1	
	TH1I * q_jet_origin_h = (TH1I*) file->Get("qcd_jet_origin");
	TH1I * q_taujet_origin_h = (TH1I*) file->Get("qcd_taujet_origin");
	TH1I * w_jet_origin_h = (TH1I*) file->Get("wjets_jet_origin");
	TH1I * w_taujet_origin_h = (TH1I*) file->Get("wjets_taujet_origin");

	TH1I * copy_q_jet_origin_h = (TH1I*) q_jet_origin_h->Clone(); // copy
	TH1I * copy_q_taujet_origin_h = (TH1I*) q_taujet_origin_h->Clone(); // copy
	TH1I * copy_w_jet_origin_h = (TH1I*) w_jet_origin_h->Clone(); // copy
	TH1I * copy_w_taujet_origin_h = (TH1I*) w_taujet_origin_h->Clone(); // copy
	copy_q_jet_origin_h->Reset();
	copy_q_taujet_origin_h->Reset();
	copy_w_jet_origin_h->Reset();
	copy_w_taujet_origin_h->Reset();

	TFile *out_f = TFile::Open(out_filename, "CREATE");
	copy_q_jets_h->Write();
	copy_q_tau_jets_h->Write();
	copy_q_dist_h->Write();
	copy_w_jets_h->Write();
	copy_w_tau_jets_h->Write();
	copy_w_dist_h->Write();

	copy_q_jet_origin_h->Write();
	copy_q_taujet_origin_h->Write();
	copy_w_jet_origin_h->Write();
	copy_w_taujet_origin_h->Write();

	out_f->Write();
	out_f->Close();

	file->Close();
	}



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
	/*
	*/

	/*
	 *   KEY: TH3F	qcd jets_distr;1	
	 *   KEY: TH3F	qcd tau_jets_distr;1	
	 *   KEY: TH1D	qcd taujet_distance;1	
	 *   KEY: TH3F	wjets jets_distr;1	
	 *   KEY: TH3F	wjets tau_jets_distr;1	
	 *   KEY: TH1D	wjets taujet_distance;1	
	 */

	TH3F * qcd_jets_h1 = (TH3F*) file1->Get("qcd_jets_distr");
	TH3F * qcd_jets_h2 = (TH3F*) file2->Get("qcd_jets_distr");

	TH3F * qcd_taus_h1 = (TH3F*) file1->Get("qcd_tau_jets_distr");
	TH3F * qcd_taus_h2 = (TH3F*) file2->Get("qcd_tau_jets_distr");

	TH3F * qcd_dist_h1 = (TH3F*) file1->Get("qcd_taujet_distance");
	TH3F * qcd_dist_h2 = (TH3F*) file2->Get("qcd_taujet_distance");

	TH3F * wjets_jets_h1 = (TH3F*) file1->Get("wjets_jets_distr");
	TH3F * wjets_jets_h2 = (TH3F*) file2->Get("wjets_jets_distr");

	TH3F * wjets_taus_h1 = (TH3F*) file1->Get("wjets_tau_jets_distr");
	TH3F * wjets_taus_h2 = (TH3F*) file2->Get("wjets_tau_jets_distr");

	TH3F * wjets_dist_h1 = (TH3F*) file1->Get("wjets_taujet_distance");
	TH3F * wjets_dist_h2 = (TH3F*) file2->Get("wjets_taujet_distance");

	TH1I * q_jet_origin_h1 = (TH1I*) file1->Get("qcd_jet_origin");
	TH1I * q_taujet_origin_h1 = (TH1I*) file1->Get("qcd_taujet_origin");
	TH1I * w_jet_origin_h1 = (TH1I*) file1->Get("wjets_jet_origin");
	TH1I * w_taujet_origin_h1 = (TH1I*) file1->Get("wjets_taujet_origin");

	TH1I * q_jet_origin_h2 = (TH1I*) file2->Get("qcd_jet_origin");
	TH1I * q_taujet_origin_h2 = (TH1I*) file2->Get("qcd_taujet_origin");
	TH1I * w_jet_origin_h2 = (TH1I*) file2->Get("wjets_jet_origin");
	TH1I * w_taujet_origin_h2 = (TH1I*) file2->Get("wjets_taujet_origin");

	qcd_jets_h1->Print();
	qcd_jets_h2->Print();
	qcd_jets_h1->Add(qcd_jets_h2, 1);
	qcd_jets_h1->Print();

	qcd_taus_h1->Print();
	qcd_taus_h2->Print();
	qcd_taus_h1->Add(qcd_taus_h2, 1);
	qcd_taus_h1->Print();

	wjets_jets_h1->Print();
	wjets_jets_h2->Print();
	wjets_jets_h1->Add(wjets_jets_h2, 1);
	wjets_jets_h1->Print();

	wjets_taus_h1->Print();
	wjets_taus_h2->Print();
	wjets_taus_h1->Add(wjets_taus_h2, 1);
	wjets_taus_h1->Print();

	//qcd_dist_h1->Print();
	//qcd_dist_h2->Print();
	qcd_dist_h1->Add(qcd_dist_h2, 1);
	//qcd_dist_h1->Print();

	//wjets_dist_h1->Print();
	//wjets_dist_h2->Print();
	wjets_dist_h1->Add(wjets_dist_h2, 1);
	//wjets_dist_h1->Print();
/*
*/

	q_jet_origin_h1->Print();
	w_jet_origin_h1->Print();
	q_jet_origin_h1->Add(q_jet_origin_h2, 1);
	q_taujet_origin_h1->Add(q_taujet_origin_h2, 1);
	w_jet_origin_h1->Add(w_jet_origin_h2, 1);
	w_taujet_origin_h1->Add(w_taujet_origin_h2, 1);
	q_jet_origin_h1->Print();
	w_jet_origin_h1->Print();

	TFile *out_f = TFile::Open(out_filename, "RECREATE");
	qcd_jets_h1->Write();
	qcd_taus_h1->Write();
	wjets_jets_h1->Write();
	wjets_taus_h1->Write();
	qcd_dist_h1->Write();
	wjets_dist_h1->Write();

	q_jet_origin_h1->Write();
	q_taujet_origin_h1->Write();
	w_jet_origin_h1->Write();
	w_taujet_origin_h1->Write();

	out_f->Write();
	out_f->Close();

	file1->Close();
	file2->Close();
}


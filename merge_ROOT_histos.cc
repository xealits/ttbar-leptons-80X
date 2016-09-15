
void merge_ROOT_histos( TString filename1, TString filename2, TString out_filename ) {
	// read histo_name from the two files
	// add the histos
	// save to out_filename

	gROOT->Reset();


	TFile *file1 = TFile::Open(filename1);
	TFile *file2 = TFile::Open(filename2);

	/*
	TH3F * jets_h1 = (TH3F*) file1->Get("jets_distr");
	TH3F * jets_h2 = (TH3F*) file2->Get("jets_distr");

	TH3F * taus_h1 = (TH3F*) file1->Get("tau_jets_distr");
	TH3F * taus_h2 = (TH3F*) file2->Get("tau_jets_distr");
	*/

	/*
	 *   KEY: TH3F	qcd jets_distr;1	
	 *   KEY: TH3F	qcd tau_jets_distr;1	
	 *   KEY: TH1D	qcd taujet_distance;1	
	 *   KEY: TH3F	wjets jets_distr;1	
	 *   KEY: TH3F	wjets tau_jets_distr;1	
	 *   KEY: TH1D	wjets taujet_distance;1	
	 */

	TH3F * qcd_jets_h1 = (TH3F*) file1->Get("qcd jets_distr");
	TH3F * qcd_jets_h2 = (TH3F*) file2->Get("qcd jets_distr");

	TH3F * qcd_taus_h1 = (TH3F*) file1->Get("qcd tau_jets_distr");
	TH3F * qcd_taus_h2 = (TH3F*) file2->Get("qcd tau_jets_distr");

	TH3F * qcd_dist_h1 = (TH3F*) file1->Get("qcd taujet_distance");
	TH3F * qcd_dist_h2 = (TH3F*) file2->Get("qcd taujet_distance");

	TH3F * wjets_jets_h1 = (TH3F*) file1->Get("wjets jets_distr");
	TH3F * wjets_jets_h2 = (TH3F*) file2->Get("wjets jets_distr");

	TH3F * wjets_taus_h1 = (TH3F*) file1->Get("wjets tau_jets_distr");
	TH3F * wjets_taus_h2 = (TH3F*) file2->Get("wjets tau_jets_distr");

	TH3F * wjets_dist_h1 = (TH3F*) file1->Get("wjets taujet_distance");
	TH3F * wjets_dist_h2 = (TH3F*) file2->Get("wjets taujet_distance");

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

	TFile *out_f = TFile::Open(out_filename, "RECREATE");
	qcd_jets_h1->Write();
	qcd_taus_h1->Write();
	wjets_jets_h1->Write();
	wjets_taus_h1->Write();
	qcd_dist_h1->Write();
	wjets_dist_h1->Write();
	out_f->Write();
	out_f->Close();

	file1->Close();
	file2->Close();
}


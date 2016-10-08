
int csv_qw_TH3Fs(TString filename, TString out_filename) {
	// read jets_distr and tau_jets_distr from the filename
	// loop through bins, printing the contents in long csv format
	// save to out_filename

	gROOT->Reset();

	// HLTjet_qcd_jets_distr;1	
	// HLTjet_qcd_tau_jets_distr;1	
	// HLTjet_wjets_jets_distr;1	
	// HLTjet_wjets_tau_jets_distr;1	
	// HLTjetmu_qcd_jets_distr;1	
	// HLTjetmu_qcd_tau_jets_distr;1	
	// HLTjetmu_wjets_jets_distr;1	
	// HLTjetmu_wjets_tau_jets_distr;1	

	Float_t bins_pt[11] = { 0, 30, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
	//Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
 	Float_t bins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
	int n_bins_eta = 7;
	Float_t bins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	 	0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges

	TFile *file = TFile::Open(filename);


	TH3F * all_h[12];
	TH3F * empty_h = (TH3F*) new TH3F("empty",      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad);

	TH3F * hm_q_h1 = (TH3F*) file->Get("HLTmu_qcd_jets_distr");
	all_h[0] = (hm_q_h1) ? hm_q_h1 : empty_h;
	TH3F * hm_q_h2 = (TH3F*) file->Get("HLTmu_qcd_tau_jets_distr");
	all_h[1] = (hm_q_h2) ? hm_q_h2 : empty_h;
	TH3F * hm_w_h1 = (TH3F*) file->Get("HLTmu_wjets_jets_distr");
	all_h[2] = (hm_w_h1) ? hm_w_h1 : empty_h;
	TH3F * hm_w_h2 = (TH3F*) file->Get("HLTmu_wjets_tau_jets_distr");
	all_h[3] = (hm_w_h2) ? hm_w_h2 : empty_h;

	TH3F * hj_q_h1 = (TH3F*) file->Get("HLTjet_qcd_jets_distr");
	all_h[4] = (hj_q_h1) ? hj_q_h1 : empty_h;
	TH3F * hj_q_h2 = (TH3F*) file->Get("HLTjet_qcd_tau_jets_distr");
	all_h[5] = (hj_q_h2) ? hj_q_h2 : empty_h;
	TH3F * hj_w_h1 = (TH3F*) file->Get("HLTjet_wjets_jets_distr");
	all_h[6] = (hj_w_h1) ? hj_w_h1 : empty_h;
	TH3F * hj_w_h2 = (TH3F*) file->Get("HLTjet_wjets_tau_jets_distr");
	all_h[7] = (hj_w_h2) ? hj_w_h2 : empty_h;

	TH3F * hjm_q_h1 = (TH3F*) file->Get("HLTjetmu_qcd_jets_distr");
	all_h[8] = (hjm_q_h1) ? hjm_q_h1 : empty_h;
	TH3F * hjm_q_h2 = (TH3F*) file->Get("HLTjetmu_qcd_tau_jets_distr");
	all_h[9] = (hjm_q_h2) ? hjm_q_h2 : empty_h;
	TH3F * hjm_w_h1 = (TH3F*) file->Get("HLTjetmu_wjets_jets_distr");
	all_h[10] = (hjm_w_h1) ? hjm_w_h1 : empty_h;
	TH3F * hjm_w_h2 = (TH3F*) file->Get("HLTjetmu_wjets_tau_jets_distr");
	all_h[11] = (hjm_w_h2) ? hjm_w_h2 : empty_h;

	TH3F * q_hd = (TH3F*) file->Get("qcd_taujet_distance");
	TH3F * hd = (TH3F*) file->Get("wjets_taujet_distance");

	for (int i = 0; i < 12; i++)
		{
		all_h[i]->Print();
		}
	//all_h[0]->Print();
	//all_h[1]->Print();
	//hj_q_h1->Print();
	//hj_q_h2->Print();
	hd->Print();


	FILE * out_f;
	out_f = fopen(out_filename.Data(), "w");
	fprintf(out_f, "i_bin,x,y,z");
	fprintf(out_f, " hm_q_jets_bin,hm_q_tau_jets_bin, hm_w_jets_bin,hm_w_tau_jets_bin");
	fprintf(out_f, " hj_q_jets_bin,hj_q_tau_jets_bin, hj_w_jets_bin,hj_w_tau_jets_bin");
	fprintf(out_f, " hjm_q_jets_bin,hjm_q_tau_jets_bin, hjm_w_jets_bin,hjm_w_tau_jets_bin\n");

	Int_t size = all_h[0]->GetSize();
	bool passed_midwork = 0;
	cout << "the histos size = " << size << " and " << all_h[0]->GetSize() << "\n";

	Int_t size_x = all_h[0]->GetNbinsX();
	Int_t size_y = all_h[0]->GetNbinsY();
	Int_t size_z = all_h[0]->GetNbinsZ();
	cout << "the histos size = " << size_x << " , " << size_y << " , " << size_z << "\n";

	for (int iz=1; iz<=size_z; ++iz)
	{
	for (int ix=1; ix<=size_x; ++ix)
	{
	for (int iy=1; iy<=size_y; ++iy)
		{
		//cout << "AAA\n";

		double xcenter1 = (double) all_h[0]->GetXaxis()->GetBinCenter(ix);
		double xcenter2 = (double) all_h[1]->GetXaxis()->GetBinCenter(ix);
		double ycenter1 = (double) all_h[0]->GetYaxis()->GetBinCenter(iy);
		double ycenter2 = (double) all_h[1]->GetYaxis()->GetBinCenter(iy);
		double zcenter1 = (double) all_h[0]->GetZaxis()->GetBinCenter(iz);
		double zcenter2 = (double) all_h[1]->GetZaxis()->GetBinCenter(iz);
		if (xcenter1 != xcenter2) return 1;

		// some brilliance from ROOT histograms (f*cking 3d array!)
		// Int_t bin = h->GetBin(binx,biny,binz);
		Int_t global_bin = all_h[0]->GetBin(ix,iy,iz);
		Int_t global_bin2 = all_h[1]->GetBin(ix,iy,iz);
		if (global_bin != global_bin2) return 2;
		// TODO: check the bins are the same in qcd jets?

		//cout << "BBB\n";

		double hm_q_n_jets     = all_h[0]->GetBinContent( global_bin );
		double hm_q_n_tau_jets = all_h[1]->GetBinContent( global_bin );
		double hm_w_n_jets     = all_h[2]->GetBinContent( global_bin );
		double hm_w_n_tau_jets = all_h[3]->GetBinContent( global_bin );

		double hj_q_n_jets     = all_h[4]->GetBinContent( global_bin );
		double hj_q_n_tau_jets = all_h[5]->GetBinContent( global_bin );
		double hj_w_n_jets     = all_h[6]->GetBinContent( global_bin );
		double hj_w_n_tau_jets = all_h[7]->GetBinContent( global_bin ); // no histo in j5.6 Data
		//double hj_w_n_tau_jets = 0;
		//double hj_q_n_tau_jets = 0;

		double hjm_q_n_jets     = all_h[8]->GetBinContent( global_bin );
		double hjm_q_n_tau_jets = all_h[9]->GetBinContent( global_bin );
		double hjm_w_n_jets     = all_h[10]->GetBinContent( global_bin );
		double hjm_w_n_tau_jets = all_h[11]->GetBinContent( global_bin );

		//cout << "CCC\n";

		fprintf(out_f, "%d, %g,%g,%g", global_bin, xcenter1, ycenter1, zcenter1);
		fprintf(out_f, " ,%g,%g,%g,%g", hm_q_n_jets,hm_q_n_tau_jets,   hm_w_n_jets,hm_w_n_tau_jets);
		fprintf(out_f, " ,%g,%g,%g,%g", hj_q_n_jets,hj_q_n_tau_jets,   hj_w_n_jets,hj_w_n_tau_jets);
		fprintf(out_f, " ,%g,%g,%g,%g", hjm_q_n_jets,hjm_q_n_tau_jets, hjm_w_n_jets,hjm_w_n_tau_jets);
		fprintf(out_f, "\n");
		//fprintf(out_f, "%d,%g,%g\n", i, n_jets, n_tau_jets);
		//cout << "DDD\n";
		}
	}
	}

	fclose(out_f);
	file->Close();


	return 0;
}



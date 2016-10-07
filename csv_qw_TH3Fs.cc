
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


	TFile *file = TFile::Open(filename);

	FILE * out_f;
	out_f = fopen(out_filename.Data(), "w");

	TH3F * hj_q_h1 = (TH3F*) file->Get("HLTjet_qcd_jets_distr");
	TH3F * hj_q_h2 = (TH3F*) file->Get("HLTjet_qcd_tau_jets_distr");
	TH3F * hj_w_h1 = (TH3F*) file->Get("HLTjet_wjets_jets_distr");
	TH3F * hj_w_h2 = (TH3F*) file->Get("HLTjet_wjets_tau_jets_distr");
	TH3F * hjm_q_h1 = (TH3F*) file->Get("HLTjetmu_qcd_jets_distr");
	TH3F * hjm_q_h2 = (TH3F*) file->Get("HLTjetmu_qcd_tau_jets_distr");
	TH3F * hjm_w_h1 = (TH3F*) file->Get("HLTjetmu_wjets_jets_distr");
	TH3F * hjm_w_h2 = (TH3F*) file->Get("HLTjetmu_wjets_tau_jets_distr");

	TH3F * q_hd = (TH3F*) file->Get("qcd_taujet_distance");
	TH3F * hd = (TH3F*) file->Get("wjets_taujet_distance");

	hj_q_h1->Print();
	hj_q_h2->Print();
	hd->Print();

	fprintf(out_f, "i_bin,x,y,z, hj_w_jets_bin,hj_w_tau_jets_bin, hj_q_jets_bin,hj_q_tau_jets_bin, hjm_w_jets_bin,hjm_w_tau_jets_bin, hjm_q_jets_bin,hjm_q_tau_jets_bin\n");

	Int_t size = hj_q_h1->GetSize();
	bool passed_midwork = 0;
	cout << "the histos size = " << size << " and " << hj_q_h2->GetSize() << "\n";

	Int_t size_x = hj_q_h1->GetNbinsX();
	Int_t size_y = hj_q_h1->GetNbinsY();
	Int_t size_z = hj_q_h1->GetNbinsZ();
	cout << "the histos size = " << size_x << " , " << size_y << " , " << size_z << "\n";

	for (int iz=1; iz<=size_z; ++iz)
	{
	for (int ix=1; ix<=size_x; ++ix)
	{
	for (int iy=1; iy<=size_y; ++iy)
		{
		/*
		if ((i / size) > 0.5 && ! passed_midwork)
			{
			passed_midwork = 1;
			cout << "Passed midwork!\n";
			}
		*/
		//cout << "AAA\n";

		double xcenter1 = (double) hj_q_h1->GetXaxis()->GetBinCenter(ix);
		double xcenter2 = (double) hj_q_h2->GetXaxis()->GetBinCenter(ix);
		double ycenter1 = (double) hj_q_h1->GetYaxis()->GetBinCenter(iy);
		double ycenter2 = (double) hj_q_h2->GetYaxis()->GetBinCenter(iy);
		double zcenter1 = (double) hj_q_h1->GetZaxis()->GetBinCenter(iz);
		double zcenter2 = (double) hj_q_h2->GetZaxis()->GetBinCenter(iz);
		if (xcenter1 != xcenter2) return 1;

		// some brilliance from ROOT histograms (f*cking 3d array!)
		// Int_t bin = h->GetBin(binx,biny,binz);
		Int_t global_bin = hj_q_h1->GetBin(ix,iy,iz);
		Int_t global_bin2 = hj_q_h2->GetBin(ix,iy,iz);
		if (global_bin != global_bin2) return 2;
		// TODO: check the bins are the same in qcd jets?

		//cout << "BBB\n";

		double hj_w_n_jets     = hj_w_h1->GetBinContent( global_bin );
		double hj_w_n_tau_jets = hj_w_h2->GetBinContent( global_bin ); // no histo in j5.6 Data
		//double hj_w_n_tau_jets = 0;
		double hj_q_n_jets     = hj_q_h1->GetBinContent( global_bin );
		double hj_q_n_tau_jets = hj_q_h2->GetBinContent( global_bin );
		//double hj_q_n_tau_jets = 0;

		double hjm_w_n_jets     = hjm_w_h1->GetBinContent( global_bin );
		double hjm_w_n_tau_jets = hjm_w_h2->GetBinContent( global_bin );
		double hjm_q_n_jets     = hjm_q_h1->GetBinContent( global_bin );
		double hjm_q_n_tau_jets = hjm_q_h2->GetBinContent( global_bin );

		//cout << "CCC\n";

		fprintf(out_f, "%d, %g,%g,%g, %g,%g,%g,%g, %g,%g", global_bin, xcenter1, ycenter1, zcenter1, hj_w_n_jets,hj_w_n_tau_jets, hj_q_n_jets,hj_q_n_tau_jets, hjm_w_n_jets,hjm_w_n_tau_jets);
		fprintf(out_f, ",%g,%g", hjm_q_n_jets,hjm_q_n_tau_jets);
		fprintf(out_f, "\n");
		//fprintf(out_f, "%d,%g,%g\n", i, n_jets, n_tau_jets);
		//cout << "DDD\n";
		}
	}
	}

	file->Close();
	fclose(out_f);


	return 0;
}



// test
int foo();

// initializators
void set_bSF_calibrators();
void set_bSF_effs_for_dtag(TString dtag);

// method
double b_taggin_SF (double jet_pt, double jet_eta, double jet_b_discr, int jet_hadronFlavour, double b_tag_WP);

// method, trying initializing directly in global space
float met_pt_recoilcor(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets);

float MTlep_met_pt_recoilcor(float lep_px, float lep_py,
	float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets);

float zPtMass_weight(float genMass, float genPt);

float met_pt_recoilcor_x(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets);

float met_pt_recoilcor_y(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets);


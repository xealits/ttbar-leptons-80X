#ifndef RECORDFUNCS_H
#define RECORDFUNCS_H


#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

//#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"

// should work in CMSSW_8_0_12 and CMSSW_8_1_0
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
// #include "CondFormats/BTauObjects/interface/BTagCalibration.h"
// #include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
// this one is for 80X -> #include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h" 
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibrator.h"  
//#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibrator.h" 

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

#include <map>
#include <string>

using namespace std;


int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight);

int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight);

//int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight)
int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight);



/* Don't expose these
string mc_decay;

// channel -> {control_point, TH}
// 1 job = 1 dtag,
// 1 dtag may have several channels
std::map<string, std::map<string, TH1D>> th1d_distr_maps_control;
std::map<string, TH1D> th1d_distr_maps_control_headers;

std::map<string, std::map<string, TH1I>> th1i_distr_maps_control;
std::map<string, TH1I> th1i_distr_maps_control_headers;

std::map<string, std::map<string, TH2D>> th2d_distr_maps_control;
std::map<string, TH2D> th2d_distr_maps_control_headers;

std::map<string, std::map<string, TH3D>> th3d_distr_maps_control;
std::map<string, TH3D> th3d_distr_maps_control_headers;
*/

//extern int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight);
//extern int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight);


int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight);


int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight);


//int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight)
int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight);


int fill_3d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, double x, double y, double z, double weight);




// TODO: I wonder where Float_t is defined?
// but probably TH3F depends on it and pulls it in anyway

/* Don't expose these
// good bins 1, 2
Float_t largebins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 55, 75, 100, 150, 500 }; // 11 bins 12 edges
static int n_largebins_pt = 11;

// 0.1 pt bins
static Float_t bins_pt[49] = { 0,    20, 21, 22, 23, 24,    25, 26, 27, 28, 29,    30, 31, 32, 33, 34,    35, 36, 37, 38, 39,    40, 41, 42, 43, 44,    45, 46, 47, 48, 49, 50, 51, 52, 53, 54,    55, 56, 57, 58, 59, 60, 65, 70,    75, 80,    100,    150,    500 }; // 48 bins 49 edges
static int n_bins_pt = 48;

// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
Float_t largebins_eta[8] = { -3, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3 }; // 7 bins, 8 edges
int n_largebins_eta = 7;

// now (2016) we've got more statistics -- thus smaller eta bins
// 0.1 eta bins
static Float_t bins_eta[53] = { -3.0,    -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6,    -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6,    -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4,    0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4,    1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4,    2.5,    3.0 }; // 52 bins 53 edges
static int n_bins_eta = 52;


Float_t largebins_rad[14] = { 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4, 1, 2 }; // 13 bins 14 edges
int n_largebins_rad = 13;

// 0.1 radius bins:
static Float_t bins_rad[36] = { 0, 0.012,   0.025, 0.032, 0.04,   0.05, 0.06,   0.075, 0.082, 0.09,   0.1, 0.112,   0.125, 0.132, 0.14,   0.15, 0.16,   0.175, 0.182, 0.19,   0.2, 0.21, 0.22, 0.23, 0.24,    0.25, 0.26, 0.27, 0.28, 0.29,    0.3, 0.31, 0.32,    0.4, 1, 2 }; // 35 bins 36 edges
static int n_bins_rad = 35;


//float tau_fake_distance = 0.1; // first-try distance
static float tau_fake_distance = 0.3; // the distance to tau for a jet to be considered tau's origin
*/


double jet_radius(pat::Jet& jet);

double jet_radius(pat::Tau& jet);



int fill_jet_distr(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius);

int fill_jet_distr_large_bins(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius)


/* Also don't expose
// btag-efficiency bins
// 0.1 pt bins
static Float_t beff_bins_pt[19] = { 0, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 65, 70, 75, 80, 100, 150, 500 }; // 18 bins 19 edges
static int beff_n_bins_pt = 18;

// 0.1 eta bins
static Float_t beff_bins_eta[51] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 }; // 50 bins 51 edges
static int beff_n_bins_eta = 50;
*/


int fill_btag_efficiency(string control_point_name, Double_t pt, Double_t eta, Double_t weight);





int record_jets_fakerate_distrs(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, vector<LorentzVector>& visible_gen_taus, double event_weight, bool isMC);

int record_jets_fakerate_distrs_large_bins(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, vector<LorentzVector>& visible_gen_taus, double event_weight, bool isMC);



#endif /* RECORDFUNCS_H */


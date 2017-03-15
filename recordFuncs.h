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

//extern int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight);
//extern int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight);


int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th1d_distr_maps_control
	if (th1d_distr_maps_control.find(mc_decay) == th1d_distr_maps_control.end() )
		{
		//
		th1d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH1D>()));
		}

	if (th1d_distr_maps_control[mc_decay].find(control_point_name) == th1d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th1d_distr_maps_control[mc_decay][control_point_name] = TH1D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_maps_control[mc_decay][control_point_name].Fill(value, weight);

	if (th1d_distr_maps_control_headers.find(control_point_name) == th1d_distr_maps_control_headers.end() )
		{
		// th1d_distr_maps_control_headers[control_point_name] = TH1D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th1d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH1D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup)));
		}

	// return success:
	return 0;
	}




int fill_1i(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, int value, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th1i_distr_maps_control
	if (th1i_distr_maps_control.find(mc_decay) == th1i_distr_maps_control.end() )
		{
		//
		th1i_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH1I>()));
		}

	if (th1i_distr_maps_control[mc_decay].find(control_point_name) == th1i_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th1i_distr_maps_control[mc_decay][control_point_name] = TH1I((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1i_distr_maps_control[mc_decay][control_point_name].Fill(value, weight);

	if (th1i_distr_maps_control_headers.find(control_point_name) == th1i_distr_maps_control_headers.end() )
		{
		// th1i_distr_maps_control_headers[control_point_name] = TH1I("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th1i_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH1I((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup)));
		}

	// return success:
	return 0;
	}





//int fill_1d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, double value, double weight)
int fill_2d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, double x, double y, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th2d_distr_maps_control
	if (th2d_distr_maps_control.find(mc_decay) == th2d_distr_maps_control.end() )
		{
		//
		th2d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH2D>()));
		}

	if (th2d_distr_maps_control[mc_decay].find(control_point_name) == th2d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th2d_distr_maps_control[mc_decay][control_point_name] = TH2D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th2d_distr_maps_control[mc_decay][control_point_name].Fill(x, y, weight);

	if (th2d_distr_maps_control_headers.find(control_point_name) == th2d_distr_maps_control_headers.end() )
		{
		// th2d_distr_maps_control_headers[control_point_name] = TH2D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th2d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH2D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup)));
		}

	// return success:
	return 0;
	}


int fill_3d(string control_point_name, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, double x, double y, double z, double weight)
	{
	// check if control point has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th3d_distr_maps_control
	if (th3d_distr_maps_control.find(mc_decay) == th3d_distr_maps_control.end() )
		{
		//
		th3d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH3D>()));
		}

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup );
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(x, y, z, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup)));
		}

	// return success:
	return 0;
	}




// TODO: I wonder where Float_t is defined?
// but probably TH3F depends on it and pulls it in anyway

// good bins 1, 2
Float_t largebins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 55, 75, 100, 150, 500 }; // 11 bins 12 edges
//Float_t bins_pt[11] = { 0, 30, 33, 37, 40, 43, 45, 48, 56, 63, 500 }; // 10 bins, 11 edges
//static Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 80, 150, 500 }; // 11 bins, 12 edges
static int n_largebins_pt = 11;
// bins 3 (exactly from AN489)
//Float_t bins_pt[12] = { 0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 500 }; // 13 bins, 14 edges
//int n_bins_pt = 13;
// bins 4 (based on j6.5 distributions)
//static Float_t bins_pt[15] = { 0, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 150, 500 }; // 14 bins 15 edges
//static int n_bins_pt = 14;

// 0.1 pt bins
static Float_t bins_pt[49] = { 0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 65, 70, 75, 80, 100, 150, 500 }; // 48 bins 49 edges
static int n_bins_pt = 48;




// Float_t bins_eta[6] = { -3, -1.5, -0.45, 0.45, 1.5, 3 }; // 5 bins, 6 edges
Float_t largebins_eta[8] = { -3, -2.5, -1.5, -0.45, 0.45, 1.5, 2.5, 3 }; // 7 bins, 8 edges
int n_largebins_eta = 7;
//static Float_t bins_eta[10] = { -3, -2.4, -2.3, -1.5, -0.45, 0.45, 1.5, 2.3, 2.4, 3 }; // 9 bins 10 edges
//static int n_bins_eta = 9;
// AN bins go in 0.5, and are done for abs eta
// I will use 0.5, but for both signs of eta:
//static Float_t bins_eta[15] = { -3, -2.4, -2.3, -2., -1.5, -1, -0.5, 0, 0.5, 1., 1.5, 2., 2.3, 2.4, 3 }; // 14 bins 15 edges
//static int n_bins_eta = 14;

// now (2016) we've got more statistics -- thus smaller eta bins
// 0.1 eta bins
static Float_t bins_eta[51] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 }; // 50 bins 51 edges
static int n_bins_eta = 50;


Float_t largebins_rad[16] = { 0, 0.06, 0.07, 0.08, 0.087, 0.093, 0.1, 0.107, 0.113, 0.12,
	0.127, 0.133, 0.14, 0.15, 0.16, 2 }; // 15 bins, 16 edges
int n_largebins_rad = 15;
// exactly from AN489:
//static Float_t bins_rad[12] = { 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.4, 1, 2 }; // 11 bins 12 edges
//static int n_bins_rad = 11;

// 0.1 radius bins:
static Float_t bins_rad[36] = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.4, 1, 2 }; // 35 bins 36 edges
static int n_bins_rad = 35;


//float tau_fake_distance = 0.1; // first-try distance
static float tau_fake_distance = 0.3; // the distance to tau for a jet to be considered tau's origin


double jet_radius(pat::Jet& jet)
	{
	//return sqrt(jet.EtaPhiMoments::etaEtaMoment + jet.EtaPhiMoments::phiPhiMoment);
	//return sqrt(jet.etaEtaMoment() + jet.phiPhiMoment());
	return sqrt(jet.etaetaMoment() + jet.phiphiMoment());
	}

double jet_radius(pat::Tau& jet)
	{
	//return sqrt(jet.EtaPhiMoments::etaEtaMoment + jet.EtaPhiMoments::phiPhiMoment);
	//return sqrt(jet.etaEtaMoment() + jet.phiPhiMoment());
	return sqrt(jet.etaetaMoment() + jet.phiphiMoment());
	}



int fill_jet_distr(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th3d_distr_maps_control
	if (th3d_distr_maps_control.find(mc_decay) == th3d_distr_maps_control.end() )
		{
		//
		th3d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH3D>()));
		}

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad);
		// th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}


	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, radius, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta, n_bins_rad, bins_rad)));
		}

	// return success:
	return 0;
	}


int fill_jet_distr_large_bins(string control_point_name, Double_t weight, Double_t pt, Double_t eta, Double_t radius)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th3d_distr_maps_control
	if (th3d_distr_maps_control.find(mc_decay) == th3d_distr_maps_control.end() )
		{
		//
		th3d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH3D>()));
		}

	if (th3d_distr_maps_control[mc_decay].find(control_point_name) == th3d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th3d_distr_maps_control[mc_decay][control_point_name] = TH3D((control_point_name + mc_decay).c_str(), ";;", n_largebins_pt, largebins_pt, n_largebins_eta, largebins_eta, n_largebins_rad, largebins_rad);
		// th3f_distr_control.insert( std::make_pair(key, TH3F(control_point_name.c_str(),      ";;", 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad)));
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}


	// fill the distribution:
	th3d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, radius, weight);

	if (th3d_distr_maps_control_headers.find(control_point_name) == th3d_distr_maps_control_headers.end() )
		{
		// th3d_distr_maps_control_headers[control_point_name] = TH3D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		// th3f_distr_control_headers.insert( std::make_pair(string("j_distr"), TH3F("Header of jets distribution",      ";;", 10, bins_pt, 5, bins_eta, 15, bins_rad)));
		th3d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH3D((string("Header of ") + control_point_name).c_str(), ";;", n_largebins_pt, largebins_pt, n_largebins_eta, largebins_eta, n_largebins_rad, largebins_rad)));
		}

	// return success:
	return 0;
	}


// btag-efficiency bins
// 0.1 pt bins
static Float_t beff_bins_pt[19] = { 0, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 65, 70, 75, 80, 100, 150, 500 }; // 18 bins 19 edges
static int beff_n_bins_pt = 18;

// 0.1 eta bins
static Float_t beff_bins_eta[51] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 }; // 50 bins 51 edges
static int beff_n_bins_eta = 50;


int fill_btag_efficiency(string control_point_name, Double_t pt, Double_t eta, Double_t weight)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	// std::pair <string,string> key (mc_decay, control_point_name);

	// create channel map in th2d_distr_maps_control
	if (th2d_distr_maps_control.find(mc_decay) == th2d_distr_maps_control.end() )
		{
		//
		th2d_distr_maps_control.insert( std::make_pair(mc_decay, std::map<string, TH2D>()));
		}

	if (th2d_distr_maps_control[mc_decay].find(control_point_name) == th2d_distr_maps_control[mc_decay].end() )
		{
		// the control point distr has not been created/initialized
		// THE HISTOGRAM NAMES HAVE TO BE DIFFERENT
		// BECAUSE O M G ROOT USES IT AS A REAL POINTER TO THE OBJECT
		//         O M G
		th2d_distr_maps_control[mc_decay][control_point_name] = TH2D((control_point_name + mc_decay).c_str(), ";;", beff_n_bins_pt, beff_bins_pt, beff_n_bins_eta, beff_bins_eta);
		// later on, when writing to the file,
		// I'll have to rename histograms on each write
		// and probably delete them along the way, so that they don't collide...
		// ROOT SUCKS
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th2d_distr_maps_control[mc_decay][control_point_name].Fill(pt, eta, weight);

	if (th2d_distr_maps_control_headers.find(control_point_name) == th2d_distr_maps_control_headers.end() )
		{
		th2d_distr_maps_control_headers.insert( std::make_pair(control_point_name, TH2D((string("Header of ") + control_point_name).c_str(), ";;", n_bins_pt, bins_pt, n_bins_eta, bins_eta)));
		}

	// return success:
	return 0;
	}





int record_jets_fakerate_distrs(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, vector<LorentzVector>& visible_gen_taus, double event_weight, bool isMC)
	{

	// for control record distr of given taus (not the matching jets -- just taus)
	for (size_t i = 0; i < selTaus.size(); ++i)
		{
		pat::Tau & tau = selTaus[i];
		// Calo Taus don't have radius?
		//fill_jet_distr(channel + selection + ("_taus_distr"), event_weight, tau.pt(), tau.eta(), jet_radius(tau));
		fill_2d(channel + selection + ("_taus_distr"), 100, 0., 300., 50, -2.5, 2.5, tau.pt(), tau.eta(), event_weight);
		}

	for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
		{
		pat::Jet& jet = selJets[ijet];
		// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
		// selection jets
		fill_jet_distr(channel + selection + ("_jets_distr"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
		//fill_3d(channel + selection + ("_jets_distr"), 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad, 300, 0, 300,   20, event_weight);

		// study jet_R ~ jet_eta -- does calculation of R depend on eta of jet?
		// (maybe it would be nice to use bins of the 3d jet distr?)
		fill_2d(channel + selection + ("_jet_radius_vs_eta"), 100, 0., 2., 120, -3., 3., jet_radius(jet), jet.eta(), event_weight);

		int jet_origin = -1; // for MC only

		// const reco::GenParticle* genParton()
		// jet parton origin
		if (isMC)
			{
			//const reco::GenParticle* jet_origin = jet.genParton();
			// the ID should be in:
			// jet_origin->pdgId();
			jet_origin = abs(jet.partonFlavour());

			// match to visible_gen_taus
			// if partonFlavour == 0 and matches to a tau (dR < tau_fake_distance (~= 0.3))
			// add it as 15 --- tau

			double minDRtj (9999.);
			//unsigned int closest_tau = -1; // not used now
			for (size_t i = 0; i < visible_gen_taus.size(); i++)
				{
				double jet_tau_distance = TMath::Min(minDRtj, reco::deltaR (jet, visible_gen_taus[i]));
				if (jet_tau_distance<minDRtj)
					{
					//closest_tau = i;
					minDRtj = jet_tau_distance;
					}
				}
			if (jet_origin == 0 && (minDRtj < tau_fake_distance)) jet_origin = 15; // BEST so far
			// maybe the true taus are lost here?
			// 1) dR can be too wide and include non-tau jets
			//    (where is the tau jet then? -- maybe mixed into a wide jet, or not reconstructed altogether?) 
			// 2) what if partonFlavour of a tau jet is not 0?
			//    maybe 0 are only wide tau jets, which have less ID efficiency,
			//    average tau jets usually get 1-4 light-quark partonFlavour etc
			// let's try 2)
			//if (minDRtj < tau_fake_distance) jet_origin = 15;
			// the fake rate in qcd is lower than it was and lower than in wjets
			// try partonFlavour != gluons
			//if (jet_origin != 21 && (minDRtj < tau_fake_distance)) jet_origin = 15;
			// 3) TODO: can also try light tau & other + light quarks (partonFlavour < 5)
			//if (jet_origin < 5 && (minDRtj < tau_fake_distance)) jet_origin = 15;

			fill_1d(channel + selection + ("_jet_partonFlavour"), 100, 0, 100,   jet_origin, event_weight);
			fill_1d(channel + selection + ("_jet_hadronFlavour"), 100, 0, 100,   abs(jet.hadronFlavour()), event_weight);

			// per-jet-origin distrs
			if (jet_origin == 21) // gluons
				fill_jet_distr(channel + selection + ("_jets_distr_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 5) // b-quarks
				fill_jet_distr(channel + selection + ("_jets_distr_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 0) // other stuff (taus probably)
				fill_jet_distr(channel + selection + ("_jets_distr_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 15) // the matched to gen_taus (dR < 1) jets
				fill_jet_distr(channel + selection + ("_jets_distr_t"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else // other (light) quarks
				fill_jet_distr(channel + selection + ("_jets_distr_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			}

		for(size_t itau=0; itau < selTaus.size(); ++itau)
			{
			// selection taus (fake taus)
			//for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			//{
			//double fake_distance = reco::deltaR(selJets[ijet], selTaus[itau]);
			double fake_distance = reco::deltaR(jet, selTaus[itau]);
			//qcd_taujet_distance->Fill(fake_distance);
			fill_1d(channel + selection + ("_taujet_distance"), 100, 0, 2,   fake_distance, event_weight);

			if (fake_distance <= tau_fake_distance)
				{
				// the tau is fake by this jet -- save distr
				//qcd_tau_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
				fill_jet_distr(channel + selection + ("_tau_jets_distr"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
				// fill_pt_pt(string control_point_name, double pt1, double pt2, double event_weight)
				//fill_pt_pt(channel + selection + ("_tau_taujet_pts"), jet.pt(), selTaus[itau].pt(), event_weight);
				fill_2d(channel + selection + ("_tau_taujet_pts"), 400, 0., 400., 400, 0., 400, jet.pt(), selTaus[itau].pt(), event_weight);

				// N tau-jets
				//increment( channel + selection + ("_selection_ntaujet"), event_weight );

				// const reco::GenParticle* genParton()
				// jet parton origin for faking jet (just in case)
				if (isMC)
					{
					//int jet_origin = abs(jet.partonFlavour());
					//const reco::GenParticle* jet_origin = selJets[ijet].genParton();
					// the ID should be in:
					// jet_origin->pdgId();
					//qcd_taujet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(channel + selection + ("_taujet_origins"), 100, 0, 100,   jet_origin, event_weight);
					//qcd_taujet_origin->Fill(abs( jet_origin->pdgId() ));
					// qcd_jet_origin

					// per-jet-origin distrs
					if (jet_origin == 21) // gluons
						fill_jet_distr(channel + selection + ("_tau_jets_distr_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 5) // b-quarks
						fill_jet_distr(channel + selection + ("_tau_jets_distr_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 0) // other stuff (taus probably)
						fill_jet_distr(channel + selection + ("_tau_jets_distr_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 15) // the matched to gen_taus (dR < 1) jets
						fill_jet_distr(channel + selection + ("_tau_jets_distr_t"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else // other (light) quarks
						fill_jet_distr(channel + selection + ("_tau_jets_distr_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					}
				continue;
				}
			//}
			}
		}

	// return success
	return 0;
	}

int record_jets_fakerate_distrs_large_bins(string channel, string selection, pat::JetCollection & selJets, pat::TauCollection & selTaus, vector<LorentzVector>& visible_gen_taus, double event_weight, bool isMC)
	{

	for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
		{
		pat::Jet& jet = selJets[ijet];
		// qcd_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
		// selection jets
		fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
		//fill_3d(channel + selection + ("_jets_distr"), 10, bins_pt, n_bins_eta, bins_eta, 15, bins_rad, 300, 0, 300,   20, event_weight);

		// study jet_R ~ jet_eta -- does calculation of R depend on eta of jet?
		// (maybe it would be nice to use bins of the 3d jet distr?)
		fill_2d(channel + selection + ("_jet_radius_vs_eta"), 100, 0., 2., 120, -3., 3., jet_radius(jet), jet.eta(), event_weight);

		int jet_origin = -1; // for MC only

		// const reco::GenParticle* genParton()
		// jet parton origin
		if (isMC)
			{
			//const reco::GenParticle* jet_origin = jet.genParton();
			// the ID should be in:
			// jet_origin->pdgId();
			jet_origin = abs(jet.partonFlavour());

			// match to visible_gen_taus
			// if partonFlavour == 0 and matches to a tau (dR < tau_fake_distance (~= 0.3))
			// add it as 15 --- tau

			double minDRtj (9999.);
			//unsigned int closest_tau = -1; // not used now
			for (size_t i = 0; i < visible_gen_taus.size(); i++)
				{
				double jet_tau_distance = TMath::Min(minDRtj, reco::deltaR (jet, visible_gen_taus[i]));
				if (jet_tau_distance<minDRtj)
					{
					//closest_tau = i;
					minDRtj = jet_tau_distance;
					}
				}
			if (jet_origin == 0 && (minDRtj < tau_fake_distance)) jet_origin = 15; // BEST so far
			// maybe the true taus are lost here?
			// 1) dR can be too wide and include non-tau jets
			//    (where is the tau jet then? -- maybe mixed into a wide jet, or not reconstructed altogether?) 
			// 2) what if partonFlavour of a tau jet is not 0?
			//    maybe 0 are only wide tau jets, which have less ID efficiency,
			//    average tau jets usually get 1-4 light-quark partonFlavour etc
			// let's try 2)
			//if (minDRtj < tau_fake_distance) jet_origin = 15;
			// the fake rate in qcd is lower than it was and lower than in wjets
			// try partonFlavour != gluons
			//if (jet_origin != 21 && (minDRtj < tau_fake_distance)) jet_origin = 15;
			// 3) TODO: can also try light tau & other + light quarks (partonFlavour < 5)
			//if (jet_origin < 5 && (minDRtj < tau_fake_distance)) jet_origin = 15;

			fill_1d(channel + selection + ("_jet_partonFlavour"), 100, 0, 100,   jet_origin, event_weight);
			fill_1d(channel + selection + ("_jet_hadronFlavour"), 100, 0, 100,   abs(jet.hadronFlavour()), event_weight);

			// per-jet-origin distrs
			if (jet_origin == 21) // gluons
				fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 5) // b-quarks
				fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 0) // other stuff (taus probably)
				fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else if (jet_origin == 15) // the matched to gen_taus (dR < 1) jets
				fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins_t"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			else // other (light) quarks
				fill_jet_distr_large_bins(channel + selection + ("_jets_distr_large_bins_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
			}

		for(size_t itau=0; itau < selTaus.size(); ++itau)
			{
			// selection taus (fake taus)
			//for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			//{
			//double fake_distance = reco::deltaR(selJets[ijet], selTaus[itau]);
			double fake_distance = reco::deltaR(jet, selTaus[itau]);
			//qcd_taujet_distance->Fill(fake_distance);
			fill_1d(channel + selection + ("_taujet_distance"), 100, 0, 2,   fake_distance, event_weight);

			if (fake_distance <= tau_fake_distance)
				{
				// the tau is fake by this jet -- save distr
				//qcd_tau_jets_distr->Fill(jet.pt(), jet.eta(), jet_radius(jet));
				fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
				// fill_pt_pt(string control_point_name, double pt1, double pt2, double event_weight)
				//fill_pt_pt(channel + selection + ("_tau_taujet_pts"), jet.pt(), selTaus[itau].pt(), event_weight);
				fill_2d(channel + selection + ("_tau_taujet_pts"), 400, 0., 400., 400, 0., 400, jet.pt(), selTaus[itau].pt(), event_weight);

				// N tau-jets
				//increment( channel + selection + ("_selection_ntaujet"), event_weight );

				// const reco::GenParticle* genParton()
				// jet parton origin for faking jet (just in case)
				if (isMC)
					{
					//int jet_origin = abs(jet.partonFlavour());
					//const reco::GenParticle* jet_origin = selJets[ijet].genParton();
					// the ID should be in:
					// jet_origin->pdgId();
					//qcd_taujet_origin->Fill(abs( jet.partonFlavour() ));
					fill_1d(channel + selection + ("_taujet_origins"), 100, 0, 100,   jet_origin, event_weight);
					//qcd_taujet_origin->Fill(abs( jet_origin->pdgId() ));
					// qcd_jet_origin

					// per-jet-origin distrs
					if (jet_origin == 21) // gluons
						fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins_g"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 5) // b-quarks
						fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins_b"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 0) // other stuff (taus probably)
						fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins_o"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else if (jet_origin == 15) // the matched to gen_taus (dR < 1) jets
						fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins_t"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					else // other (light) quarks
						fill_jet_distr_large_bins(channel + selection + ("_tau_jets_distr_large_bins_q"), event_weight, jet.pt(), jet.eta(), jet_radius(jet));
					}
				continue;
				}
			//}
			}
		}

	// return success
	return 0;
	}






#endif /* RECORDFUNCS_H */

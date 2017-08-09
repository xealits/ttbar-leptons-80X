#ifndef PROCESSINGELECTRONS_H
#define PROCESSINGELECTRONS_H

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/ttbar-leptons-80X/interface/recordFuncs.h"



int processElectrons_ID_ISO_Kinematics(pat::ElectronCollection& electrons, reco::Vertex goodPV, double rho, double weight, // input
	patUtils::llvvElecId::ElecId el_ID, patUtils::llvvElecId::ElecId veto_el_ID,                           // config/cuts
	patUtils::llvvElecIso::ElecIso el_ISO, patUtils::llvvElecIso::ElecIso veto_el_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	pat::ElectronCollection& selElectrons, LorentzVector& elDiff, unsigned int& nVetoE,                    // output
	bool record, bool debug); // more output

int processElectrons_MatchHLT(
	pat::ElectronCollection& electrons,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR,
	pat::ElectronCollection& electrons_matched
	);

bool processElectron_matchesHLT(
	pat::Electron& electron,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR
	);

#endif /* PROCESSINGELECTRONS_H */


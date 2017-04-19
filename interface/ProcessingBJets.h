#ifndef PROCESSINGBJETS_H
#define PROCESSINGBJETS_H

#include <string>

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TH2F.h"

#include "UserCode/ttbar-leptons-80X/interface/SystematicShifts.h"

using namespace std;

struct bTaggingEfficiencyHistograms {
	TH2F* b_alljet   ;
	TH2F* b_tagged   ;
	TH2F* c_alljet   ;
	TH2F* c_tagged   ;
	TH2F* udsg_alljet;
	TH2F* udsg_tagged;
	};

double bTagging_b_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta);

double bTagging_c_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta);

double bTagging_udsg_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta);


int processBJets_BTag(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
	BTagCalibrationReader& btagCal,
	struct bTaggingEfficiencyHistograms& bEffs,
	string& b_tagger_label, float b_tag_WP,
	pat::JetCollection& selBJets,                          // output
	bool record, bool debug); // more output

int processBJets_BTag_with_systematics(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
	BTagCalibrationReader& btagCal, // BTagSFUtil& btsfutil, old b-tag SF weighting, done with bEffs now
	struct bTaggingEfficiencyHistograms& bEffs,
	string& b_tagger_label, float b_tag_WP,
	//pat::JetCollection& selBJets,                          // output
	map<systematic_shift, pat::JetCollection>& selBJets,                          // output
	bool record, bool debug); // more output

#endif /* PROCESSINGBJETS_H */


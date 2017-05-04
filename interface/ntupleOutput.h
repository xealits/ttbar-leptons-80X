#ifndef NTUPLEOUTPUT_H
#define NTUPLEOUTPUT_H

/*
 * NTuple is:
 *   bunch of Float_t parameters
 *   and their names for the definition string
 *
 * need:
 *   obtain the def string from the ntuple structure
 *   the ntuple contains bunch of the same parameters for different particles (5 jets etc)
 *   initialize an empty ntuple (-1 init values)
 *   reset to an empty ntuple
 *   fill it with stuff
 *   Fill the ntuple with this output
 *
 * do separate maps (string -> Float_t and string -> N -> Float_t)
 * and then merge them
 */


#include <map>
#include <string>

#include "TSystem.h"
#include "TNtuple.h"


//using namespace std;


/* 
 * group maps as:
 *
 * aMCatNLO_weight
 * gen_t_pt
 * gen_tb_pt
 * NUP_gen
 * nvtx_gen
 * nvtx
 * nvtx_good
 * fixedGridRhoFastjetAll
 * fixedGridRhoFastjetCentral
 * fixedGridRhoFastjetCentralNeutral
 * fixedGridRhoFastjetCentralChargedPileUp
 * HLT_el
 * HLT_mu
 * leps_ID
 * nleps
 * njets
 * nbjets
 * ntaus
 * met_init
 * met_uncorrected
 * met_corrected
 * lj_peak_distance
 * lj_taumatched_peak_distance

 * tau_decay
 * tau_hasSecondaryVertex
 * tau_hcalEnergy
 * tau_hcalEnergyLeadChargedHadrCand

 * tau_secondaryVertexCov_00
 * tau_secondaryVertexCov_01
 * tau_secondaryVertexCov_02
 * tau_secondaryVertexCov_10
 * tau_secondaryVertexCov_11
 * tau_secondaryVertexCov_12
 * tau_secondaryVertexCov_20
 * tau_secondaryVertexCov_21
 * tau_secondaryVertexCov_22

 * lep1_id
 * lep1_eta
 * lep1_phi
 * lep1_pt
 * lep1_p
 * lep2_id
 * lep2_eta
 * lep2_phi
 * lep2_pt
 * lep2_p

 * jet1_id
 * jet1_eta
 * jet1_phi
 * jet1_pt
 * jet1_p
 * jet1_rad
 * jet1_b_discr
 * jet1_hadronFlavour
 * jet1_partonFlavour
 * jet2_id
 * jet2_eta
 * jet2_phi
 * jet2_pt
 * jet2_p
 * jet2_rad
 * jet2_b_discr
 * jet2_hadronFlavour
 * jet2_partonFlavour
 * jet3_id
 * jet3_eta
 * jet3_phi
 * jet3_pt
 * jet3_p
 * jet3_rad
 * jet3_b_discr
 * jet3_hadronFlavour
 * jet3_partonFlavour
 * jet4_id
 * jet4_eta
 * jet4_phi
 * jet4_pt
 * jet4_p
 * jet4_rad
 * jet4_b_discr
 * jet4_hadronFlavour
 * jet4_partonFlavour
 * jet5_id
 * jet5_eta
 * jet5_phi
 * jet5_pt
 * jet5_p
 * jet5_rad
 * jet5_b_discr
 * jet5_hadronFlavour
 * jet5_partonFlavour

 * tau1_id
 * tau1_eta
 * tau1_phi
 * tau1_pt
 * tau1_p
 * tau1_IDlev
 * tau2_id
 * tau2_eta
 * tau2_phi
 * tau2_pt
 * tau2_p
 * tau2_IDlev
 */


/*
 * these maps are ntuple output
 * they in principle should be stored in a struct or class
 * -- leaving this for later
 */
std::map<TString, Float_t> NT_common_event = {
	{"aMCatNLO_weight", -1},
	{"gen_t_pt", -1},
	{"gen_tb_pt", -1},
	{"NUP_gen", -1},
	{"nvtx_gen", -1},
	{"nvtx", -1},
	{"nvtx_good", -1},
	{"fixedGridRhoFastjetAll", -1},
	{"fixedGridRhoFastjetCentral", -1},
	{"fixedGridRhoFastjetCentralNeutral", -1},
	{"fixedGridRhoFastjetCentralChargedPileUp", -1},
	{"HLT_el", -1},
	{"HLT_mu", -1},
	{"leps_ID", -1},
	{"nleps", -1},
	{"njets", -1},
	{"nbjets", -1},
	{"ntaus", -1},
	{"met_init", -1},
	{"met_uncorrected", -1},
	{"met_corrected", -1},
	{"lj_peak_distance", -1},
	{"lj_taumatched_peak_distance", -1},
	{"gen_t_w_decay_id", -1}, // = id of lepton (but the sign means which product is lepton: minus=1, plus=2) or 1 for quarks
	{"gen_t_w_p1_eta",  -1},
	{"gen_t_w_p1_phi",  -1},
	{"gen_t_w_p2_eta",  -1},
	{"gen_t_w_p2_phi",  -1},
	{"gen_tb_w_decay_id", -1}, // = id of lepton (same sign thing) or 1 for quarks
	{"gen_tb_w_p1_eta", -1},
	{"gen_tb_w_p1_phi", -1},
	{"gen_tb_w_p2_eta", -1},
	{"gen_tb_w_p2_phi", -1} };


/*
 * These are on tau decay
 * More in class reference:
 * http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_8_0_25/doc/html/d1/de9/classpat_1_1Tau.html
 *
 * the EnergyLead/all_Energy should provide something similar to R variable
 *   -- a (not so great) discriminant for true/fake taus, which is also correlated with tau polarization (to check?)
 * the secondary vertex and decay should provide selection of the 3-pion decay of hadronic taus
 * there are more available, tracks etc
 */
std::map<TString, Float_t> NT_leading_tau_decay_info = {
	{"tau_decay", -1},
	{"tau_hasSecondaryVertex", -1},
	{"tau_hcalEnergy", -1},
	{"tau_hcalEnergyLeadChargedHadrCand", -1}};

Float_t NT_leading_tau_secondaryVertexCov[3][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};


std::map<TString, vector<Float_t>> NT_leptons = {
	{"lep_id_", {-1, -1}},
	{"lep_eta_", {-1, -1}},
	{"lep_phi_", {-1, -1}},
	{"lep_pt_", {-1, -1}},
	{"lep_p_", {-1, -1}}};

std::map<TString, vector<Float_t>> NT_jets = {
	{"jet_id_", {-1, -1, -1, -1, -1}},
	{"jet_eta_", {-1, -1, -1, -1, -1}},
	{"jet_phi_", {-1, -1, -1, -1, -1}},
	{"jet_pt_", {-1, -1, -1, -1, -1}},
	{"jet_p_", {-1, -1, -1, -1, -1}},
	{"jet_rad_", {-1, -1, -1, -1, -1}},
	{"jet_b_discr_", {-1, -1, -1, -1, -1}},
	{"jet_hadronFlavour_", {-1, -1, -1, -1, -1}},
	{"jet_partonFlavour_", {-1, -1, -1, -1, -1}}};


std::map<TString, vector<Float_t>> NT_taus = {
	{"tau_id_", {-1, -1}},
	{"tau_eta_", {-1, -1}},
	{"tau_phi_", {-1, -1}},
	{"tau_pt_", {-1, -1}},
	{"tau_p_", {-1, -1}},
	{"tau_IDlev_", {-1, -1}}};
// tau IDlev is the highest ID of the tau:
// 3 -- Tight
// 2 -- Medium
// 1 -- loose
// 0 -- it somehow passed into the tau collection but doesn't pass any of the 3 IDs

#define N_OUTPUT 113+1 // placeholder for nevent
// for now: define N ouf output values (N of -1 in the above)



/*
 * construct the name and return TString
 * 
 * access with:
 * 	const char* TString::Data() const
 */
TString return_ntuple_names()
{
TString s("nevent");

for(std::map<TString, Float_t>::iterator it = NT_common_event.begin(); it != NT_common_event.end(); ++it)
	{
	TString controlpoint_name = it->first;
	s += ":" + controlpoint_name;
	}

for(std::map<TString, Float_t>::iterator it = NT_leading_tau_decay_info.begin(); it != NT_leading_tau_decay_info.end(); ++it)
	{
	TString controlpoint_name = it->first;
	s += ":" + controlpoint_name;
	}

for (Long_t r=0; r<3; r++)
	for (Long_t c=0; c<3; c++)
		{
		s += ":" + TString("tau_secondaryVertexCov_") + r + c;
		}

for(std::map<TString, vector<Float_t>>::iterator it = NT_leptons.begin(); it != NT_leptons.end(); ++it)
	{
	TString controlpoint_name = it->first;
	for (Long_t i=0; i<it->second.size(); i++)
		{
		s += ":" + controlpoint_name + i;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_jets.begin(); it != NT_jets.end(); ++it)
	{
	TString controlpoint_name = it->first;
	for (Long_t i=0; i<it->second.size(); i++)
		{
		s += ":" + controlpoint_name + i;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_taus.begin(); it != NT_taus.end(); ++it)
	{
	TString controlpoint_name = it->first;
	for (Long_t i=0; i<it->second.size(); i++)
		{
		s += ":" + controlpoint_name + i;
		}
	}

return s;
}




void fill_ntuple(TNtuple& ntuple, bool debug)
{
// main problem: construct a correct array const Float_t * to feed into ntuple.Fill();
// need N of floats...
Float_t output[N_OUTPUT];
int i = 0;

if (debug) cout << "ntuple has " << ntuple.GetNbranches() << " branches, filling it with " << N_OUTPUT << " array, got ";

output[i] = 0; // the first element is placeholder for Nlumi/Nevent stuff
i++;

for(std::map<TString, Float_t>::iterator it = NT_common_event.begin(); it != NT_common_event.end(); ++it)
	{
	output[i] = it->second;
	i++;
	}

for(std::map<TString, Float_t>::iterator it = NT_leading_tau_decay_info.begin(); it != NT_leading_tau_decay_info.end(); ++it)
	{
	output[i] = it->second;
	i++;
	}

for (Long_t r=0; r<3; r++)
	for (Long_t c=0; c<3; c++)
		{
		//s += ":" + TString("tau_secondaryVertexCov_") + r + c;
		output[i] = NT_leading_tau_secondaryVertexCov[r][c];
		i++;
		}

for(std::map<TString, vector<Float_t>>::iterator it = NT_leptons.begin(); it != NT_leptons.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		output[i] = v[j];
		i++;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_jets.begin(); it != NT_jets.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		output[i] = v[j];
		i++;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_taus.begin(); it != NT_taus.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		output[i] = v[j];
		i++;
		}
	}

if (debug) cout << i << " records" << endl;

ntuple.Fill(output);
}




void reset_ntuple_output()
{
// loop and reset all to -1
Float_t output[N_OUTPUT];
//int i = 0;

for(std::map<TString, Float_t>::iterator it = NT_common_event.begin(); it != NT_common_event.end(); ++it)
	{
	NT_common_event[it->first] = -1;
	}

for(std::map<TString, Float_t>::iterator it = NT_leading_tau_decay_info.begin(); it != NT_leading_tau_decay_info.end(); ++it)
	{
	NT_leading_tau_decay_info[it->first] = -1;
	}

for (Long_t r=0; r<3; r++)
	for (Long_t c=0; c<3; c++)
		{
		NT_leading_tau_secondaryVertexCov[r][c] = -1;
		}

for(std::map<TString, vector<Float_t>>::iterator it = NT_leptons.begin(); it != NT_leptons.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		v[j] = -1;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_jets.begin(); it != NT_jets.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		v[j] = -1;
		}
	}

for(std::map<TString, vector<Float_t>>::iterator it = NT_taus.begin(); it != NT_taus.end(); ++it)
	{
	vector<Float_t>& v = it->second;
	for (int j=0; j<v.size(); j++)
		{
		v[j] = -1;
		}
	}
}



#endif /* NTUPLEOUTPUT_H */


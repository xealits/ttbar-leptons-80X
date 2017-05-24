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


//#include <map>
//#include <string>
//
//#include "TSystem.h"
//#include "TNtuple.h"

/* macro declaring the object and setting a branch with its' pointer --- all in current, __not_global__ space (in main space)
 *
 * Notice the protocol:
 *    1) the object name in current namespace is `NT_Name`
 *    2) the branch name in the ntuple is `Name`
 */
#define OBJECT_in_NTuple(NTuple, Class, Name) Class NT_##Name; NTuple.Branch(#Name, #Class, &NT_##Name)

//using namespace std;

/*
 * So, this file ties the interface to our ntuple in the current namespace
 *
 * Usage
 *
 * declare your ntuple:
 *
 *     TNtuple output;
 *     // if the name is not `ntuple` (which is assumed here in the interface)
 *     // define your name for preprocessor:
 *     #define OUTNTUPLE output
 * 
 *     // load this interface:
 *     #include "ntupleOutput.h"
 *
 * now you have NT_Name objects in the name space and the ntuple has branches "Name" with pointers to these objects
 * copy objects from event (pseudocode):
 *
 *     NT_foo = events[i]["foo"])
 * 
 * when done fill the ntuple:
 *
 *     output.Fill();
 *
 * clearing/reseting of the objects for each event -- currently responsibility of the programmer, there is a sketchy macro for this now, it is in development
 */

// default name of the output ntuple
#ifndef OUTNTUPLE
	#define OUTNTUPLE ntuple
#endif

/* This works for only 1 ntuple.
 * But not for many in case of different decays of MC (tt).
 * A nice way is to add generation info to the ntuple output and split the decays then later, instead of having different files.
 * It's not large, 1 more int or whatever, and should work fine.
 * --- taking this
 */

// the interface (all floats, since it is really just TNtuple
OBJECT_in_NTuple(OUTNTUPLE, Float_t, aMCatNLO_weight);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, NUP_gen); // TODO: add gen info from TTbar
OBJECT_in_NTuple(OUTNTUPLE, Float_t, nvtx_gen);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, nvtx);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, nvtx_good);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, fixedGridRhoFastjetAll);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, fixedGridRhoFastjetCentral);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, fixedGridRhoFastjetCentralNeutral);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, fixedGridRhoFastjetCentralChargedPileUp);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, HLT_el);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, HLT_mu);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, leps_ID);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, nleps);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, njets);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, nbjets);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, ntaus);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, met_init);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, met_uncorrected);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, met_corrected);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lj_peak_distance);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lj_taumatched_peak_distance);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_decay);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_hasSecondaryVertex);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_hcalEnergy);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_hcalEnergyLeadChargedHadrCand);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_00);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_01);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_02);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_10);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_11);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_12);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_20);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_21);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau_secondaryVertexCov_22);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep2_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep2_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep2_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep2_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep2_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet1_partonFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet2_partonFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet3_partonFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet4_partonFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet5_partonFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_IDlev);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau2_IDlev);

// convenience:
#define NT_tau_secondaryVertexCov_(i, j) NT_tau_secondaryVertexCov_ ##i ##j

#define NT_LEPTONS_N 2
#define NT_lep_id_(i)     NT_lep_id_##i
#define NT_lep_eta_(i)    NT_lep_eta_##i
#define NT_lep_phi_(i)    NT_lep_phi_##i
#define NT_lep_pt_(i)     NT_lep_pt_##i
#define NT_lep_p_(i)      NT_lep_p_##i

#define NT_JETS_N 5
#define NT_jet_id_(i)                NT_jet_id_##i
#define NT_jet_eta_(i)               NT_jet_eta_##i
#define NT_jet_phi_(i)               NT_jet_phi_##i
#define NT_jet_pt_(i)                NT_jet_pt_##i
#define NT_jet_p_(i)                 NT_jet_p_##i
#define NT_jet_rad_(i)               NT_jet_rad_##i
#define NT_jet_b_discr_(i)           NT_jet_b_discr_##i
#define NT_jet_hadronFlavour_(i)     NT_jet_hadronFlavour_##i
#define NT_jet_partonFlavour_(i)     NT_jet_partonFlavour_##i

#define NT_TAUS_N 2
#define NT_tau_id_(i)     NT_tau_id_##i
#define NT_tau_eta_(i)    NT_tau_eta_##i
#define NT_tau_phi_(i)    NT_tau_phi_##i
#define NT_tau_pt_(i)     NT_tau_pt_##i
#define NT_tau_p_(i)      NT_tau_p_##i
#define NT_tau_IDlev_(i)  NT_tau_IDlev_##i



// the automatic reset of all parameters for now
#define RESET_NTUPLE_PARAMETERS \ \
NT_aMCatNLO_weight = -1; \
NT_gen_t_pt = -1; \
NT_gen_tb_pt = -1; \
NT_NUP_gen = -1; \
NT_nvtx_gen = -1; \
NT_nvtx = -1; \
NT_nvtx_good = -1; \
NT_fixedGridRhoFastjetAll = -1; \
NT_fixedGridRhoFastjetCentral = -1; \
NT_fixedGridRhoFastjetCentralNeutral = -1; \
NT_fixedGridRhoFastjetCentralChargedPileUp = -1; \
NT_HLT_el = -1; \
NT_HLT_mu = -1; \
NT_leps_ID = -1; \
NT_nleps = -1; \
NT_njets = -1; \
NT_nbjets = -1; \
NT_ntaus = -1; \
NT_met_init = -1; \
NT_met_uncorrected = -1; \
NT_met_corrected = -1; \
NT_lj_peak_distance = -1; \
NT_lj_taumatched_peak_distance = -1; \
NT_tau_decay = -1; \
NT_tau_hasSecondaryVertex = -1; \
NT_tau_hcalEnergy = -1; \
NT_tau_hcalEnergyLeadChargedHadrCand = -1; \
NT_tau_secondaryVertexCov_00 = -1; \
NT_tau_secondaryVertexCov_01 = -1; \
NT_tau_secondaryVertexCov_02 = -1; \
NT_tau_secondaryVertexCov_10 = -1; \
NT_tau_secondaryVertexCov_11 = -1; \
NT_tau_secondaryVertexCov_12 = -1; \
NT_tau_secondaryVertexCov_20 = -1; \
NT_tau_secondaryVertexCov_21 = -1; \
NT_tau_secondaryVertexCov_22 = -1; \
NT_lep1_id = -1; \
NT_lep1_eta = -1; \
NT_lep1_phi = -1; \
NT_lep1_pt = -1; \
NT_lep1_p = -1; \
NT_lep2_id = -1; \
NT_lep2_eta = -1; \
NT_lep2_phi = -1; \
NT_lep2_pt = -1; \
NT_lep2_p = -1; \
NT_jet1_id = -1; \
NT_jet1_eta = -1; \
NT_jet1_phi = -1; \
NT_jet1_pt = -1; \
NT_jet1_p = -1; \
NT_jet1_rad = -1; \
NT_jet1_b_discr = -1; \
NT_jet1_hadronFlavour = -1; \
NT_jet1_partonFlavour = -1; \
NT_jet2_id = -1; \
NT_jet2_eta = -1; \
NT_jet2_phi = -1; \
NT_jet2_pt = -1; \
NT_jet2_p = -1; \
NT_jet2_rad = -1; \
NT_jet2_b_discr = -1; \
NT_jet2_hadronFlavour = -1; \
NT_jet2_partonFlavour = -1; \
NT_jet3_id = -1; \
NT_jet3_eta = -1; \
NT_jet3_phi = -1; \
NT_jet3_pt = -1; \
NT_jet3_p = -1; \
NT_jet3_rad = -1; \
NT_jet3_b_discr = -1; \
NT_jet3_hadronFlavour = -1; \
NT_jet3_partonFlavour = -1; \
NT_jet4_id = -1; \
NT_jet4_eta = -1; \
NT_jet4_phi = -1; \
NT_jet4_pt = -1; \
NT_jet4_p = -1; \
NT_jet4_rad = -1; \
NT_jet4_b_discr = -1; \
NT_jet4_hadronFlavour = -1; \
NT_jet4_partonFlavour = -1; \
NT_jet5_id = -1; \
NT_jet5_eta = -1; \
NT_jet5_phi = -1; \
NT_jet5_pt = -1; \
NT_jet5_p = -1; \
NT_jet5_rad = -1; \
NT_jet5_b_discr = -1; \
NT_jet5_hadronFlavour = -1; \
NT_jet5_partonFlavour = -1; \
NT_tau1_id = -1; \
NT_tau1_eta = -1; \
NT_tau1_phi = -1; \
NT_tau1_pt = -1; \
NT_tau1_p = -1; \
NT_tau1_IDlev = -1; \
NT_tau2_id = -1; \
NT_tau2_eta = -1; \
NT_tau2_phi = -1; \
NT_tau2_pt = -1; \
NT_tau2_p = -1; \
NT_tau2_IDlev = -1; \

#endif /* NTUPLEOUTPUT_H */


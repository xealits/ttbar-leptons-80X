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
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_w_decay_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_w_p1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_w_p1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_w_p2_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_t_w_p2_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_w_decay_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_w_p1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_w_p1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_w_p2_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, gen_tb_w_p2_phi);

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
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep0_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep0_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep0_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep0_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep0_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, lep1_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_rad);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_b_discr);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_hadronFlavour);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, jet0_partonFlavour);
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
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau0_IDlev);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_id);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_eta);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_phi);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_pt);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_p);
OBJECT_in_NTuple(OUTNTUPLE, Float_t, tau1_IDlev);

// convenience:
#define NT_tau_secondaryVertexCov_(i, j) NT_tau_secondaryVertexCov_ ##i ##j

#define NT_LEPTONS_N 2
#define NT_JETS_N 5
#define NT_TAUS_N 2

// this should be static, no run-time i
// but it's in the loops now...
// thus making the case for now
#define NT_lep(i, id, eta, phi, pt, p) case i: { \
NT_lep ##i ##_id   = id;  \
NT_lep ##i ##_eta  = eta; \
NT_lep ##i ##_phi  = phi; \
NT_lep ##i ##_pt   = pt;  \
NT_lep ##i ##_p    = p;   \
break; }

#define NT_jet(i, id, eta, phi, pt, p, rad, b_discr, hadronFlavour, partonFlavour) case i: { \
NT_jet ##i ##_id                = id; \
NT_jet ##i ##_eta               = eta; \
NT_jet ##i ##_phi               = phi; \
NT_jet ##i ##_pt                = pt; \
NT_jet ##i ##_p                 = p; \
NT_jet ##i ##_rad               = rad; \
NT_jet ##i ##_b_discr           = b_discr; \
NT_jet ##i ##_hadronFlavour     = hadronFlavour; \
NT_jet ##i ##_partonFlavour     = partonFlavour; \
break; }

#define NT_tau(i, id, eta, phi, pt, p, IDlev) case i: { \
NT_tau ##i ##_id     = id; \
NT_tau ##i ##_eta    = eta; \
NT_tau ##i ##_phi    = phi; \
NT_tau ##i ##_pt     = pt; \
NT_tau ##i ##_p      = p; \
NT_tau ##i ##_IDlev  = IDlev; \
break; }



// the automatic reset of all parameters for now
#define RESET_NTUPLE_PARAMETERS \
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
NT_lep0_id = -1; \
NT_lep0_eta = -1; \
NT_lep0_phi = -1; \
NT_lep0_pt = -1; \
NT_lep0_p = -1; \
NT_lep1_id = -1; \
NT_lep1_eta = -1; \
NT_lep1_phi = -1; \
NT_lep1_pt = -1; \
NT_lep1_p = -1; \
NT_jet0_id = -1; \
NT_jet0_eta = -1; \
NT_jet0_phi = -1; \
NT_jet0_pt = -1; \
NT_jet0_p = -1; \
NT_jet0_rad = -1; \
NT_jet0_b_discr = -1; \
NT_jet0_hadronFlavour = -1; \
NT_jet0_partonFlavour = -1; \
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
NT_tau0_id = -1; \
NT_tau0_eta = -1; \
NT_tau0_phi = -1; \
NT_tau0_pt = -1; \
NT_tau0_p = -1; \
NT_tau0_IDlev = -1; \
NT_tau1_id = -1; \
NT_tau1_eta = -1; \
NT_tau1_phi = -1; \
NT_tau1_pt = -1; \
NT_tau1_p = -1; \
NT_tau1_IDlev = -1; \
NT_gen_t_w_decay_id = -1; \
NT_gen_t_w_p1_eta = -1; \
NT_gen_t_w_p1_phi = -1; \
NT_gen_t_w_p2_eta = -1; \
NT_gen_t_w_p2_phi = -1; \
NT_gen_tb_w_decay_id = -1; \
NT_gen_tb_w_p1_eta = -1; \
NT_gen_tb_w_p1_phi = -1; \
NT_gen_tb_w_p2_eta = -1; \
NT_gen_tb_w_p2_phi = -1;

#endif /* NTUPLEOUTPUT_H */


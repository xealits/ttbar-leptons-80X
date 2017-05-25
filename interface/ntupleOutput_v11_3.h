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


/* macro declaring the object and setting a branch with its' pointer --- all in current, __not_global__ space (in main space)
 *
 * Notice the protocol:
 *    1) the object name in current namespace is `NT_Name`
 *    2) the branch name in the ntuple is `Name`
 *
 * Also (TODO check) NTuple.Branch should open existing branch or create new one in empty ntuple.
 */
#if defined(NTUPLE_INTERFACE_CREATE)
	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.Branch(#Name, #Class, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         Float_t NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/F")
	#define Int_t_in_NTuple(NTuple, Name)           Int_t   NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/I")
#elif defined(NTUPLE_INTERFACE_OPEN)
	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.SetBranchAddress(#Name, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         OBJECT_in_NTuple(NTuple, Float_t, Name)
	#define Int_t_in_NTuple(NTuple, Name)           OBJECT_in_NTuple(NTuple, Int_t, Name)
#else
	error: set ntuple interface mode
#endif



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

// default name of the output
#ifndef OUTNTUPLE
	#define OUTNTUPLE NT_output_ttree
#endif

/* This works for only 1 ntuple.
 * But not for many in case of different decays of MC (tt).
 * A nice way is to add generation info to the ntuple output and split the decays then later, instead of having different files.
 * It's not large, 1 more int or whatever, and should work fine.
 * --- taking this
 */

// the interface (all Float_ts, since it is really just TNtuple
Float_t_in_NTuple(OUTNTUPLE, aMCatNLO_weight);
Float_t_in_NTuple(OUTNTUPLE, gen_t_pt);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_pt);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_decay_id); // = id of lepton (11/13/15, but the sign means which product is lepton: minus=1, plus=2) or 1 for quarks
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_decay_id);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_phi);
Float_t_in_NTuple(OUTNTUPLE, NUP_gen); // TODO: add gen info from TTbar
Float_t_in_NTuple(OUTNTUPLE, nvtx_gen);
Float_t_in_NTuple(OUTNTUPLE, nvtx);
Float_t_in_NTuple(OUTNTUPLE, nvtx_good);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetAll);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralNeutral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralChargedPileUp);
Float_t_in_NTuple(OUTNTUPLE, HLT_el);
Float_t_in_NTuple(OUTNTUPLE, HLT_mu);
Float_t_in_NTuple(OUTNTUPLE, leps_ID);
Float_t_in_NTuple(OUTNTUPLE, nleps);
Float_t_in_NTuple(OUTNTUPLE, njets);
Float_t_in_NTuple(OUTNTUPLE, nbjets);
Float_t_in_NTuple(OUTNTUPLE, ntaus);
Float_t_in_NTuple(OUTNTUPLE, met_init);
Float_t_in_NTuple(OUTNTUPLE, met_uncorrected);
Float_t_in_NTuple(OUTNTUPLE, met_corrected);
Float_t_in_NTuple(OUTNTUPLE, lj_peak_distance);
Float_t_in_NTuple(OUTNTUPLE, lj_taumatched_peak_distance);
Float_t_in_NTuple(OUTNTUPLE, tau_decay);
Float_t_in_NTuple(OUTNTUPLE, tau_hasSecondaryVertex);
Float_t_in_NTuple(OUTNTUPLE, tau_hcalEnergy);
Float_t_in_NTuple(OUTNTUPLE, tau_hcalEnergyLeadChargedHadrCand);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_00);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_01);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_02);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_10);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_11);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_12);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_20);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_21);
Float_t_in_NTuple(OUTNTUPLE, tau_secondaryVertexCov_22);
Float_t_in_NTuple(OUTNTUPLE, lep_id_0);
Float_t_in_NTuple(OUTNTUPLE, lep_eta_0);
Float_t_in_NTuple(OUTNTUPLE, lep_phi_0);
Float_t_in_NTuple(OUTNTUPLE, lep_pt_0);
Float_t_in_NTuple(OUTNTUPLE, lep_p_0);
Float_t_in_NTuple(OUTNTUPLE, lep_id_1);
Float_t_in_NTuple(OUTNTUPLE, lep_eta_1);
Float_t_in_NTuple(OUTNTUPLE, lep_phi_1);
Float_t_in_NTuple(OUTNTUPLE, lep_pt_1);
Float_t_in_NTuple(OUTNTUPLE, lep_p_1);
Float_t_in_NTuple(OUTNTUPLE, jet_id_0);
Float_t_in_NTuple(OUTNTUPLE, jet_eta_0);
Float_t_in_NTuple(OUTNTUPLE, jet_phi_0);
Float_t_in_NTuple(OUTNTUPLE, jet_pt_0);
Float_t_in_NTuple(OUTNTUPLE, jet_p_0);
Float_t_in_NTuple(OUTNTUPLE, jet_rad_0);
Float_t_in_NTuple(OUTNTUPLE, jet_b_discr_0);
Float_t_in_NTuple(OUTNTUPLE, jet_hadronFlavour_0);
Float_t_in_NTuple(OUTNTUPLE, jet_partonFlavour_0);
Float_t_in_NTuple(OUTNTUPLE, jet_id_1);
Float_t_in_NTuple(OUTNTUPLE, jet_eta_1);
Float_t_in_NTuple(OUTNTUPLE, jet_phi_1);
Float_t_in_NTuple(OUTNTUPLE, jet_pt_1);
Float_t_in_NTuple(OUTNTUPLE, jet_p_1);
Float_t_in_NTuple(OUTNTUPLE, jet_rad_1);
Float_t_in_NTuple(OUTNTUPLE, jet_b_discr_1);
Float_t_in_NTuple(OUTNTUPLE, jet_hadronFlavour_1);
Float_t_in_NTuple(OUTNTUPLE, jet_partonFlavour_1);
Float_t_in_NTuple(OUTNTUPLE, jet_id_2);
Float_t_in_NTuple(OUTNTUPLE, jet_eta_2);
Float_t_in_NTuple(OUTNTUPLE, jet_phi_2);
Float_t_in_NTuple(OUTNTUPLE, jet_pt_2);
Float_t_in_NTuple(OUTNTUPLE, jet_p_2);
Float_t_in_NTuple(OUTNTUPLE, jet_rad_2);
Float_t_in_NTuple(OUTNTUPLE, jet_b_discr_2);
Float_t_in_NTuple(OUTNTUPLE, jet_hadronFlavour_2);
Float_t_in_NTuple(OUTNTUPLE, jet_partonFlavour_2);
Float_t_in_NTuple(OUTNTUPLE, jet_id_3);
Float_t_in_NTuple(OUTNTUPLE, jet_eta_3);
Float_t_in_NTuple(OUTNTUPLE, jet_phi_3);
Float_t_in_NTuple(OUTNTUPLE, jet_pt_3);
Float_t_in_NTuple(OUTNTUPLE, jet_p_3);
Float_t_in_NTuple(OUTNTUPLE, jet_rad_3);
Float_t_in_NTuple(OUTNTUPLE, jet_b_discr_3);
Float_t_in_NTuple(OUTNTUPLE, jet_hadronFlavour_3);
Float_t_in_NTuple(OUTNTUPLE, jet_partonFlavour_3);
Float_t_in_NTuple(OUTNTUPLE, jet_id_4);
Float_t_in_NTuple(OUTNTUPLE, jet_eta_4);
Float_t_in_NTuple(OUTNTUPLE, jet_phi_4);
Float_t_in_NTuple(OUTNTUPLE, jet_pt_4);
Float_t_in_NTuple(OUTNTUPLE, jet_p_4);
Float_t_in_NTuple(OUTNTUPLE, jet_rad_4);
Float_t_in_NTuple(OUTNTUPLE, jet_b_discr_4);
Float_t_in_NTuple(OUTNTUPLE, jet_hadronFlavour_4);
Float_t_in_NTuple(OUTNTUPLE, jet_partonFlavour_4);
Float_t_in_NTuple(OUTNTUPLE, tau_id_0);
Float_t_in_NTuple(OUTNTUPLE, tau_eta_0);
Float_t_in_NTuple(OUTNTUPLE, tau_phi_0);
Float_t_in_NTuple(OUTNTUPLE, tau_pt_0);
Float_t_in_NTuple(OUTNTUPLE, tau_p_0);
Float_t_in_NTuple(OUTNTUPLE, tau_IDlev_0);
Float_t_in_NTuple(OUTNTUPLE, tau_id_1);
Float_t_in_NTuple(OUTNTUPLE, tau_eta_1);
Float_t_in_NTuple(OUTNTUPLE, tau_phi_1);
Float_t_in_NTuple(OUTNTUPLE, tau_pt_1);
Float_t_in_NTuple(OUTNTUPLE, tau_p_1);
Float_t_in_NTuple(OUTNTUPLE, tau_IDlev_1);

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
NT_lep_ID_s = -1; \
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
NT_tau_secondaryVertexCov_00 = -1; \
NT_tau_secondaryVertexCov_01 = -1; \
NT_tau_secondaryVertexCov_02 = -1; \
NT_tau_secondaryVertexCov_10 = -1; \
NT_tau_secondaryVertexCov_11 = -1; \
NT_tau_secondaryVertexCov_12 = -1; \
NT_lep_id_0 = -1; \
NT_lep_eta_0 = -1; \
NT_lep_phi_0 = -1; \
NT_lep_pt_0 = -1; \
NT_lep_p_0 = -1; \
NT_lep_id_1 = -1; \
NT_lep_eta_1 = -1; \
NT_lep_phi_1 = -1; \
NT_lep_pt_1 = -1; \
NT_lep_p_1 = -1; \
NT_jet_id_0 = -1; \
NT_jet_eta_0 = -1; \
NT_jet_phi_0 = -1; \
NT_jet_pt_0 = -1; \
NT_jet_p_0 = -1; \
NT_jet_rad_0 = -1; \
NT_jet_b_discr_0 = -1; \
NT_jet_hadronFlavour_0 = -1; \
NT_jet_partonFlavour_0 = -1; \
NT_jet_id_1 = -1; \
NT_jet_eta_1 = -1; \
NT_jet_phi_1 = -1; \
NT_jet_pt_1 = -1; \
NT_jet_p_1 = -1; \
NT_jet_rad_1 = -1; \
NT_jet_b_discr_1 = -1; \
NT_jet_hadronFlavour_1 = -1; \
NT_jet_partonFlavour_1 = -1; \
NT_jet_id_2 = -1; \
NT_jet_eta_2 = -1; \
NT_jet_phi_2 = -1; \
NT_jet_pt_2 = -1; \
NT_jet_p_2 = -1; \
NT_jet_rad_2 = -1; \
NT_jet_b_discr_2 = -1; \
NT_jet_hadronFlavour_2 = -1; \
NT_jet_partonFlavour_2 = -1; \
NT_jet_id_3 = -1; \
NT_jet_eta_3 = -1; \
NT_jet_phi_3 = -1; \
NT_jet_pt_3 = -1; \
NT_jet_p_3 = -1; \
NT_jet_rad_3 = -1; \
NT_jet_b_discr_3 = -1; \
NT_jet_hadronFlavour_3 = -1; \
NT_jet_partonFlavour_3 = -1; \
NT_jet_id_4 = -1; \
NT_jet_eta_4 = -1; \
NT_jet_phi_4 = -1; \
NT_jet_pt_4 = -1; \
NT_jet_p_4 = -1; \
NT_jet_rad_4 = -1; \
NT_jet_b_discr_4 = -1; \
NT_jet_hadronFlavour_4 = -1; \
NT_jet_partonFlavour_4 = -1; \
NT_tau_id_0 = -1; \
NT_tau_eta_0 = -1; \
NT_tau_phi_0 = -1; \
NT_tau_pt_0 = -1; \
NT_tau_p_0 = -1; \
NT_tau_IDlev_0 = -1; \
NT_tau_id_1 = -1; \
NT_tau_eta_1 = -1; \
NT_tau_phi_1 = -1; \
NT_tau_pt_1 = -1; \
NT_tau_p_1 = -1; \
NT_tau_IDlev_1 = -1; \
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


#ifndef NTUPLEOUTPUT_H
#define NTUPLEOUTPUT_H

/*
 * Useful info
 *
 * TNtuple is:
 *   branches of Float_t parameters, named somehow
 *
 * TTree is:
 *   branches of simple Float_t, Int_t parameters or complex classes, known to ROOT
 *   (new classes are given to ROOT by call gROOT->ProcessLine("#include <vector>"))
 *   the creation method for these 2 types is different,
 *   also these methods are different from from TNtuple ones (in better way)
 *
 */

/*
 * What is wanted:
 *   keep the output TTree interface (i.e. the definitions of Branches, their classes and names) in 1 file
 *   easily create or open TTree of this interface in main process
 *   loop over Entires
 *   and have full access to all the branches
 *
 * To do it:
 *   there are macro creating this interface
 *   and the list of them with current definition of the interface is in this file
 *   there are 2 bunches of macro -- for creating Branches of new TTree or for opening existing one
 *   which unfold into commands like
 *      Class_X NT_branchFoo; outputTTreeObject.Branch("branchFoo", "Class_X", &NT_branchFoo);
 *      or
 *      Class_X NT_branchFoo; outputTTreeObject.SetBranchAddress("branchFoo", &NT_branchFoo);
 *
 * the outputTTreeObject is defined in OUTNTUPLE
 * the mode of the interface (create or open ttree) is defined with NTUPLE_INTERFACE_CREATE or NTUPLE_INTERFACE_OPEN
 *
 * there are also a bunch of convenience macro for handling the TNtuple legacy bunches of Float_t parameters
 * -- they mostly should go away when propper objects are used
 *  and there is a macro reseting all the branch parameters -- it's ad-hoc, TODO: do it somehow in more automated, convenient way
 *
 * branch object names are prepended with NT_
 * so a branch named "foo" in the program namespace will have the object named NT_foo
 *
 * also default name of the TTree is NT_output_ttree
 *
 * there is a usage example in a comment further
 */

// default name of the output
#ifndef OUTNTUPLE
	#define OUTNTUPLE NT_output_ttree
#endif

/* macro declaring the object and setting a branch with its' pointer --- all in current, __not_global__ space (in main space)
 *
 * Notice the protocol:
 *    1) the object name in current namespace is `NT_Name`
 *    2) the branch name in the ntuple is `Name`
 */
#if defined(NTUPLE_INTERFACE_CREATE)
	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.Branch(#Name, #Class, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         Float_t NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/F")
	#define Int_t_in_NTuple(NTuple, Name)           Int_t   NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/I")
	#define Bool_t_in_NTuple(NTuple, Name)          Bool_t  NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/O")
#elif defined(NTUPLE_INTERFACE_OPEN)
	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.SetBranchAddress(#Name, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         OBJECT_in_NTuple(NTuple, Float_t, Name)
	#define Int_t_in_NTuple(NTuple, Name)           OBJECT_in_NTuple(NTuple, Int_t, Name)
	#define Bool_t_in_NTuple(NTuple, Name)          OBJECT_in_NTuple(NTuple, Bool_t, Name)
#else
	error: set ntuple interface mode
#endif


/*
 * this file ties the interface to our ntuple (TTree) in the current namespace
 *
 * Usage
 * 
 * actual files using it are:
 *   bin/ntupleEvents.cc           (creates the TTree, at line 1542 in main)
 *   test/likelihood_regions.cc    (uses existing TTree, at line 51 in pull_likelihood_regions)
 *
 * Roughly the idea is as follows
 *
 * declare your ntuple:
 *
 *     TTree output("reduced_ttree", "TTree with reduced event data");
 *     // if the name is not `NT_output_ttree` (which is assumed here in the interface)
 *     // define your name for preprocessor:
 *     #define OUTNTUPLE output
 *
 *     // set the mode of the interface to branches (create branches in a new TTree or open branches of existing TTree):
 *     #define NTUPLE_INTERFACE_CREATE
 * 
 *     // load this interface:
 *     #include "ntupleOutput.h"
 *
 * now you have NT_Name objects in the name space and the ntuple has branches "Name" with pointers to these objects
 * you can loop over TTree:
 *
 *     for (Long64_t i=0; i<NT_output_ttree.GetEntries(); i++)
 *         {
 *         NT_output_ttree.GetEntry(i);
 *         ...
 *         }
 *
 * copy objects from event (pseudocode):
 *
 *     NT_foo = events[i]["foo"])
 *
 * or actual example from ntupleEvents:
 *     NT_aMCatNLO_weight = evt->weight();
 * or from likelihood_regions:
 *     if (NT_tau_IDlev_0 != 3. && NT_tau_IDlev_1 != 3.) continue;
 *     
 * when done fill the ntuple:
 *
 *     output.Fill();
 *
 * clearing/reseting of the objects for each event -- currently it is responsibility of the programmer
 * but there is a sketchy macro for this now, it is in development
 * used as in ntupleEvents:
 *     RESET_NTUPLE_PARAMETERS // defaults all parameters
 *
 * -- with no ;
 * that's how sketchy it is
 *
 */


/* TODO: add what's missing from eventSelection
 *       and use proper TLorentzVector etc classes
 *
 * MET filters and lumisection certificate are done on the fly at ntuple production
 * lumi passes after MET filters -- to properly account for it in luminosity
 */

// the interface (all Float_ts, compatibility to first runs with TNtuple)
Float_t_in_NTuple(OUTNTUPLE, aMCatNLO_weight);
Float_t_in_NTuple(OUTNTUPLE, gen_t_pt);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_pt);
//Int_t_in_NTuple(OUTNTUPLE, gen_tt_decay_id);  // simple ID of TTbar: pdgId lep 1 * pdgId lep2 or * 1 for quark-decay
// it = product of the following two without signs
Int_t_in_NTuple(OUTNTUPLE, gen_t_w_decay_id); // = id of lepton (11/13/15, but the sign means which product is lepton: minus=1, plus=2) or 1 for quarks
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_phi);
Int_t_in_NTuple(OUTNTUPLE, gen_tb_w_decay_id);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_phi);
Int_t_in_NTuple(OUTNTUPLE, NUP_gen); // TODO: add gen info from TTbar
Int_t_in_NTuple(OUTNTUPLE, nvtx_gen);
Int_t_in_NTuple(OUTNTUPLE, nvtx);
Int_t_in_NTuple(OUTNTUPLE, nvtx_good);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetAll);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralNeutral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralChargedPileUp);

Bool_t_in_NTuple(OUTNTUPLE, HLT_el);
Bool_t_in_NTuple(OUTNTUPLE, HLT_mu);

Int_t_in_NTuple(OUTNTUPLE, leps_ID);
Int_t_in_NTuple(OUTNTUPLE, nleps);
Int_t_in_NTuple(OUTNTUPLE, njets);
Int_t_in_NTuple(OUTNTUPLE, nbjets);
Int_t_in_NTuple(OUTNTUPLE, ntaus);

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

OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_init);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_uncorrected);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_corrected);

Int_t_in_NTuple(OUTNTUPLE, lep0_id);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, lep0_p4);
Int_t_in_NTuple(OUTNTUPLE, lep1_id);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, lep1_p4);

#define JET_OUTPUT(num) \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_id); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_initial_p4); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_p4);  \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_uncorrected_p4); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_matched_genjet_p4); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_resolution); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf_up); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf_down); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_rad); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_pu_discr); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_b_discr); \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_hadronFlavour); \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_partonFlavour);

JET_OUTPUT(0)
JET_OUTPUT(1)
JET_OUTPUT(2)
JET_OUTPUT(3)
JET_OUTPUT(4)

#define TAU_OUTPUT(num) \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_id); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, tau##num##_p4); \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_IDlev);

TAU_OUTPUT(0)
TAU_OUTPUT(1)

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

#define NT_jet(i, jet, id_jet_p4, matched_genjet_p4, jet_radius_func, btagger_label) case i: { \
NT_jet ##i ##_id                = jet.pdgId(); \
NT_jet ##i ##_initial_p4        = id_jet_p4; \
NT_jet ##i ##_p4                = jet.p4(); \
NT_jet ##i ##_uncorrected_p4    = jet.correctedP4("Uncorrected"); \
NT_jet ##i ##_matched_genjet_p4 = matched_genjet_p4; \
NT_jet ##i ##_resolution  = jet.userFloat("jet_resolution"); \
NT_jet ##i ##_sf          = jet.userFloat("jer_sf"); \
NT_jet ##i ##_sf_up       = jet.userFloat("jer_sf_up"); \
NT_jet ##i ##_sf_down     = jet.userFloat("jer_sf_down"); \
NT_jet ##i ##_rad               = jet_radius_func(jet); \
NT_jet ##i ##_pu_discr          = jet.userFloat("pileupJetId:fullDiscriminant"); \
NT_jet ##i ##_b_discr           = jet.bDiscriminator(btagger_label); \
NT_jet ##i ##_hadronFlavour     = jet.hadronFlavour(); \
NT_jet ##i ##_partonFlavour     = jet.partonFlavour(); \
break; }

#define NT_tau(i, tau, IDlev) case i: { \
NT_tau ##i ##_id    = id; \
NT_tau ##i ##_p4    = tau.p4(); \
NT_tau ##i ##_IDlev = IDlev; \
break; }



// the automatic reset of all parameters for now
#define RESET_JET(num) \
jet##num##_id = -1; \
jet##num##_initial_p4.SetXYZT(0,0,0,0); \
jet##num##_p4.SetXYZT(0,0,0,0);  \
jet##num##_uncorrected_p4.SetXYZT(0,0,0,0); \
jet##num##_matched_genjet_p4.SetXYZT(0,0,0,0); \
jet##num##_resolution = -1; \
jet##num##_sf = -1; \
jet##num##_sf_up = -1; \
jet##num##_sf_down = -1; \
jet##num##_rad = -1; \
jet##num##_pu_discr = -1; \
jet##num##_b_discr = -1; \
jet##num##_hadronFlavour = -1; \
jet##num##_partonFlavour = -1;

#define RESET_TAU(num) \
tau##num##_id = -1; \
tau##num##_p4.SetXYZT(0,0,0,0); \
tau##num##_IDlev = -1;

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
NT_HLT_el = false; \
NT_HLT_mu = false; \
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
NT_lep0_p4.SetXYZT(0, 0, 0, 0); \
NT_lep1_id = -1; \
NT_lep1_p4.SetXYZT(0, 0, 0, 0); \
RESET_JET(0) \
RESET_JET(1) \
RESET_JET(2) \
RESET_JET(3) \
RESET_JET(4) \
RESET_TAU(0) \
RESET_TAU(1) \
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


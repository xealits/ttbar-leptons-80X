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
#elif defined(NTUPLE_INTERFACE_OPEN)
	#define OBJECT_in_NTuple(NTuple, Class, Name)   Class   NT_##Name; NTuple.SetBranchAddress(#Name, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         OBJECT_in_NTuple(NTuple, Float_t, Name)
	#define Int_t_in_NTuple(NTuple, Name)           OBJECT_in_NTuple(NTuple, Int_t, Name)
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

// the interface (all Float_ts, compatibility to first runs with TNtuple)
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
Float_t_in_NTuple(OUTNTUPLE, lep0_id);
Float_t_in_NTuple(OUTNTUPLE, lep0_eta);
Float_t_in_NTuple(OUTNTUPLE, lep0_phi);
Float_t_in_NTuple(OUTNTUPLE, lep0_pt);
Float_t_in_NTuple(OUTNTUPLE, lep0_p);
Float_t_in_NTuple(OUTNTUPLE, lep1_id);
Float_t_in_NTuple(OUTNTUPLE, lep1_eta);
Float_t_in_NTuple(OUTNTUPLE, lep1_phi);
Float_t_in_NTuple(OUTNTUPLE, lep1_pt);
Float_t_in_NTuple(OUTNTUPLE, lep1_p);
Float_t_in_NTuple(OUTNTUPLE, jet0_id);
Float_t_in_NTuple(OUTNTUPLE, jet0_eta);
Float_t_in_NTuple(OUTNTUPLE, jet0_phi);
Float_t_in_NTuple(OUTNTUPLE, jet0_pt);
Float_t_in_NTuple(OUTNTUPLE, jet0_p);
Float_t_in_NTuple(OUTNTUPLE, jet0_rad);
Float_t_in_NTuple(OUTNTUPLE, jet0_b_discr);
Float_t_in_NTuple(OUTNTUPLE, jet0_hadronFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet0_partonFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet1_id);
Float_t_in_NTuple(OUTNTUPLE, jet1_eta);
Float_t_in_NTuple(OUTNTUPLE, jet1_phi);
Float_t_in_NTuple(OUTNTUPLE, jet1_pt);
Float_t_in_NTuple(OUTNTUPLE, jet1_p);
Float_t_in_NTuple(OUTNTUPLE, jet1_rad);
Float_t_in_NTuple(OUTNTUPLE, jet1_b_discr);
Float_t_in_NTuple(OUTNTUPLE, jet1_hadronFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet1_partonFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet2_id);
Float_t_in_NTuple(OUTNTUPLE, jet2_eta);
Float_t_in_NTuple(OUTNTUPLE, jet2_phi);
Float_t_in_NTuple(OUTNTUPLE, jet2_pt);
Float_t_in_NTuple(OUTNTUPLE, jet2_p);
Float_t_in_NTuple(OUTNTUPLE, jet2_rad);
Float_t_in_NTuple(OUTNTUPLE, jet2_b_discr);
Float_t_in_NTuple(OUTNTUPLE, jet2_hadronFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet2_partonFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet3_id);
Float_t_in_NTuple(OUTNTUPLE, jet3_eta);
Float_t_in_NTuple(OUTNTUPLE, jet3_phi);
Float_t_in_NTuple(OUTNTUPLE, jet3_pt);
Float_t_in_NTuple(OUTNTUPLE, jet3_p);
Float_t_in_NTuple(OUTNTUPLE, jet3_rad);
Float_t_in_NTuple(OUTNTUPLE, jet3_b_discr);
Float_t_in_NTuple(OUTNTUPLE, jet3_hadronFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet3_partonFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet4_id);
Float_t_in_NTuple(OUTNTUPLE, jet4_eta);
Float_t_in_NTuple(OUTNTUPLE, jet4_phi);
Float_t_in_NTuple(OUTNTUPLE, jet4_pt);
Float_t_in_NTuple(OUTNTUPLE, jet4_p);
Float_t_in_NTuple(OUTNTUPLE, jet4_rad);
Float_t_in_NTuple(OUTNTUPLE, jet4_b_discr);
Float_t_in_NTuple(OUTNTUPLE, jet4_hadronFlavour);
Float_t_in_NTuple(OUTNTUPLE, jet4_partonFlavour);
Float_t_in_NTuple(OUTNTUPLE, tau0_id);
Float_t_in_NTuple(OUTNTUPLE, tau0_eta);
Float_t_in_NTuple(OUTNTUPLE, tau0_phi);
Float_t_in_NTuple(OUTNTUPLE, tau0_pt);
Float_t_in_NTuple(OUTNTUPLE, tau0_p);
Float_t_in_NTuple(OUTNTUPLE, tau0_IDlev);
Float_t_in_NTuple(OUTNTUPLE, tau1_id);
Float_t_in_NTuple(OUTNTUPLE, tau1_eta);
Float_t_in_NTuple(OUTNTUPLE, tau1_phi);
Float_t_in_NTuple(OUTNTUPLE, tau1_pt);
Float_t_in_NTuple(OUTNTUPLE, tau1_p);
Float_t_in_NTuple(OUTNTUPLE, tau1_IDlev);

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


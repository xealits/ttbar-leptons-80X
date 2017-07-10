/*
 * definitions for attaching branches to existing ttree
 * or declaring branches for new ttree
 * etc
 */
#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"

#ifndef NTUPLEOUTPUT_H
#define NTUPLEOUTPUT_H

// the full interface (all Float_ts, compatibility to first runs with TNtuple)
#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_common.h"
//Int_t_in_NTuple(OUTNTUPLE, gen_tt_decay_id);  // simple ID of TTbar: pdgId lep 1 * pdgId lep2 or * 1 for quark-decay
// it = product of the following two without signs
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_t_w_p2_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p1_phi);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_eta);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_w_p2_phi);
Int_t_in_NTuple(OUTNTUPLE, nvtx_gen);
Int_t_in_NTuple(OUTNTUPLE, nvtx_good);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetAll);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralNeutral);
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetCentralChargedPileUp);

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

Int_t_in_NTuple(OUTNTUPLE, lep0_id);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, lep0_p4);
Float_t_in_NTuple(OUTNTUPLE, lep0_dxy);
Float_t_in_NTuple(OUTNTUPLE, lep0_dz);
Float_t_in_NTuple(OUTNTUPLE, lep0_relIso);
Float_t_in_NTuple(OUTNTUPLE, lep0_dB);
Int_t_in_NTuple(OUTNTUPLE, lep1_id);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, lep1_p4);
Float_t_in_NTuple(OUTNTUPLE, lep1_dxy);
Float_t_in_NTuple(OUTNTUPLE, lep1_dz);
Float_t_in_NTuple(OUTNTUPLE, lep1_relIso);
Float_t_in_NTuple(OUTNTUPLE, lep1_dB);

#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_jets.h"

#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_taus.h"

// convenience:
#define NT_tau_secondaryVertexCov_(i, j) NT_tau_secondaryVertexCov_ ##i ##j

#define NT_LEPTONS_N 2



#define RESET_NTUPLE_PARAMETERS \
RESET_COMMON \
NT_nvtx_gen = -1; \
NT_nvtx_good = -1; \
NT_fixedGridRhoFastjetAll = -1; \
NT_fixedGridRhoFastjetCentral = -1; \
NT_fixedGridRhoFastjetCentralNeutral = -1; \
NT_fixedGridRhoFastjetCentralChargedPileUp = -1; \
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
NT_lep0_dxy = -1; \
NT_lep0_dz  = -1; \
NT_lep0_relIso  = -1; \
NT_lep0_dB  = -1; \
NT_lep1_id = -1; \
NT_lep1_p4.SetXYZT(0, 0, 0, 0); \
NT_lep1_dxy = -1; \
NT_lep1_dz  = -1; \
NT_lep1_relIso  = -1; \
NT_lep1_dB  = -1; \
RESET_JETS \
RESET_TAUS \
NT_gen_t_w_p1_eta = -1; \
NT_gen_t_w_p1_phi = -1; \
NT_gen_t_w_p2_eta = -1; \
NT_gen_t_w_p2_phi = -1; \
NT_gen_tb_w_p1_eta = -1; \
NT_gen_tb_w_p1_phi = -1; \
NT_gen_tb_w_p2_eta = -1; \
NT_gen_tb_w_p2_phi = -1;

#endif /* NTUPLEOUTPUT_H */


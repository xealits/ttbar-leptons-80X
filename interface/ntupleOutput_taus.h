#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"

#ifndef NTUPLEOUTPUT_TAUS_H
#define NTUPLEOUTPUT_TAUS_H

#define NT_TAUS_N 2

#define TAU_OUTPUT(num) \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_id); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, tau##num##_p4); \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_IDlev); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leading_track_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadChargedHadrCand_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadNeutralCand_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadCand_pt);

TAU_OUTPUT(0)
TAU_OUTPUT(1)

#define NT_tau(i, tau, IDlev) case i: { \
NT_tau ##i ##_id    = tau.pdgId(); \
NT_tau ##i ##_p4    = tau.p4(); \
NT_tau ##i ##_IDlev = IDlev; \
NT_tau ##i ##_leading_track_pt = tau.userFloat("leading_track_pt"); \
NT_tau ##i ##_leadChargedHadrCand_pt = tau.userFloat("leadChargedHadrCand_pt"); \
NT_tau ##i ##_leadNeutralCand_pt     = tau.userFloat("leadNeutralCand_pt"); \
NT_tau ##i ##_leadCand_pt            = tau.userFloat("leadCand_pt"); \
break; }

#define RESET_TAU(num) \
NT_tau##num##_id = -1; \
NT_tau##num##_p4.SetXYZT(0,0,0,0); \
NT_tau##num##_IDlev = -1; \
NT_tau##num##_leading_track_pt = -1; \
NT_tau##num##_leadChargedHadrCand_pt = -1; \
NT_tau##num##_leadNeutralCand_pt = -1; \
NT_tau##num##_leadCand_pt = -1;

#define RESET_TAUS \
RESET_TAU(0) \
RESET_TAU(1)

#endif /* NTUPLEOUTPUT_TAUS_H */


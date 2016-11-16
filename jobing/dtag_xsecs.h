#ifndef DTAG_XSECS_H_   /* Include guard */
#define DTAG_XSECS_H_


#include <map>

#include "TString.h"

using namespace std;


const double W_lep_br = 0.108;
const double W_qar_br = 0.676;

const double W_lep_br2 = 0.108*0.108;
const double W_qar_br2 = 0.676*0.676;

const double ttbar_xsec = 831.76;

std::map<TString, double> xsecs = {
{"Data13TeV_SingleElectron2016B_PromptRecoV2" , 1.0},
{"Data13TeV_SingleMuon2016B_PromptRecoV2" , 1.0},
{"Data13TeV_SingleElectron2016D_PromptRecoV2_", 1},
{"Data13TeV_SingleMuon2016D_PromptRecoV2_", 1},

{ "MC2016_noHLT_TTJets_powheg_scaleup_elelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elmubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_eltaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mumubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mutaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qqbar"           , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauelbar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_taumubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tautaubar"       , ttbar_xsec * W_lep_br2 / 2},

{ "MC2016_noHLT_TTJets_powheg_scaledown_elelbar"       , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elmubar"       , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elqbar"        , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_eltaubar"      , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muelbar"       , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mumubar"       , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muqbar"        , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mutaubar"      , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qelbar"        , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qmubar"        , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qqbar"         , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qtaubar"       , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauelbar"      , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_taumubar"      , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauqbar"       , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tautaubar"     , ttbar_xsec * W_lep_br2 / 2},

{ "MC2016_withHLT_TTJets_up_powheg_elelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_elmubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_eltaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_muelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_mumubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_mutaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_qqbar"           , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_tauelbar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_taumubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_up_powheg_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_up_powheg_tautaubar"       , ttbar_xsec * W_lep_br2 / 2},

{ "MC2016_withHLT_TTJets_down_powheg_elelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_elmubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_eltaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_muelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_mumubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_mutaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_qqbar"           , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_tauelbar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_taumubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_withHLT_TTJets_down_powheg_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_withHLT_TTJets_down_powheg_tautaubar"       , ttbar_xsec * W_lep_br2 / 2},

{ "MC2016_reHLT_TTJets_powheg_elelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_elmubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_eltaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_muelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_mumubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_mutaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_qqbar"           , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_tauelbar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_taumubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_powheg_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_powheg_tautaubar"       , ttbar_xsec * W_lep_br2 / 2},

{ "MC2016_reHLT_TTJets_amcatnlo_elelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_elmubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_eltaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_muelbar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_mumubar"         , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_mutaubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_qqbar"           , ttbar_xsec * W_qar_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_tauelbar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_taumubar"        , ttbar_xsec * W_lep_br2 / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br / 2},
{ "MC2016_reHLT_TTJets_amcatnlo_tautaubar"       , ttbar_xsec * W_lep_br2 / 2},

{"MC2016_noHLT_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W1Jets_madgraph", 9493},
{"MC2016_noHLT_W2Jets_madgraph", 3120},
{"MC2016_noHLT_W3Jets_madgrapg", 942.3},
{"MC2016_noHLT_W4Jets_madgraph", 524.2},
{"MC2016_reHLT_WJets_amcatnlo", 61526.7},
{"MC2016_reHLT_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_reHLT_DYJetsToLL_50toInf_amcatnlo", 6025.2},
{"MC2016_noHLT_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_noHLT_DYJetsToLL_50toInf_amcatnlo", 6025.2},
{"MC2016_noHLT_WW", 113.89},
{"MC2016_noHLT_WZ", 47.13},
{"MC2016_noHLT_ZZ", 16.52},
{"MC2016_noHLT_SingleT_tW_5FS_powheg",    35.6},
{"MC2016_noHLT_SingleTbar_tW_5FS_powheg", 35.6},
{"MC2016_noHLT_schannel_4FS_leptonicDecays_amcatnlo", 3.36},
{"MC2016_noHLT_tchannel_antitop_4f_leptonicDecays_powheg", 70.69/2},
{"MC2016_noHLT_tchannel_top_4f_leptonicDecays_powheg", 70.69/2},

{"MC2016_noHLT_QCD_HT-100-200",  27540000},
{"MC2016_noHLT_QCD_HT-200-300",  1717000},
{"MC2016_noHLT_QCD_HT-300-500",  351300},
{"MC2016_noHLT_QCD_HT-500-700",  31630},
{"MC2016_noHLT_QCD_HT-700-1000",  6802},
{"MC2016_noHLT_QCD_HT-1000-1500",  1206},
{"MC2016_noHLT_QCD_HT-1500-2000",  120.4},
{"MC2016_noHLT_QCD_HT-2000-Inf",  25.25},
};

#endif // DTAG_XSECS_H_
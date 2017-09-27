#ifndef DTAG_XSECS_H_   /* Include guard */
#define DTAG_XSECS_H_

#include <map>
#include "TString.h"

using namespace std;


const double W_lep_br = 0.108;
const double W_qar_br = 0.676;

const double W_lep_br2 = W_lep_br*W_lep_br;
const double W_qar_br2 = W_qar_br*W_qar_br;

const double br_tau_electron = 0.1785;
const double br_tau_muon     = 0.1736;
const double br_tau_lepton   = br_tau_electron + br_tau_muon;
const double br_tau_hadronic = 1 - br_tau_lepton;

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

{ "MC2016_reHLT_TTJets_powheg_elelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_elmubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_eltaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_muelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_mumubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_mutaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_qqbar"           , ttbar_xsec * W_qar_br2 },
{ "MC2016_reHLT_TTJets_powheg_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_tauelbar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_taumubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_powheg_tautaubar"       , ttbar_xsec * W_lep_br2 },

{ "MC2016_reHLT_TTJets_amcatnlo_elelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_elmubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_eltaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_muelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_mumubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_mutaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_qqbar"           , ttbar_xsec * W_qar_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_tauelbar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_taumubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_tautaubar"       , ttbar_xsec * W_lep_br2 },

{ "MC2016_reHLT_TTJets_amcatnlo_elelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_elmubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_elqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_eltaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_muelbar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_mumubar"         , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_muqbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_mutaubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_qelbar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_qmubar"          , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_qqbar"           , ttbar_xsec * W_qar_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_qtaubar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_tauelbar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_taumubar"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_amcatnlo_tauqbar"         , ttbar_xsec * W_qar_br*W_lep_br },
{ "MC2016_reHLT_TTJets_amcatnlo_tautaubar"       , ttbar_xsec * W_lep_br2 },

{ "MC2016_reHLT_TTJets_powheg_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_reHLT_TTJets_powheg_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_reHLT_TTJets_powheg_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_reHLT_TTJets_powheg_eell"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_elmu"        , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_reHLT_TTJets_powheg_elq"         , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_reHLT_TTJets_powheg_mmuu"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_reHLT_TTJets_powheg_mqu"         , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_reHLT_TTJets_powheg_qq"          , ttbar_xsec * W_qar_br2 },

// mcm xsecs
{"MC2016_noHLT_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_noHLT_W1Jets_madgraph", 9493},
{"MC2016_noHLT_W2Jets_madgraph", 3120},
{"MC2016_noHLT_W3Jets_madgrapg", 942.3},
{"MC2016_noHLT_W3Jets_madgraph", 942.3},
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
{"MC2016_noHLT_QCD_HT-200-300",   1717000},
{"MC2016_noHLT_QCD_HT-300-500",    351300},
{"MC2016_noHLT_QCD_HT-500-700",     31630},
{"MC2016_noHLT_QCD_HT-700-1000",     6802},
{"MC2016_noHLT_QCD_HT-1000-1500",    1206},
{"MC2016_noHLT_QCD_HT-1500-2000",   120.4},
{"MC2016_noHLT_QCD_HT-2000-Inf",    25.25},

// mcm xsecs
{"MC2016_Summer16_WJets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_Summer16_WJets_amcatnlo_full", 61526.7},
{"MC2016_Summer16_W0Jets_amcatnlo", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_Summer16_WJets_madgraph_075", 50690 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_Summer16_W0Jets_madgraph", 50690 - 9493 - 3120 - 942.3 - 524.2},

{"MC2016_Summer16_WJets_madgraph_perc", 37713.4 },
{"MC2016_Summer16_W1Jets_madgraph_perc", 8921.44},
{"MC2016_Summer16_W2Jets_madgraph_perc", 2838.64},
{"MC2016_Summer16_W3Jets_madgraph_perc",  861.73},
{"MC2016_Summer16_W4Jets_madgraph_perc",  354.83},

// since _perc ~ worked doing this way
// no NUP, but reduced xsec
{"MC2016_Summer16_WJets_madgraph", 61526.7}, // 50690 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_Summer16_W1Jets_madgraph", 9493},
{"MC2016_Summer16_W2Jets_madgraph", 3120},
{"MC2016_Summer16_W3Jets_madgraph", 942.3},
{"MC2016_Summer16_W4Jets_madgraph", 524.2},

// m-i-t scales
{"MC2016_Summer16_mScale_WJets_full_madgraph",  61526},
{"MC2016_Summer16_mScale_W0Jets_full_madgraph", 61526},
{"MC2016_Summer16_mScale_WJets_madgraph", 34273.6},
{"MC2016_Summer16_mScale_W0Jets_madgraph", 34273.6},
{"MC2016_Summer16_mScale_W1Jets_madgraph", 18455.98},
{"MC2016_Summer16_mScale_W2Jets_madgraph", 6035.9},
{"MC2016_Summer16_mScale_W3Jets_madgraph", 1821.5},
{"MC2016_Summer16_mScale_W4Jets_madgraph", 939.7 },

// xsecs computed from inclusive set divided by NUP
{"MC2016_Summer16_fixedXsec_W0Jets_amcatnlo",  1355.7  },
{"MC2016_Summer16_fixedXsec_W1Jets_madgraph",  4769.65 },
{"MC2016_Summer16_fixedXsec_W2Jets_madgraph", 36689.6  },
{"MC2016_Summer16_fixedXsec_W3Jets_madgraph", 18711.7  },
{"MC2016_Summer16_fixedXsec_W4Jets_madgraph",   524.2  }, // left the same

{"MC2016_Summer16_WJets_HT70to100_madgraph"    , 1319 },
{"MC2016_Summer16_WJets_HT100to200_madgraph"   , 1345 },
{"MC2016_Summer16_WJets_HT200to400_madgraph"   , 359.7},
{"MC2016_Summer16_WJets_HT400to600_madgraph"   , 48.91},
{"MC2016_Summer16_WJets_HT600to800_madgraph"   , 12.05},
{"MC2016_Summer16_WJets_HT800to1200_madgraph"  , 5.501},
{"MC2016_Summer16_WJets_HT1200to2500_madgraph" , 1.329},
{"MC2016_Summer16_WJets_HT2500toInf_madgraph"  , 0.03216},

{"MC2016_Summer16_WJetsToLNu_HT-70To100"    , 1319 },
{"MC2016_Summer16_WJetsToLNu_HT-100To200"   , 1345 },
{"MC2016_Summer16_WJetsToLNu_HT-200To400"   , 359.7},
{"MC2016_Summer16_WJetsToLNu_HT-400To600"   , 48.91},
{"MC2016_Summer16_WJetsToLNu_HT-600To800"   , 12.05},
{"MC2016_Summer16_WJetsToLNu_HT-800To1200"  , 5.501},
{"MC2016_Summer16_WJetsToLNu_HT-1200To2500" , 1.329},
{"MC2016_Summer16_WJetsToLNu_HT-2500ToInf"  , 0.03216},

{"DYJetsToLL_M-5to50_HT-100to200", 224.2 },
{"DYJetsToLL_M-5to50_HT-200to400", 37.2 },
{"DYJetsToLL_M-5to50_HT-400to600", 3.581 },
{"DYJetsToLL_M-5to50_HT-600toInf", 1.124 },

{"DYJetsToLL_M-50_HT-70to100", 175.3 },
{"DYJetsToLL_M-50_HT-100to200", 147.4 },
{"DYJetsToLL_M-50_HT-200to400", 40.99 },
{"DYJetsToLL_M-50_HT-400to600", 5.678 },
{"DYJetsToLL_M-50_HT-600to800", 1.367 },
{"DYJetsToLL_M-50_HT-800to1200", 0.6304 },
{"DYJetsToLL_M-50_HT-1200to2500", 0.1514 },
{"DYJetsToLL_M-50_HT-2500toInf", 0.003565 },

{"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo", 18610},
{"MC2016_Summer16_DYJetsToLL_50toInf_madgraph", 5765.4}, // but McM shows 4970
{"MC2016_Summer16_DYJetsToTauTau_50toInf_amcatnlo", 1867}, // from McM TODO: check

{ "MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin"      , 3.22   },
//{ "MC2016_Summer16_ZZTo2L2Nu_powheg"               , 0.564  },
{ "MC2016_Summer16_ZZTo2L2Nu_powheg"               , 1.256  },
{ "MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin"      , 5.595  },
{ "MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin"     , 3.033  },
{ "MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin"   , 10.71  },
{ "MC2016_Summer16_WZTo3LNu_powheg"                , 4.42965  },
{ "MC2016_Summer16_WWToLNuQQ_powheg"               ,  49.997  },
{ "MC2016_Summer16_WWTo2L2Nu_powheg"               ,  12.178  },

{"MC2016_Summer16_QCD_HT-2000-Inf",  25.25},
{"MC2016_Summer16_QCD_HT-1500-2000",  120.4},
{"MC2016_Summer16_QCD_HT-1000-1500",  1206},
{"MC2016_Summer16_QCD_HT-700-1000",  6802},
{"MC2016_Summer16_QCD_HT-500-700",  31630},
{"MC2016_Summer16_QCD_HT-300-500",  351300},
{"MC2016_Summer16_QCD_HT-200-300",  1717000},
{"MC2016_Summer16_QCD_HT-100-200",  27540000},
{"MC2016_Summer16_QCD_HT-50-100", 0. },      // FIXME: xsection for 50-100 QCD (McM shows xsec = 1.0 /pb)
//{"MC2016_Summer16_QCD_HT-50-100", 27540000.*8 }, // FIXME: just trying out the effect from this dset
{"MC2016_Summer16_QCD_HT-50-100", 0. }, // FIXME: just trying out the effect from this dset

{"MC2016_Summer16_QCD_HT-2000-Inf_",  25.25},
{"MC2016_Summer16_QCD_HT-1500-2000_",  120.4},
{"MC2016_Summer16_QCD_HT-1000-1500_",  1206},
{"MC2016_Summer16_QCD_HT-700-1000_",  6802},
{"MC2016_Summer16_QCD_HT-500-700_",  31630},
{"MC2016_Summer16_QCD_HT-300-500_",  351300},
{"MC2016_Summer16_QCD_HT-200-300_",  1717000},
{"MC2016_Summer16_QCD_HT-100-200_",  27540000},
//{"MC2016_Summer16_QCD_HT-50-100_", 27540000.*8 }, // FIXME: just trying out the effect from this dset
{"MC2016_Summer16_QCD_HT-50-100_", 0. }, // FIXME: xsection for 50-100 QCD

{"MC2016_Summer16_W0Jets_amcatnlo_", 61526.7 - 9493 - 3120 - 942.3 - 524.2},
{"MC2016_Summer16_W1Jets_madgraph_", 9493},
{"MC2016_Summer16_W2Jets_madgraph_", 3120},
{"MC2016_Summer16_W3Jets_madgrapg_", 942.3},
{"MC2016_Summer16_W3Jets_madgraph_", 942.3},
{"MC2016_Summer16_W4Jets_madgraph_", 524.2},

{"MC2016_Summer16_WJets_amcatnlo_", 61526.7},
{"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_", 18610},
//{"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_", 6025.2}, // FIXME: the DY x-section was updated by about 5%
{"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_", 4970}, // FIXME: the DY x-section was updated by about 5%

{"MC2016_Summer16_SingleT_tW_5FS_powheg_",    35.6},
{"MC2016_Summer16_SingleTbar_tW_5FS_powheg_", 35.6},
{"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo_", 10.11}, //3.36},
{"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg_", 80.95}, //70.69/2},
{"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg_", 136.02}, //70.69/2},

{"MC2016_Summer16_SingleT_tW_5FS_powheg",    35.6},
{"MC2016_Summer16_SingleTbar_tW_5FS_powheg", 35.6},
{"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo", 10.11}, //3.36},
{"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg", 80.95}, //70.69/2},
{"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg", 136.02}, //70.69/2},

// inclusive sample
{ "MC2016_Summer16_TTJets_powheg"      , ttbar_xsec },
{ "MC2016_Summer16_TTJets_powheg_1"      , ttbar_xsec },
{ "MC2016_Summer16_TTJets_powheg_2"      , ttbar_xsec },
{ "MC2016_Summer16_TTJets_powheg_3"      , ttbar_xsec },
{ "MC2016_Summer16_TTJets_powheg_4"      , ttbar_xsec },
{ "MC2016_Summer16_TTJets_powheg_5"      , ttbar_xsec },

/* elmu fakerates 2
 */
{ "MC2016_Summer16_TTJets_powheg_2_aelmtuu" ,  ttbar_xsec * W_lep_br2 * 2 * br_tau_lepton },
{ "MC2016_Summer16_TTJets_powheg_2_elmu"    ,  ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_2_other"   ,  ttbar_xsec * (1 - 2*W_lep_br2 - 2*W_lep_br2*br_tau_lepton) },

/* inclusive tau decays
*/
{ "MC2016_Summer16_TTJets_powheg_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},

// 2 tau decays
{ "MC2016_Summer16_TTJets_powheg_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TTJets_powheg_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TTJets_powheg_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TTJets_powheg_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TTJets_powheg_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TTJets_powheg_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},

// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TTJets_powheg_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },

{ "MC2016_Summer16_TTJets_powheg_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TTJets_powheg_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },

{ "MC2016_Summer16_TTJets_powheg_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TTJets_powheg_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },

// hadronic tau decays
{ "MC2016_Summer16_TTJets_powheg_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TTJets_powheg_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TTJets_powheg_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},

/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TTJets_powheg_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TTJets_powheg_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TTJets_powheg_qq"        , ttbar_xsec * W_qar_br2 },

// herwigpp
// inclusive sample
{ "MC2016_Summer16_TT_powheg_herwigpp"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_herwigpp_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_herwigpp_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_herwigpp_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_herwigpp_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_herwigpp_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_herwigpp_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_herwigpp_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_herwigpp_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_herwigpp_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_herwigpp_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_herwigpp_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_herwigpp_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_herwigpp_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_herwigpp_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_herwigpp_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_herwigpp_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_herwigpp_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_herwigpp_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_herwigpp_qq"        , ttbar_xsec * W_qar_br2 },


// fsrdown
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_fsrdown_qq"        , ttbar_xsec * W_qar_br2 },

// fsrup
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_fsrup"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_fsrup_qq"        , ttbar_xsec * W_qar_br2 },

// isrdown
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_isrdown"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_isrdown_qq"        , ttbar_xsec * W_qar_br2 },

// isr up
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_isrup"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_isrup_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_isrup_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_isrup_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_isrup_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_isrup_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_isrup_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_isrup_qq"        , ttbar_xsec * W_qar_br2 },

// hdumpdown
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_hdampDOWN_qq"        , ttbar_xsec * W_qar_br2 },

// hdampUP
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_hdampUP_qq"        , ttbar_xsec * W_qar_br2 },

// GluonMoveCRTune
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_GluonMoveCRTune_qq"        , ttbar_xsec * W_qar_br2 },

// QCDbasedCRTune
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_QCDbasedCRTune_qq"        , ttbar_xsec * W_qar_br2 },

// TuneCUETP8M2T4down
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4down_qq"        , ttbar_xsec * W_qar_br2 },

// TuneCUETP8M2T4up
// inclusive sample
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up"      , ttbar_xsec },
/* inclusive tau decays
*/
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br * 2},
// 2 tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aaelmttuuu",  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aaeellttuu",  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aaehlttuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aammttuuuu",  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aahmttuuu" ,  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aahhttuu"  ,  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic},
// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aelmtuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aeelltu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_ammtuuu"   ,  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aelqtu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_amqtuu"    ,  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon },
// hadronic tau decays
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_aehltu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_ahmtuu"     , ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_ahqtu"      , ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic},
/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_eell"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_elmu"      , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_elq"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_mmuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_mqu"       , ttbar_xsec * W_lep_br*W_qar_br * 2},
{ "MC2016_Summer16_TT_powheg_pythia_TuneCUETP8M2T4up_qq"        , ttbar_xsec * W_qar_br2 },

};

std::map<TString, double> xsecs_inclusive_tt = {
/* inclusive tau decays
*/
{ "MC2016_Summer16_TTJets_powheg_aattuu"      , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aeltu"       , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_amtuu"       , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aqtu"        , ttbar_xsec},

// 2 tau decays
{ "MC2016_Summer16_TTJets_powheg_aaelmttuuu",  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aaeellttuu",  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aaehlttuu" ,  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aammttuuuu",  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aahmttuuu" ,  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_aahhttuu"  ,  ttbar_xsec},

// 1 tau decays
// leptonic tau decays
{ "MC2016_Summer16_TTJets_powheg_aelmtuu"   ,  ttbar_xsec},

{ "MC2016_Summer16_TTJets_powheg_aeelltu"   ,  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_ammtuuu"   ,  ttbar_xsec},

{ "MC2016_Summer16_TTJets_powheg_aelqtu"    ,  ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_amqtuu"    ,  ttbar_xsec},

// hadronic tau decays
{ "MC2016_Summer16_TTJets_powheg_aehltu"     , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_ahmtuu"     , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_ahqtu"      , ttbar_xsec},

/* no taus, taken care by inclusive tau dtags
 */
{ "MC2016_Summer16_TTJets_powheg_eell"      , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_elmu"      , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_elq"       , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_mmuu"      , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_mqu"       , ttbar_xsec},
{ "MC2016_Summer16_TTJets_powheg_qq"        , ttbar_xsec}
};



// nick and colour for dtags
std::pair<TString, Color_t> dtag_nick_colour_old(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("dyjets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("wjets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("singletop", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qqbar")) return std::make_pair("tt_jj", kGreen+4);
		else if (dtag.Contains("elqbar") || dtag.Contains("qelbar") ||dtag.Contains("muqbar") || dtag.Contains("qmubar") || dtag.Contains("tauqbar") || dtag.Contains("qtaubar")) return std::make_pair("tt_lj", kGreen+3);
		else if (dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt_em", kYellow-7);
		else if (dtag.Contains("elelbar")) return std::make_pair("tt_ee", kAzure-9);
		else if (dtag.Contains("mumubar")) return std::make_pair("tt_mm", kGreen-9);
		else if (dtag.Contains("eltaubar") || dtag.Contains("tauelbar")) return std::make_pair("tt_et", kOrange+4);
		else if (dtag.Contains("mutaubar") || dtag.Contains("taumubar")) return std::make_pair("tt_mt", kOrange+1);
		else return std::make_pair("tt_other", kYellow+1);
		}
	else if(dtag.Contains("QCD")) return std::make_pair("qcd", kViolet);
	else return std::make_pair("other", kBlack);

	}

// nick and colour for dtags
std::pair<TString, Color_t> dtag_nick_colour(TString dtag)
	{
	if (dtag.Contains("Data")) return std::make_pair("data", kWhite);
	else if(dtag.Contains("DYJets")) return std::make_pair("dyjets", kGray);
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return std::make_pair("w-jets", kRed+1);
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return std::make_pair("dibosons", kCyan);
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return std::make_pair("single top", kAzure);
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qqbar")) return std::make_pair("tt_jj", kGreen+4);
		else if (dtag.Contains("elq") || dtag.Contains("mqu") || dtag.Contains("elqbar") || dtag.Contains("qelbar") ||dtag.Contains("muqbar") || dtag.Contains("qmubar") || dtag.Contains("tauqbar") || dtag.Contains("qtaubar")) return std::make_pair("tt_lj", kGreen+3);
		else if (dtag.Contains("elmu") || dtag.Contains("elmubar") || dtag.Contains("muelbar")) return std::make_pair("tt_em", kYellow-7);
		else if (dtag.Contains("elelbar")) return std::make_pair("tt_ee", kAzure-9);
		else if (dtag.Contains("mumubar")) return std::make_pair("tt_mm", kGreen-9);
		else if (dtag.Contains("aeltu"))
			return std::make_pair("tt_{e\\tau}", kOrange+2);
		else if (dtag.Contains("amtuu"))
			return std::make_pair("tt_{\\mu\\tau}", kOrange+1);
		else if (dtag.Contains("aelmtuu"))
			return std::make_pair("tt_{l\\tau-l}", kOrange+3);
		else return std::make_pair("tt_{other}", kCyan-5);
		}
	else if(dtag.Contains("QCD")) return std::make_pair("qcd", kViolet);
	else return std::make_pair("other", kBlack);
	}

// nick and colour for dtags
TString dtag_nick(TString dtag)
	{
	if (dtag.Contains("Data")) return "data";
	else if(dtag.Contains("DYJets")) return "dyjets";
	else if(dtag.Contains("W0Jets") ||dtag.Contains("W4Jets") ||dtag.Contains("W1Jets") ||dtag.Contains("W2Jets") ||dtag.Contains("W3Jets") ||dtag.Contains("WJets") ) return "w-jets";
	else if(dtag.Contains("WW") ||dtag.Contains("WZ") ||dtag.Contains("ZZ")) return "dibosons";
	else if(dtag.Contains("Single") || dtag.Contains("schannel") ||dtag.Contains("tchannel")) return "single top";
	else if(dtag.Contains("TT"))
		{
		if (dtag.Contains("qqbar")) return "tt_jj";
		else if (dtag.Contains("elq") || dtag.Contains("mqu") || dtag.Contains("elqbar") || dtag.Contains("qelbar") ||dtag.Contains("muqbar") || dtag.Contains("qmubar") || dtag.Contains("tauqbar") || dtag.Contains("qtaubar")) return "tt_lj";
		else if (dtag.Contains("elmu") || dtag.Contains("elmubar") || dtag.Contains("muelbar")) return "tt_em";
		else if (dtag.Contains("elelbar")) return "tt_ee";
		else if (dtag.Contains("mumubar")) return "tt_mm";
		else if (dtag.Contains("aeltu"))
			return "tt_{e\\tau}";
		else if (dtag.Contains("amtuu"))
			return "tt_{\\mu\\tau}";
		else if (dtag.Contains("aelmtuu"))
			return "tt_{l\\tau-l}";
		else return "tt_{other}";
		}
	else if(dtag.Contains("QCD")) return "qcd";
	else return "other";
	}

// nick and colour for nicks
Color_t nick_colour(TString nick)
	{
	if (nick.Contains("data")) return kWhite;
	else if(nick.Contains("dyjets")) return kGray;
	else if(nick.Contains("w-jets")) return kRed+1;
	else if(nick.Contains("dibosons")) return kCyan;
	else if(nick.Contains("single top")) return kAzure;

	else if (nick.Contains("tt_jj")) return kGreen+4;
	else if (nick.Contains("tt_lj")) return kGreen+3;
	else if (nick.Contains("tt_em")) return kYellow-7;
	else if (nick.Contains("tt_{em}")) return kYellow-7;
	else if (nick.Contains("tt_ee")) return kAzure-9;
	else if (nick.Contains("tt_mm")) return kGreen-9;
	else if (nick.Contains("tt_{e\\tau}")) return kOrange+2;
	else if (nick.Contains("tt_{\\mu\\tau}")) return kOrange+1;
	else if (nick.Contains("tt_{l\\tau-l}"))  return kOrange+3;
	else if (nick.Contains("tt_{other}"))     return kCyan-5;

	else if(nick.Contains("qcd")) return kViolet;
	else return kBlack;
	}

/* Include it to colors somehow, to preserve inclusiveness of ttbar in jet fake rates scripts
{ "MC2016_Summer16_TTJets_powheg_aattuu"      , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_aeltu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_amtuu"       , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_aqtu"        , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_Summer16_TTJets_powheg_eell"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_elmu"        , ttbar_xsec * W_lep_br2 * 2 },
{ "MC2016_Summer16_TTJets_powheg_elq"         , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_Summer16_TTJets_powheg_mmuu"        , ttbar_xsec * W_lep_br2 },
{ "MC2016_Summer16_TTJets_powheg_mqu"         , ttbar_xsec * W_lep_br*W_qar_br },
{ "MC2016_Summer16_TTJets_powheg_qq"          , ttbar_xsec * W_qar_br2 },
*/


#endif // DTAG_XSECS_H_

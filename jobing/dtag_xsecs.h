#ifndef DTAG_XSECS_H_   /* Include guard */
#define DTAG_XSECS_H_

#include <map>

using namespace std;


const double W_lep_br = 0.108;
const double W_qar_br = 0.676;

double W_lep_br2 = 0.108*0.108;
double W_qar_br2 = 0.676*0.676;

const double ttbar_xsec = 831.76;

std::map<string, double> xsecs = {
{"Data13TeV_SingleElectron2016B_PromptRecoV2" , 1.0},
{"Data13TeV_SingleMuon2016B_PromptRecoV2" , 1.0},
{"Data13TeV_SingleElectron2016D_PromptRecoV2_", 1},
{"Data13TeV_SingleMuon2016D_PromptRecoV2_", 1},

{ "MC2016_noHLT_TTJets_powheg_scaleup_elelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elmubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_elqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_eltaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muelbar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mumubar"         , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_muqbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_mutaubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qelbar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qmubar"          , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qqbar"           , 831.76 * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_qtaubar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauelbar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_taumubar"        , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tauqbar"         , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaleup_tautaubar"       , 831.76 * W_lep_br2 / 2},

{ "MC2016_noHLT_TTJets_powheg_scaledown_elelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elmubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_elqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_eltaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muelbar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mumubar"       , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_muqbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_mutaubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qelbar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qmubar"        , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qqbar"         , 831.76 * W_qar_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_qtaubar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauelbar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_taumubar"      , 831.76 * W_lep_br2 / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tauqbar"       , 831.76 * W_qar_br*W_lep_br / 2},
{ "MC2016_noHLT_TTJets_powheg_scaledown_tautaubar"     , 831.76 * W_lep_br2 / 2},

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
}

#endif // DTAG_XSECS_H_

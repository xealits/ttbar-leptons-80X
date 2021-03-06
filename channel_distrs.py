import logging
import argparse
from os import path
from collections import OrderedDict

# TODO: somehow move it to __main__
lep, channel_lep_id = 'el', 11 # or 'mu', 13 for mutau
lep, channel_lep_id = 'mu', 13






tt_el_processes = {
  'tt_eltauh': (["MC2016_Summer16_TTJets_powheg.root"],
                            ["((abs(gen_t_w_decay_id) == 11 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 11 && abs(gen_t_w_decay_id) > 15*15))"], 0),
  'tt_taueltauh': (["MC2016_Summer16_TTJets_powheg.root"],
                               ["((abs(gen_t_w_decay_id) == 15*11 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*11 && abs(gen_t_w_decay_id) > 15*15))"], 1),
  'tt_eljets': (["MC2016_Summer16_TTJets_powheg.root"],
                            [" (abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) "], 1),
  'tt_other':   (["MC2016_Summer16_TTJets_powheg.root"],
                 [" !(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && !((abs(gen_t_w_decay_id) == 15*11 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*11 && abs(gen_t_w_decay_id) > 15*15)) && !((abs(gen_t_w_decay_id) == 11 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 11 && abs(gen_t_w_decay_id) > 15*15))"], 1)
}

tt_mu_processes = {
  'tt_mutauh': (["MC2016_Summer16_TTJets_powheg.root"],
                            ["((abs(gen_t_w_decay_id) == 13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 13 && abs(gen_t_w_decay_id) > 15*15))"], 0),
  'tt_taumutauh': (["MC2016_Summer16_TTJets_powheg.root"],
                               ["((abs(gen_t_w_decay_id) == 15*13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*13 && abs(gen_t_w_decay_id) > 15*15))"], 1),
  'tt_mujets': (["MC2016_Summer16_TTJets_powheg.root"],
                            [" (abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13) "], 1),
  'tt_other':   (["MC2016_Summer16_TTJets_powheg.root"],
                 [" !(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13) && !((abs(gen_t_w_decay_id) == 15*13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*13 && abs(gen_t_w_decay_id) > 15*15)) && !((abs(gen_t_w_decay_id) == 13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 13 && abs(gen_t_w_decay_id) > 15*15))"], 1)
}

tt_elmu_processes = {'tt_elmu': (["MC2016_Summer16_TTJets_powheg.root"],
                                 ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*13)"], 0),
                     'tt_taueltaumu': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*13*15*15)"], 1),
                     'tt_other': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["!(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*13*15*15)"], 1) }

tt_elel_processes = {'tt_elel': (["MC2016_Summer16_TTJets_powheg.root"],
                                 ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*11)"], 0),
                     'tt_taueltauel': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*11*15*15)"], 1),
                     'tt_other': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["!(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*11 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11*11*15*15)"], 1) }

tt_mumu_processes = {'tt_mumu': (["MC2016_Summer16_TTJets_powheg.root"],
                                 ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*13)"], 0),
                     'tt_taumutaumu': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*13*15*15)"], 1),
                     'tt_other': (["MC2016_Summer16_TTJets_powheg.root"],
                                  ["!(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13*13*15*15)"], 1) }

tt_mujets_processes = {'tt': (["MC2016_Summer16_TTJets_powheg.root"], [''], 0)}


st_tW_el_processes = {
  'st_tW_eltauh':    (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == 11 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 11 && abs(gen_wdecays_IDs[0]) > 15*15))"], 2),
  'st_tW_taueltauh': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == 15*11 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*11 && abs(gen_wdecays_IDs[0]) > 15*15))"], 2),
  'st_tW_eljets':  (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                ["(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 11)"], 2),
  'st_tW_other':   (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                    ["!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 11) && !(abs(gen_wdecays_IDs[1]) == 11 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 11 && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*11 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*11 && abs(gen_wdecays_IDs[1]) > 15*15)"], 2)
}

st_tW_mu_processes = {
  'st_tW_mutauh':    (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15))"], 2),
  'st_tW_taumutauh': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15))"], 2),
  'st_tW_mujets':  (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                ["(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13)"], 2),
  'st_tW_other':   (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                    ["!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13) && !(abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15)"], 2)
}

st_tW_elmu_processes = {'st_tW_elmu': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13*11)'], 2),
                        'st_tW_taueltaumu': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*13*11)'], 2),
                        'st_tW_other': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*13*11 || abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13*11)'], 2) }

st_tW_mumu_processes = {'st_tW_mumu': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13*13)'], 2),
                        'st_tW_taumutaumu': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*13*13)'], 2),
                        'st_tW_other': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*13*13 || abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13*13)'], 2) }

st_tW_elel_processes = {'st_tW_elel': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 11*11)'], 2),
                        'st_tW_taueltauel': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*11*11)'], 2),
                        'st_tW_other': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                                       ['!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 15*15*11*11 || abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 11*11)'], 2) }

st_tW_mujets_processes = {'st_tW': (["MC2016_Summer16_SingleT_tW_5FS_powheg.root", "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"], [''], 2)}

dy_el_processes = {
  'dy_taueltauh':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && ((abs(gen_pythia8_prompt_leptons_IDs[1]) == 11*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 11*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))",
              "gen_N_zdecays > 0 && ((abs(gen_zdecays_IDs[1]) == 11*15 && abs(gen_zdecays_IDs[0]) > 15*15) || (abs(gen_zdecays_IDs[0]) == 11*15 && abs(gen_zdecays_IDs[1]) > 15*15))"], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !((abs(gen_pythia8_prompt_leptons_IDs[1]) == 11*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 11*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))",
              "gen_N_zdecays > 0 && !((abs(gen_zdecays_IDs[1]) == 11*15 && abs(gen_zdecays_IDs[0]) > 15*15) || (abs(gen_zdecays_IDs[0]) == 11*15 && abs(gen_zdecays_IDs[1]) > 15*15))"], 5)
}

dy_mu_processes = {
  'dy_taumutauh':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && ((abs(gen_pythia8_prompt_leptons_IDs[1]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))",
              "gen_N_zdecays > 0 && ((abs(gen_zdecays_IDs[1]) == 13*15 && abs(gen_zdecays_IDs[0]) > 15*15) || (abs(gen_zdecays_IDs[0]) == 13*15 && abs(gen_zdecays_IDs[1]) > 15*15))"], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !((abs(gen_pythia8_prompt_leptons_IDs[1]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))",
              "gen_N_zdecays > 0 && !((abs(gen_zdecays_IDs[1]) == 13*15 && abs(gen_zdecays_IDs[0]) > 15*15) || (abs(gen_zdecays_IDs[0]) == 13*15 && abs(gen_zdecays_IDs[1]) > 15*15))"], 5)
}

dy_elmu_processes = {
  'dy_taumutauel':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && (abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*11*15*15)",
              "gen_N_zdecays > 0 && (abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*13*15*15)"], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !(abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*11*15*15)",
              "gen_N_zdecays > 0 && !(abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*13*15*15)"], 5)
}

dy_mumu_processes = {
  'dy_mumu': (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
              ["gen_N_zdecays < 1 && (abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*13)",
               "gen_N_zdecays > 0 && (abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 13*13)"], 5),
  'dy_taumutaumu': (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && (abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*13*15*15)",
              "gen_N_zdecays > 0 && (abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 13*13*15*15)"], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !(abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*13*15*15 || abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 13*13)",
              "gen_N_zdecays > 0 && !(abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 13*13*15*15 || abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 13*13)"], 5)
}

dy_elel_processes = {
  'dy_elel': (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
              ["gen_N_zdecays < 1 && (abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 11*11)",
               "gen_N_zdecays > 0 && (abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*11)"], 5),
  'dy_taueltauel': (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && (abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 11*11*15*15)",
              "gen_N_zdecays > 0 && (abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*11*15*15)"], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !(abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 11*11*15*15 || abs(gen_pythia8_prompt_leptons_IDs[1] * gen_pythia8_prompt_leptons_IDs[0]) == 11*11)",
              "gen_N_zdecays > 0 && !(abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*11*15*15 || abs(gen_zdecays_IDs[1]*gen_zdecays_IDs[0]) == 11*11)"], 5)
}

dy_mujets_processes = {'dy': (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"], [''], 5)}





all_processes = {
  #'data': ([data_files[lep]], [""], 0),

  'st_schan':         (["MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo.root"],      [""], 3),
  'st_tchan_antitop': (["MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg.root"], [""], 3),
  'st_tchan_top':     (["MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg.root"],     [""], 3),

  #'wjets' : ([#'MC2016_Summer16_W0Jets_amcatnlo.root',
  #            'MC2016_Summer16_WJets_madgraph.root',
  #            'MC2016_Summer16_W1Jets_madgraph.root',
  #            'MC2016_Summer16_W2Jets_madgraph.root',
  #            'MC2016_Summer16_W3Jets_madgraph.root',
  #            'MC2016_Summer16_W4Jets_madgraph.root'], [""], 4),

  'wjets'  : (['MC2016_Summer16_WJets_madgraph.root'], [""], 4),
  'w1jets' : (['MC2016_Summer16_W1Jets_madgraph.root'], [""], 4),
  'w2jets' : (['MC2016_Summer16_W2Jets_madgraph.root'], [""], 4),
  'w3jets' : (['MC2016_Summer16_W3Jets_madgraph.root'], [""], 4),
  'w4jets' : (['MC2016_Summer16_W4Jets_madgraph.root'], [""], 4),

  'qcd': (
    ["MC2016_Summer16_QCD_HT-2000-Inf.root",
     "MC2016_Summer16_QCD_HT-1500-2000.root",
     "MC2016_Summer16_QCD_HT-1000-1500.root",
     "MC2016_Summer16_QCD_HT-700-1000.root",
     "MC2016_Summer16_QCD_HT-500-700.root",
     "MC2016_Summer16_QCD_HT-300-500.root",
     "MC2016_Summer16_QCD_HT-200-300.root",
     "MC2016_Summer16_QCD_HT-100-200.root"], [''], 7
  )
}

#event_selection += "HLT_mu == 1 && leps_ID == 13" if channel_lep_id == 13 else "HLT_el == 1 && abs(leps_ID) == 11"
#event_selection = "nbjets > 0 && tau_IDlev[0] > 2 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) > 0 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) < 250 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4"

data_files = {
    'el': ['Data13TeV_SingleElectron2016.root'],
    'mu': ['Data13TeV_SingleMuon2016.root',
           'Data13TeV_SingleMuon2016B_03Feb2017_ver2.root',
           'Data13TeV_SingleMuon2016C_03Feb2017_v1.root',
           'Data13TeV_SingleMuon2016D_03Feb2017_v1.root',
           'Data13TeV_SingleMuon2016E_03Feb2017_v1.root',
           'Data13TeV_SingleMuon2016F_03Feb2017_v1.root',
           'Data13TeV_SingleMuon2016H_03Feb2017_ver2.root',
           'Data13TeV_SingleMuon2016H_03Feb2017_ver3.root']
}

selection_specific = {'el': {'common': 'HLT_el && abs(leps_ID) == 11 && lep_matched_HLT[0] && met_corrected.pt() > 40 && nbjets > 0 && tau_IDlev[0] > 2 && tau_id[0]*lep_id[0] < 0 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4',
                             'processes': [{'data': (data_files['el'], [''], 0)}, tt_el_processes, dy_el_processes, st_tW_el_processes]},
                      'mu': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && met_corrected.pt() > 40 && nbjets > 0 && tau_IDlev[0] > 2 && tau_id[0]*lep_id[0] < 0 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4',
                             'processes': [{'data': (data_files['mu'], [''], 0)},     tt_mu_processes, dy_mu_processes, st_tW_mu_processes]},
                  'mujets': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && tau_IDlev[0] > 0 && transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',
                  #Float_t transverse_mass(Float_t v1_pt, Float_t v2_pt, Float_t v1_phi, Float_t v2_phi)
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mujets_processes, st_tW_mujets_processes, dy_mujets_processes]},
                  'mujets_b': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && nbjets > 0 && transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',
                  #Float_t transverse_mass(Float_t v1_pt, Float_t v2_pt, Float_t v1_phi, Float_t v2_phi)
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mujets_processes, st_tW_mujets_processes, dy_mujets_processes]},
               # HLT + match + tau ID medium + tau pt > 30 + 0 b tags + mT < 40 + 45 < (mu + tau).M < 85
               # HLT + match + tau ID medium + tau pt > 30 + 0 b tags + mT < 40 + 45 < (mu + tau).M < 85
               'taumutauh': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && nbjets == 0 && tau_IDlev[0] > 2 && tau_p4[0].pt() > 30 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) < 40 && divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 45 && divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 85',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mu_processes, st_tW_mujets_processes, dy_mu_processes]},
               'taumutauh_antiMt': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && nbjets == 0 && tau_IDlev[0] > 2 && tau_p4[0].pt() > 30 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) > 50 && (divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 45 || divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 85)',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mu_processes, st_tW_mujets_processes, dy_mu_processes]},
               # to be sure lets' add a requirement that some tau candidates are recorded,
               # even though there is a record scheme tauID which records only these events
               'taumutauh_antiMt_pretau': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && tau_IDlev[0] > -1 && nbjets == 0 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) > 50 && (divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 45 || divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 85)',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mu_processes, st_tW_mujets_processes, dy_mu_processes]},
               'taumutauh_antiMt_pretau_allb': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && tau_IDlev[0] > -1 && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) > 50 && (divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 45 || divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 85)',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mu_processes, st_tW_mujets_processes, dy_mu_processes]},




              'elmu_hltmu': {'common': 'HLT_mu && leps_ID == -13*11 && ((abs(lep_id[0]) == 13 && lep_matched_HLT[0]) || (abs(lep_id[1]) == 13 && lep_matched_HLT[1])) && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && lep_p4[1].pt() > 30 && fabs(lep_p4[1].eta()) < 2.4',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_elmu_processes, st_tW_elmu_processes, dy_elmu_processes]},
              'elmu_hltel': {'common': 'HLT_el && leps_ID == -13*11 && ((abs(lep_id[0]) == 11 && lep_matched_HLT[0]) || (abs(lep_id[1]) == 11 && lep_matched_HLT[1])) && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && lep_p4[1].pt() > 30 && fabs(lep_p4[1].eta()) < 2.4',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_elmu_processes, st_tW_elmu_processes, dy_elmu_processes]},
            'elmu_hltelmu': {'common': 'HLT_el && HLT_mu && leps_ID == -13*11 && lep_matched_HLT[0] && lep_matched_HLT[1] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && lep_p4[1].pt() > 30 && fabs(lep_p4[1].eta()) < 2.4',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_elmu_processes, st_tW_elmu_processes, dy_elmu_processes]},

                    'elel': {'common': 'HLT_el && abs(leps_ID) == 11*11 && ((abs(lep_id[0]) == 11 && lep_matched_HLT[0]) || (abs(lep_id[1]) == 11 && lep_matched_HLT[1])) && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && lep_p4[1].pt() > 30 && fabs(lep_p4[1].eta()) < 2.4',
                             'processes': [{'data': (data_files['el'], [''], 0)}, tt_elel_processes, st_tW_elel_processes, dy_elel_processes]},
                    'mumu': {'common': 'HLT_mu && abs(leps_ID) == 13*13 && ((abs(lep_id[0]) == 13 && lep_matched_HLT[0]) || (abs(lep_id[1]) == 13 && lep_matched_HLT[1])) && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 && lep_p4[1].pt() > 30 && fabs(lep_p4[1].eta()) < 2.4',
                             'processes': [{'data': (data_files['mu'], [''], 0)}, tt_mumu_processes, st_tW_mumu_processes, dy_mumu_processes]},

                     }


processes_histos = dict()

# info on MC cross-sections
# should move it to a separate module?
W_lep_br = 0.108;
W_qar_br = 0.676;

W_lep_br2 = W_lep_br*W_lep_br;
W_qar_br2 = W_qar_br*W_qar_br;

br_tau_electron = 0.1785;
br_tau_muon     = 0.1736;
br_tau_lepton   = br_tau_electron + br_tau_muon;
br_tau_hadronic = 1 - br_tau_lepton;

ttbar_xsec = 831.76;

wjets_normalization = 1.3 / 2.078 #1 / 2.078

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "Produce ROOT file with channels/processes/[NOMINAL,SYSTEMATIC] distributions for HiggsCombine fit from the Ntupler TTrees.",
        epilog = "Example:\n$ python channel_distrs.py directory_name"
        )

    parser.add_argument("ntupler_directory", help="the name of the directory with the Ntupler outputs")
    parser.add_argument("channel", help="channel to select: el, mu, elmu, mumu, mujets")
    parser.add_argument("--lumi",  help="the data luminosity", default=1, type=float)
    parser.add_argument("--debug", help="set DEBUG level in logging output", action = "store_true")
    parser.add_argument("--small", help="run on 1 small (DY) dataset for test", action = "store_true")
    parser.add_argument("--nup",   help="run with NUP cut on WJets (small effect)", action = "store_true")
    parser.add_argument("--wjet-xsec", help="run with reduced WJets xsec", action = "store_true")
    parser.add_argument("--wjets-normalization", help="factor from control sample", default=1., type=float)
    parser.add_argument("--process", help="run only on given process (data, tt_elmu)")
    parser.add_argument("--draw", help="draw distr")

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)
    logging.debug("ntupler_dir, lumi, small_test = %s, %r, %r" % (args.ntupler_directory, args.lumi, args.small))

    wjets_normalization = args.wjets_normalization

    mc_factors = {
     "MC2016_Summer16_TTJets_powheg.root": (ttbar_xsec, 1.077e-05),  # (xsec, old_full_mc_factor) maybe should add Nevents (from mc datasets

     "MC2016_Summer16_SingleT_tW_5FS_powheg.root":    (35.6,5.12022e-06),
     "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root": (35.6,5.13479e-06),
     "MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo.root":      ( 10.11,1.011e-05),
     "MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg.root": ( 80.95,2.08575e-06),
     "MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg.root":     (136.02,2.02288e-06),

     "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root":      (18610  ,0.000601864), # same
     "MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root":     ( 5765.4,0.000117316), # 6025.2
     "MC2016_Summer16_DYJetsToTauTau_50toInf_amcatnlo.root": ( 1867  ,6.95408e-05),

     "MC2016_Summer16_WJets_madgraph.root" : (61526.7 * wjets_normalization, 0.0019671),
     #"MC2016_Summer16_WJets_madgraph.root" : (61526.7 - 9493 - 3120 - 942.3 - 524.2, 0.0019671),
     #"MC2016_Summer16_WJets_madgraph.root" : (50690 - 9493 - 3120 - 942.3 - 524.2, 0.0019671),
     "MC2016_Summer16_WJets_amcatnlo.root" : (61526.7 * wjets_normalization, 0.0019671),
     "MC2016_Summer16_W0Jets_amcatnlo.root": ((61526.7 - 9493 - 3120 - 942.3 - 524.2) * wjets_normalization, 0.0019671),
     "MC2016_Summer16_W1Jets_madgraph.root": (9493  * wjets_normalization, 0.000209249),
     "MC2016_Summer16_W2Jets_madgraph.root": (3120  * wjets_normalization, 0.000104423),
     "MC2016_Summer16_W3Jets_madgraph.root": (942.3 * wjets_normalization, 4.75954e-05),
     "MC2016_Summer16_W4Jets_madgraph.root": (524.2 * wjets_normalization, 5.71611e-05),

     "MC2016_Summer16_QCD_HT-2000-Inf.root":  (   25.25,1.2678e-05 ),
     "MC2016_Summer16_QCD_HT-1500-2000.root": (  120.4 ,3.03212e-05),
     "MC2016_Summer16_QCD_HT-1000-1500.root": ( 1206   ,0.000252984),
     "MC2016_Summer16_QCD_HT-700-1000.root":  ( 6802   ,0.00043521 ),
     "MC2016_Summer16_QCD_HT-500-700.root":   (31630   ,0.0016709  ),
     "MC2016_Summer16_QCD_HT-300-500.root":   (351300  ,0.0206212  ),
     "MC2016_Summer16_QCD_HT-200-300.root":   (1717000 ,0.0917082  ),
     "MC2016_Summer16_QCD_HT-100-200.root":   (27540000,0.34133    )
    }

    if args.wjet_xsec:
        wjets_xsec, old_wjets_mc_factor = mc_factors["MC2016_Summer16_WJets_madgraph.root"]
        the_njets_xsec = (9493 + 3120 + 942.3 + 524.2) * wjets_normalization
        #mc_factors["MC2016_Summer16_WJets_madgraph.root"] = (wjets_xsec - the_njets_xsec, old_wjets_mc_factor)
        new_wjets_xsec = 50690. * wjets_normalization
        mc_factors["MC2016_Summer16_WJets_madgraph.root"] = (new_wjets_xsec, old_wjets_mc_factor)
    logging.info("set WJets xsec = %f" % mc_factors["MC2016_Summer16_WJets_madgraph.root"][0])

    for d in selection_specific[args.channel]['processes']:
        all_processes.update(d)

    if args.small:
        test_process = {'el': 'dy_taueltauh', 'elel': 'dy_elel', 'mu': 'dy_taumutauh', 'mujets': 'dy', 'elmu': 'dy_taumutauel', 'mumu': 'dy_mumu', 'taumutauh': 'dy_taumutauh'}[args.channel]
        small_test_process = {test_process: all_processes[test_process]}
        small_test_process.update(dy_elmu_processes)
        processes = small_test_process
    else:
        processes = all_processes

    if args.process:
        processes = {args.process: all_processes[args.process]}
        logging.info("running on %s" % repr(processes))

    # ROOT is heavy, init it after comline
    from ROOT import TFile, TTree, TH1D, gROOT, gSystem, TCanvas
    gROOT.SetBatch(True)

    gROOT.ProcessLine(".L pu_weight.C+")
    gSystem.Load("libUserCodettbar-leptons-80X.so")

    #draw_distr = "tau_flightLengthSignificance[0]>>h(100,-5,50)"
    #event_selection = "HLT_mu == 1 && leps_ID == 13 && tau_IDlev[0] > 1 && tau_flightLengthSignificance[0] > -5 && tau_flightLengthSignificance[0] < 50"

    #draw_distr = "lep_p4[0].pt()>>h(25,0,250)"
    #event_selection = "HLT_mu == 1 && leps_ID == 13" if channel_lep_id == 13 else "HLT_el == 1 && leps_ID == 11" + " && tau_IDlev[0] > 1 && lep_p4[0].pt() > 0 && lep_p4[0].pt() < 250"

    if args.draw:
        draw_distr = args.draw
        event_selection = ''
    else:
        draw_distr = "sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi())))>>h(25,0,250)"
        event_selection = "sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) > 0 && sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) < 250"
        event_selection += ' && '
    event_selection += selection_specific[args.channel]['common']
    #"HLT_mu == 1 && leps_ID == 13" if channel_lep_id == 13 else "HLT_el == 1 && abs(leps_ID) == 11"

    logging.debug("Draw    " + draw_distr)
    logging.debug("Select  " + event_selection)

    out_file = TFile("test_channels_output_%s.root" % args.channel, "recreate")
    channel1 = out_file.mkdir("channel1")

    # name: (parameters)
    # parameters = PU, 
    class NominalParameters:
        known_parameters = 'PU', 'TOP_PT'
        def __init__(self, PU=0, TOP_PT=0):
            '''__init__(self, PU=0, TOP_PT=0)

            PU      0 = nominal, 1 = up, -1 = down
            TOP_PT  0 = no weight, 1 = reweight (2nd Run parameters)
            '''
            args = {k:v for k, v in locals().items() if k is not 'self'}
            # save the input parameters in dict for unpacking
            self.parameters_dict = OrderedDict()
            for p in self.known_parameters:
                self.parameters_dict[p] = args.pop(p)
            if args: # still some input left
                print 'WARNING: these parameters are not known:', args

        # parameters
        def __getattr__(self, param_name):
            if param_name in self.parameters_dict:
               return self.parameters_dict[param_name]
            else:
               raise AttributeError("Parameter %s doesn't exist" % param_name)
        def keys(self):
           return self.parameters_dict.keys()
        # for tuple-like unpacking
        def __getitem__(self, index):
            return self.parameters_dict[index] if index in self.parameters_dict else self.parameters_dict[self.known_parameters[index]]

    systematics = {'NOMINAL': NominalParameters(), 'PU_UP': NominalParameters(PU=1), 'PU_DOWN': NominalParameters(PU=-1)}

    n_events_control, n_weight_control = {}, {}

    for name, (filenames, proc_defs, proc_color) in processes.items():
        logging.info(name)
        chan_proc_dir = channel1.mkdir(name)
        chan_proc_dir.cd()

        # for a process (= nickname) go over all given files
        # and extract this process according to given definitions
        # in case of DY there are 2 definitions, thus extract both ways
        # add up all the distributions for this process
        for filename in filenames:
          if path.exists(args.ntupler_directory + '/' + filename):
              logging.debug(filename)
              # DY processes have 2 definitions, since there are 2 ways to split MC
              for proc_def in proc_defs:
                  f = TFile(args.ntupler_directory + '/' + filename, 'read')
                  t = f.Get('ntupler/reduced_ttree')
                  subchan_selection = event_selection + (' && ' + proc_def if proc_def else '')
                  if 'WJe' in filename and args.nup:
                      logging.debug('W0Jets NUP cut')
                      subchan_selection += ' && NUP_gen == 5'

                  if 'data' not in name: # the weights for MC
                      # the aMCatNLO weight
                      logging.debug('adding MC weights to selection')
                      weight  = '(pu_weight(nvtx_gen, {}))* '.format(0) # 0 for nominal PU weights
                      if 'amcatnlo' in filename:
                          logging.debug('aMCatNLO -1 weights')
                          weight += '(aMCatNLO_weight > 0 ? 1 : -1)* ' 
                      if 'TT' in filename:
                          logging.debug('top pT')
                          weight += '(ttbar_pT_SF(gen_t_w_decay_id, gen_tb_w_decay_id))* '
                  else:
                      weight = ''

                  selection = weight + '({})'.format(subchan_selection)
                  logging.debug("draw_distr selection = {} {}".format(draw_distr, selection))
                  nentries_processed = t.Draw(draw_distr, selection)
                  h = t.GetHistogram().Clone()
                  logging.debug("N processed, N entries, Integral = %d %d %d" % (nentries_processed, h.GetEntries(), h.Integral()))
                  h.SetLineColor = proc_color
                  h.SetName("NOMINAL")

                  events_h = f.Get('ntupler/events_counter')
                  n_events = events_h.GetBinContent(2)

                  n_events_control[name] = n_events_control.get(name, 0) + n_events

                  if filename in mc_factors: # scaling with the xsec*lumi/Nevents factors
                      xsec, old_factor = mc_factors[filename]

                      # new factor = xsec / N original weight
                      # but I shouldn't rely on weight of counters -- theprocessing may use different weights
                      # therefore N original weight = N orig events (from the counter) * h_weight / h_events
                      # in debug compare it with N weights from the counter
                      weight_h = f.Get('ntupler/weight_counter')
                      n_weight = weight_h.GetBinContent(2)
                      n_weight_control[name] = n_weight_control.get(name, 0) + n_weight

                      n_new_weight = n_events * h.Integral() / (h.GetEntries() if h.GetEntries() > 0 else 1)
                      if n_new_weight == 0: # division by zero fix
                          n_new_weight = 1
                      # the mc factor is its' xsec / N events before cuts
                      # it normalizes N events of MC -- multiply by lumi to get expected amount
                      logging.debug("old_weight, new_weight = {} {}".format(n_weight, n_new_weight))
                      logging.debug("old_factor, new_factor = {} {}".format(old_factor, xsec / n_new_weight))
                      logging.debug("Scale %f by %f * %f = %f" % (h.Integral(), xsec / n_new_weight, args.lumi, args.lumi * xsec / n_new_weight))
                      if n_new_weight > 0:
                          h.Scale(args.lumi * xsec / n_new_weight)
                      #h.Scale(args.lumi * old_factor)
                      logging.debug("Scaled to %f" % h.Integral())
                  elif 'MC' in filename:
                      logging.info(filename + ' is not scaled, since it is not in the array of MC factors')

                  #h.SetDirectory(0)
                  #chan_proc_dir.Append(h)
                  #h.Write()
                  if name not in processes_histos:
                      h.SetDirectory(chan_proc_dir)
                      processes_histos[name] = h
                      logging.debug("new histo   %f" % processes_histos[name].Integral())
                  else:
                      processes_histos[name].Add(h)
                  logging.debug("added histo %f" % processes_histos[name].Integral())

    logging.debug("writing the file")
    out_file.Write()

    """
    gROOT.Reset()
    c1 = TCanvas( 'c1', 'Example', 10, 10, 700, 700 )
    c1.cd()

    c1.Update()
    c1.SaveAs( "test_channel_distrs.png" )
    """

    print args.channel
    mc_sum = sum(histo.Integral() for n, histo in processes_histos.items() if 'data' not in n)
    data_sum = processes_histos['data'].Integral() if 'data' in processes_histos else 1
    for i, (n, histo) in enumerate(processes_histos.items()):
        if 'data' in n: continue
        print '%-20s  %10.2f  %10.2f : %10.2f  %12.2f' % (n, histo.Integral(), histo.Integral() / mc_sum, histo.Integral() / data_sum, n_weight_control.get(n, 0))

    print "mc_sum", mc_sum
    if 'data' in processes_histos:
        print 'data', processes_histos['data'].Integral()
        print 'ratio', processes_histos['data'].Integral() / mc_sum
        print 'ratio', mc_sum / data_sum

    if 'data' in n_events_control:
        print 'data N events', n_events_control['data']


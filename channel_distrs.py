
channels = {
  'tt_mutauh': ("MC2016_Summer16_TTJets_powheg.root", "((abs(gen_t_w_decay_id) == 13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 13 && abs(gen_t_w_decay_id) > 15*15))", 0),
  'tt_taumutauh': ("MC2016_Summer16_TTJets_powheg.root", "((abs(gen_t_w_decay_id) == 15*13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*13 && abs(gen_t_w_decay_id) > 15*15))", 1),
  'tt_mujets':  ("MC2016_Summer16_TTJets_powheg.root", " abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 ", 1),
  'tt_other':   ("MC2016_Summer16_TTJets_powheg.root", " !abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 && !((abs(gen_t_w_decay_id) == 15*13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*13 && abs(gen_t_w_decay_id) > 15*15)) && !((abs(gen_t_w_decay_id) == 13 && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 13 && abs(gen_t_w_decay_id) > 15*15))", 1),
  'st_tW_mutauh':    ("MC2016_Summer16_SingleT_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15)", 2),
  'st_tW_taumutauh': ("MC2016_Summer16_SingleT_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15)", 2),
  'st_tW_mujets':  ("MC2016_Summer16_SingleT_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13)", 2),
  'st_tW_other':   ("MC2016_Summer16_SingleT_tW_5FS_powheg.root", "!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13) && !(abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15)", 2),
  'stbar_tW_mutauh':    ("MC2016_Summer16_SingleTbar_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15)", 2),
  'stbar_tW_taumutauh': ("MC2016_Summer16_SingleTbar_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15)", 2),
  'stbar_tW_mujets':    ("MC2016_Summer16_SingleTbar_tW_5FS_powheg.root", "(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13)", 2),
  'stbar_tW_other':     ("MC2016_Summer16_SingleTbar_tW_5FS_powheg.root", "!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == 13) && !(abs(gen_wdecays_IDs[1]) == 13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 13 && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*13 && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*13 && abs(gen_wdecays_IDs[1]) > 15*15)", 2),

  'st_schan': ("MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo.root", "", 3),
  'st_tchan_antitop': ("MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg.root", "", 3),
  'st_tchan_top':     ("MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg.root", "", 3),

  'wjets' : ("MC2016_Summer16_W4Jets_madgraph.root", "", 4),

  'dy_taumutauh':  ("MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "(abs(gen_pythia8_prompt_leptons_IDs[1]) == 13*15 && gen_pythia8_prompt_leptons_IDs[0] > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 13*15 && gen_pythia8_prompt_leptons_IDs[1] > 15*15) || (abs(gen_zdecays_IDs[1]) == 13*15 && gen_zdecays_IDs[0] > 15*15) || (abs(gen_zdecays_IDs[0]) == 13*15 && gen_zdecays_IDs[1] > 15*15)", 5),
# for some reason this doesn't work:
#root [23] t->Draw("HLT_mu", "(abs(gen_pythia8_prompt_leptons_IDs[1]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == 13*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15) || (abs(gen_zdecays_IDs[1]) == 13*15 && abs(gen_zdecays_IDs[0]) > 15*15) || (abs(gen_zdecays_IDs[0]) == 13*15 && abs(gen_zdecays_IDs[1] > 15*15))")


  'dy_other'    :  ("MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "", 5),
}



for name, (filename, ch_def, ch_color) in channels.items():
    f = TFile(filename, 'read')
    t = f.Get('ntupler/reduced_ttree')
    t.Draw("tau_flightLengthSignificance[0]", "tau_IDlev[0] > 1 && tau_flightLengthSignificance[0] > -5 && tau_flightLengthSignificance[0] < 50" + (' && ' + ch_def if ch_def else ''))
    h = t.GetHistogram().Clone()
    h.SetDirectory(0)
    channel_histos[name] = h        



# >>> [type(channel_histos[i]) for i in channel_histos]
# [<class '__main__.TH1F'>, <class '__main__.TH1F'>, <class '__main__.TH1F'>]


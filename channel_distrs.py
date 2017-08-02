import logging
import argparse


lep, channel_lep_id = 'el', 11 # or 'mu', 13 for mutau

parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "Produce ROOT file with channels/processes/[NOMINAL,SYSTEMATIC] distributions for HiggsCombine fit from the Ntupler TTrees.",
    epilog = "Example:\n$ python channel_distrs.py directory_name"
    )

parser.add_argument("ntupler_directory", help="the name of the directory with the Ntupler outputs")
parser.add_argument("--lumi",  help="the data luminosity", default=1, type=float)
parser.add_argument("--debug", help="set DEBUG level in logging output", action = "store_true")
test_process = 'dy_tau{}tauh'.format(lep)
parser.add_argument("--small", help="run on 1 small (%s) dataset for test" % test_process, action = "store_true")

args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)

logging.debug("ntupler_dir, lumi, small_test = %s, %r, %r" % (args.ntupler_directory, args.lumi, args.small))

# ROOT is heavy, init it after comline
from ROOT import TFile, TTree, TH1D, gROOT, TCanvas
gROOT.SetBatch(True)


all_processes = {
  'data': (['Data13TeV_SingleElectron2016.root'], [""], 0),
  'tt_{}tauh'.format(lep): (["MC2016_Summer16_TTJets_powheg.root"],
                            ["((abs(gen_t_w_decay_id) == {lep} && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == {lep} && abs(gen_t_w_decay_id) > 15*15))".format(lep=channel_lep_id)], 0),
  'tt_tau{}tauh'.format(lep): (["MC2016_Summer16_TTJets_powheg.root"],
                               ["((abs(gen_t_w_decay_id) == 15*{lep} && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*{lep} && abs(gen_t_w_decay_id) > 15*15))".format(lep=channel_lep_id)], 1),
  'tt_{}jets'.format(lep): (["MC2016_Summer16_TTJets_powheg.root"],
                            [" (abs(gen_t_w_decay_id * gen_tb_w_decay_id) == {lep}) ".format(lep=channel_lep_id)], 1),
  'tt_other':   (["MC2016_Summer16_TTJets_powheg.root"],
                 [" !(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == {lep}) && !((abs(gen_t_w_decay_id) == 15*{lep} && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == 15*{lep} && abs(gen_t_w_decay_id) > 15*15)) && !((abs(gen_t_w_decay_id) == {lep} && abs(gen_tb_w_decay_id) > 15*15) || (abs(gen_tb_w_decay_id) == {lep} && abs(gen_t_w_decay_id) > 15*15))".format(lep=channel_lep_id)], 1),

  'st_tW_{}tauh'.format(lep):    (["MC2016_Summer16_SingleT_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == {lep} && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == {lep} && abs(gen_wdecays_IDs[0]) > 15*15))".format(lep=channel_lep_id)], 2),
  'st_tW_tau{}tauh'.format(lep): (["MC2016_Summer16_SingleT_tW_5FS_powheg.root"],
                                  ["((abs(gen_wdecays_IDs[0]) == 15*{lep} && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*{lep} && abs(gen_wdecays_IDs[0]) > 15*15))".format(lep=channel_lep_id)], 2),
  'st_tW_{}jets'.format(lep):  (["MC2016_Summer16_SingleT_tW_5FS_powheg.root"],
                                ["(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == {lep})".format(lep=channel_lep_id)], 2),
  'st_tW_other':   (["MC2016_Summer16_SingleT_tW_5FS_powheg.root"],
                    ["!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == {lep}) && !(abs(gen_wdecays_IDs[1]) == {lep} && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == {lep} && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*{lep} && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*{lep} && abs(gen_wdecays_IDs[1]) > 15*15)".format(lep=channel_lep_id)], 2),
  'stbar_tW_{}tauh'.format(lep):    (["MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                         ["((abs(gen_wdecays_IDs[0]) == {lep} && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == {lep} && abs(gen_wdecays_IDs[0]) > 15*15))".format(lep=channel_lep_id)], 2),
  'stbar_tW_tau{}tauh'.format(lep): (["MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                         ["((abs(gen_wdecays_IDs[0]) == 15*{lep} && abs(gen_wdecays_IDs[1]) > 15*15) || (abs(gen_wdecays_IDs[1]) == 15*{lep} && abs(gen_wdecays_IDs[0]) > 15*15))".format(lep=channel_lep_id)], 2),
  'stbar_tW_{}jets'.format(lep):    (["MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                         ["(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == {lep})".format(lep=channel_lep_id)], 2),
  'stbar_tW_other':     (["MC2016_Summer16_SingleTbar_tW_5FS_powheg.root"],
                         ["!(abs(gen_wdecays_IDs[0] * gen_wdecays_IDs[1]) == {lep}) && !(abs(gen_wdecays_IDs[1]) == {lep} && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == {lep} && abs(gen_wdecays_IDs[1]) > 15*15) && !(abs(gen_wdecays_IDs[1]) == 15*{lep} && abs(gen_wdecays_IDs[0]) > 15*15) && !(abs(gen_wdecays_IDs[0]) == 15*{lep} && abs(gen_wdecays_IDs[1]) > 15*15)".format(lep=channel_lep_id)], 2),

  'st_schan':         (["MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo.root"],      [""], 3),
  'st_tchan_antitop': (["MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg.root"], [""], 3),
  'st_tchan_top':     (["MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg.root"],     [""], 3),

  'wjets' : (['MC2016_Summer16_W0Jets_amcatnlo.root',
              'MC2016_Summer16_W1Jets_madgraph.root',
              'MC2016_Summer16_W2Jets_madgraph.root',
              'MC2016_Summer16_W3Jets_madgraph.root',
              'MC2016_Summer16_W4Jets_madgraph.root'], [""], 4),

  'dy_tau{}tauh'.format(lep):  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
                    ["gen_N_zdecays < 1 && ((abs(gen_pythia8_prompt_leptons_IDs[1]) == {lep}*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == {lep}*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))".format(lep=channel_lep_id),
              "((abs(gen_zdecays_IDs[1]) == {lep}*15 && gen_zdecays_IDs[0] > 15*15) || (abs(gen_zdecays_IDs[0]) == {lep}*15 && gen_zdecays_IDs[1] > 15*15))".format(lep=channel_lep_id)], 5),
  'dy_other':  (["MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root", "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"],
             ["gen_N_zdecays < 1 && !((abs(gen_pythia8_prompt_leptons_IDs[1]) == {lep}*15 && abs(gen_pythia8_prompt_leptons_IDs[0]) > 15*15) || (abs(gen_pythia8_prompt_leptons_IDs[0]) == {lep}*15 && abs(gen_pythia8_prompt_leptons_IDs[1]) > 15*15))".format(lep=channel_lep_id),
              " !((abs(gen_zdecays_IDs[1]) == {lep}*15 && gen_zdecays_IDs[0] > 15*15) || (abs(gen_zdecays_IDs[0]) == {lep}*15 && gen_zdecays_IDs[1] > 15*15))".format(lep=channel_lep_id)], 5)
}

small_test_process = {test_process: all_processes[test_process]}

processes = small_test_process if args.small else all_processes

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
 
mc_factors = {
 "MC2016_Summer16_TTJets_powheg.root": (ttbar_xsec, 1.077e-05),  # (xsec, old_full_mc_factor) maybe should add Nevents (from mc datasets

 "MC2016_Summer16_SingleT_tW_5FS_powheg.root":    (35.6,5.12022e-06),
 "MC2016_Summer16_SingleTbar_tW_5FS_powheg.root": (35.6,5.13479e-06),
 "MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo.root":      ( 10.11,1.011e-05),
 "MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg.root": ( 80.95,2.08575e-06),
 "MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg.root":     (136.02,2.02288e-06),

 "MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root":      (18610  ,0.000601864),
 "MC2016_Summer16_DYJetsToLL_50toInf_madgraph.root":     ( 5765.4,0.000117316),
 "MC2016_Summer16_DYJetsToTauTau_50toInf_amcatnlo.root": ( 1867  ,6.95408e-05),

 "MC2016_Summer16_WJets_amcatnlo.root" : (61526.7,0.0019671),
 "MC2016_Summer16_W0Jets_amcatnlo.root": (61526.7 - 9493 - 3120 - 942.3 - 524.2,0.0019671),
 "MC2016_Summer16_W1Jets_madgraph.root": (9493 ,0.000209249),
 "MC2016_Summer16_W2Jets_madgraph.root": (3120 ,0.000104423),
 "MC2016_Summer16_W3Jets_madgraph.root": (942.3,4.75954e-05),
 "MC2016_Summer16_W4Jets_madgraph.root": (524.2,5.71611e-05),

 "MC2016_Summer16_QCD_HT-2000-Inf.root":  (   25.25,1.2678e-05 ),
 "MC2016_Summer16_QCD_HT-1500-2000.root": (  120.4 ,3.03212e-05),
 "MC2016_Summer16_QCD_HT-1000-1500.root": ( 1206   ,0.000252984),
 "MC2016_Summer16_QCD_HT-700-1000.root":  ( 6802   ,0.00043521 ),
 "MC2016_Summer16_QCD_HT-500-700.root":   (31630   ,0.0016709  ),
 "MC2016_Summer16_QCD_HT-300-500.root":   (351300  ,0.0206212  ),
 "MC2016_Summer16_QCD_HT-200-300.root":   (1717000 ,0.0917082  ),
 "MC2016_Summer16_QCD_HT-100-200.root":   (27540000,0.34133    )

}


#draw_distr = "tau_flightLengthSignificance[0]>>h(100,-5,50)"
#event_selection = "HLT_mu == 1 && leps_ID == 13 && tau_IDlev[0] > 1 && tau_flightLengthSignificance[0] > -5 && tau_flightLengthSignificance[0] < 50"

draw_distr = "sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi())))>>h(25,0,250)"
event_selection = "HLT_mu == 1 && leps_ID == 13" if channel_lep_id == 13 else "HLT_el == 1 && leps_ID == 11" + " && tau_IDlev[0] > 1"

logging.debug("Draw    " + draw_distr)
logging.debug("Select  " + event_selection)

out_file = TFile("test_channels_output.root", "recreate")
channel1 = out_file.mkdir("channel1")

for name, (filenames, proc_defs, proc_color) in processes.items():
    logging.info(name)
    chan_proc_dir = channel1.mkdir(name)
    chan_proc_dir.cd()

    # for a process (= nickname) go over all given files
    # and extract this process according to given definitions
    # in case of DY there are 2 definitions, thus extract both ways
    # add up all the distributions for this process
    for filename in filenames:
        # DY processes have 2 definitions, since there are 2 ways to split MC
        for proc_def in proc_defs:
            f = TFile(args.ntupler_directory + '/' + filename, 'read')
            t = f.Get('ntupler/reduced_ttree')
            subchan_selection = event_selection + (' && ' + proc_def if proc_def else '')
            logging.debug("Select    " + subchan_selection)
            nentries_processed = t.Draw(draw_distr, subchan_selection)
            logging.debug("N entries = %d" % nentries_processed)
            h = t.GetHistogram().Clone()
            h.SetLineColor = proc_color
            h.SetName("NOMINAL")

            if filename in mc_factors: # scaling with the xsec*lumi/Nevents factors
                xsec, old_factor = mc_factors[filename]
                # the mc factor is its' xsec / N events before cuts
                # it normalizes N events of MC -- multiply by lumi to get expected amount
                logging.debug("Scale %f by %f * %f" % (h.Integral(), old_factor, args.lumi))
                h.Scale(old_factor*args.lumi)
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

mc_sum = 0
for i, (n, histo) in enumerate(processes_histos.items()):
    if n is 'data': continue
    print n, histo.Integral()
    mc_sum += histo.Integral()

print "mc_sum", mc_sum
print 'data', processes_histos['data'].Integral()


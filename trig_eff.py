import argparse, logging



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "extract the HLT efficiencies and save them",
    epilog = "Example:\n$ python trig_eff.py data_input.root data"
    )

parser.add_argument("ntuple_file", help="the name of the file with the ntuple")
parser.add_argument("file_nick",   help="the nick name of the file with the ntuple")
parser.add_argument("--debug", help="set DEBUG level in logging output", action = "store_true")
parser.add_argument("-p", "--pt-max", type=float, default=150, help="cut at pT of lepton")
parser.add_argument("-l", "--lep-id", type=float, default=11,  help="PDG ID of the selected: 11 (electron) or 13 (muon)")

args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)



from ROOT import TH1D, TFile, TTree, gROOT, gStyle
gROOT.SetBatch(True)
gStyle.SetOptStat(0)


f = TFile(args.ntuple_file)
t = f.Get("ntupler/reduced_ttree")

weight_selection = ''
if 'MC' in args.ntuple_file:
    gROOT.ProcessLine(".L pu_weight.C+")
    weight_selection = '(pu_weight(nvtx_gen, {}))* '.format(0) # 0 is nominal

draw_com = "lep_p4[0].pt()>>h(100, 0, 200)"

draw_selection = "(HLT_lepMonitor && lep_id[0]*lep_id[1] == -11*11 && lep_p4[0].pt() > 20 && lep_p4[0].pt() < {} && lep_p4[1].pt() > 20 && lep_p4[1].pt() < {} && (lep_matched_HLT[0] || lep_matched_HLT[1]))".format(args.pt_max, args.pt_max)
#draw_selection = "(HLT_lepMonitor && abs(lep_id[0]) == 11 && lep_p4[0].pt() > 20 && lep_p4[0].pt() < {} && lep_matched_HLT[0])".format(args.pt_max)
#draw_selection = "(HLT_lepMonitor && abs(lep_id[0]) == {} && lep_p4[0].pt() > 20 && lep_p4[0].pt() < {} {})"
draw_selection = "(HLT_lepMonitor && abs(leps_ID) == {} && lep_p4[0].pt() > 20 && lep_p4[0].pt() < {} {})"

t.Draw(draw_com, weight_selection + draw_selection.format(args.lep_id, args.pt_max, "&& lep_matched_HLT[0]"))

fout = TFile("trigEff_" + args.file_nick + '.root', "recreate")

h_ = t.GetHistogram()

h_match = h_.Clone()
h_match.SetDirectory(0)
h_match.SetName("trig_eff_" + args.file_nick)

t.Draw(draw_com, weight_selection + draw_selection.format(args.lep_id, args.pt_max, ''))

h_ = t.GetHistogram()
h_all = h_.Clone()
h_all.SetName("trig_eff_" + args.file_nick + "_all")
h_match.Divide(h_all)

h_match.SetDirectory(fout)

fout.Write()
fout.Close()


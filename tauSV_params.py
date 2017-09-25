import logging
from sys import argv

logging.basicConfig(level=logging.INFO)

logging.debug("loading ROOT...")
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, gStyle, TCanvas, gPad
from ROOT.TMath import Sqrt
gROOT.SetBatch(True)
#gStyle.SetOptStat(0)

filename_tt = 'test_full_loop_on_taus2_matchQuality_runMillion_tt'
filename_tt = 'test_full_loop_on_taus2_tt'
filename_wjets = 'test_full_loop_on_taus2_wjets_'
filename_qcd = 'test_full_loop_on_taus2_qcd_'
filename_data = 'test_full_loop_on_taus2_data_'

params_tt_1D = [
  'refit_flightLen_to_other_PV_tt',	
  "refit_flightLen_to_other_PV_large_tt",
  'refit_flightLen_small_tt',	
  'refit_flightLen_tt',	
  'refit_flightSign_tt',	
  'refit_flightLen_lt',	
  'refit_flightLen_lj',	
  'refit_flightSign_lt',	
  'refit_flightSign_lj',	
  'pat_flightLen_tt',	
  'pat_flightLen_lt',	
  'pat_flightLen_lj',	
  'pat_flightSign_tt',	
  'pat_flightSign_lt',	
  'pat_flightSign_lj',	
  'refit_flightSign_lt_Btag',	
  'refit_flightSign_lj_Btag',	
  'pat_flightSign_lt_Btag',	
  'pat_flightSign_lj_Btag',	
  'refit_flightSign_lt_noBtag',	
  'refit_flightSign_lj_noBtag',	
  'pat_flightSign_lt_noBtag',	
  'pat_flightSign_lj_noBtag',	
  'n_goodPV_tt',	
  'refit_PV_z_tt',
  'refit_flightErr_lt',
  'refit_flightErr_lj',
  'refit_flightErr_tt',
  'refit_PV_to_all_goodPV_small_tt',
  'refit_PV_to_all_goodPV_tt',
  'refit_PV_to_all_goodPV_large_tt',
  'refit_PV_to_other_goodPV_small_tt',
  'refit_PV_to_other_goodPV_tt',
  'refit_PV_to_other_goodPV_large_tt',
]

params_tt_2D = [
  'pat_flightLen_refit_flightLen_tt',	
  'pat_flightSign_refit_flightSign_tt',	
  'pat_flightLen_flightSign_tt',	
  'pat_flightLen_flightSign_lt',	
  'pat_flightLen_flightSign_lj',	
  'pat_bTag_flightSign_tt',	
  'pat_bTag_flightSign_lt',	
  'pat_bTag_flightSign_lj',	
  'refit_m1m2_lt',	
  'refit_m1m2_lj',	
  'refit_flightLen_flightSign_tt',
  'refit_flightLen_flightSign_lt',
  'refit_flightLen_flightSign_lj',
  'refit_bTag_flightSign_tt',
  'refit_bTag_flightSign_lt',
  'refit_bTag_flightSign_lj',
  'refit_flightLen_Energy_tt',
  'refit_flightLen_Energy_lt',
  'refit_flightLen_Energy_lt_large',
  'refit_flightLen_Energy_lj',
  'pat_flightLen_Energy_tt',
  'pat_flightLen_Energy_lt',
  'pat_flightLen_Energy_lj',
  'refit_PV_xy_tt',
]

params_wjets_1D = [
  'refit_flightLen_to_other_PV_wjets',	
  'refit_flightLen_to_other_PV_large_wjets',
  'refit_flightLen_small_wjets',	
  'refit_flightLen_wjets',	
  'refit_flightSign_wjets',	
  'pat_flightLen_wjets',	
  'pat_flightSign_wjets',	
  'n_goodPV_wjets',	
  'refit_PV_z_wjets',
  'refit_flightErr_wjets',
  'refit_PV_to_all_goodPV_small_wjets',
  'refit_PV_to_all_goodPV_wjets',
  'refit_PV_to_all_goodPV_large_wjets',
  'refit_PV_to_other_goodPV_small_wjets',
  'refit_PV_to_other_goodPV_wjets',
  'refit_PV_to_other_goodPV_large_wjets',
]

params_wjets_2D = [
  'pat_flightLen_refit_flightLen_wjets',	
  'pat_flightSign_refit_flightSign_wjets',	
  'pat_flightLen_flightSign_wjets',	
  'pat_bTag_flightSign_wjets',	
  'refit_flightLen_flightSign_wjets',
  'refit_bTag_flightSign_wjets',
  'refit_flightLen_Energy_wjets',
  'pat_flightLen_Energy_wjets',
  'refit_PV_xy_wjets',
  'refit_m1m2_wjets',
]

params_data_1D = [
  'refit_flightLen_to_other_PV_data',	
  'refit_flightLen_to_other_PV_large_data',
  'refit_flightLen_small_data',	
  'refit_flightLen_data',	
  'refit_flightSign_data',	
  'pat_flightLen_data',	
  'pat_flightSign_data',	
  'n_goodPV_data',	
  'refit_PV_z_data',
  'refit_flightErr_data',
  'refit_PV_to_all_goodPV_small_data',
  'refit_PV_to_all_goodPV_data',
  'refit_PV_to_all_goodPV_large_data',
  'refit_PV_to_other_goodPV_small_data',
  'refit_PV_to_other_goodPV_data',
  'refit_PV_to_other_goodPV_large_data',
]

params_data_2D = [
  'pat_flightLen_refit_flightLen_data',	
  'pat_flightSign_refit_flightSign_data',	
  'pat_flightLen_flightSign_data',	
  'pat_bTag_flightSign_data',	
  'refit_flightLen_flightSign_data',
  'refit_bTag_flightSign_data',
  'refit_flightLen_Energy_data',
  'pat_flightLen_Energy_data',
  'refit_PV_xy_data',
  'refit_m1m2_data',
]

params_qcd_1D = [
  'refit_flightLen_to_other_PV_qcd',	
  'refit_flightLen_to_other_PV_large_qcd',
  'refit_flightLen_small_qcd',	
  'refit_flightLen_qcd',	
  'refit_flightSign_qcd',	
  'pat_flightLen_qcd',	
  'pat_flightSign_qcd',	
  'n_goodPV_qcd',	
  'refit_PV_z_qcd',
  'refit_flightErr_qcd',
  'refit_PV_to_all_goodPV_small_qcd',
  'refit_PV_to_all_goodPV_qcd',
  'refit_PV_to_all_goodPV_large_qcd',
  'refit_PV_to_other_goodPV_small_qcd',
  'refit_PV_to_other_goodPV_qcd',
  'refit_PV_to_other_goodPV_large_qcd',
]

params_qcd_2D = [
  'pat_flightLen_refit_flightLen_qcd',	
  'pat_flightSign_refit_flightSign_qcd',	
  'pat_flightLen_flightSign_qcd',	
  'pat_bTag_flightSign_qcd',	
  'refit_flightLen_flightSign_qcd',
  'refit_bTag_flightSign_qcd',
  'refit_flightLen_Energy_qcd',
  'pat_flightLen_Energy_qcd',
  'refit_PV_xy_qcd',
  'refit_m1m2_qcd',
]

params_1D = params_tt_1D
params_2D = params_tt_2D
filename = filename_tt

if '-w' in argv:
    params_1D = params_wjets_1D
    params_2D = params_wjets_2D
    filename = filename_wjets

if '-q' in argv:
    params_1D = params_qcd_1D
    params_2D = params_qcd_2D
    filename = filename_qcd

if '-d' in argv:
    params_1D = params_data_1D
    params_2D = params_data_2D
    filename = filename_data

in_file = TFile(filename + '.root', "read")

c1 = TCanvas("c", "canvas", 800, 800) #

for p in params_1D:
    #logging.debug('drawing ' + p)
    #if p == "refit_flightLen_tt": c1.SetLogy()
    in_file.Get(p).Draw()
    c1.SaveAs('plots/track-reco-investigation/v1/' + filename + '/tauSV_' + p + ".png")
    #if p == "refit_flightLen_tt": c1.SetLogy(0) # unset

for p in params_2D:
    in_file.Get(p).Draw("col")
    c1.SaveAs('plots/track-reco-investigation/v1/' + filename + '/tauSV_' + p + ".png")


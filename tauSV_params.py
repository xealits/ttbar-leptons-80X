import logging
from sys import argv

logging.basicConfig(level=logging.INFO)

logging.debug("loading ROOT...")
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, gStyle, TCanvas, gPad
from ROOT.TMath import Sqrt
gROOT.SetBatch(True)
#gStyle.SetOptStat(0)

filename = 'test_full_loop_on_taus2_matchQuality_runMillion_tt'
filename = 'test_full_loop_on_taus2_tt'
in_file = TFile(filename + '.root', "read")

params_1D = [
  'refit_flightLen_to_other_PV_tt',	
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
]

params_2D = [
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
  'refit_flightLen_Energy_lj',	
  'pat_flightLen_Energy_tt',	
  'pat_flightLen_Energy_lt',	
  'pat_flightLen_Energy_lj',	
  'refit_PV_xy_tt',	
]

c1 = TCanvas("c", "canvas", 800, 800) #

for p in params_1D:
    logging.debug('drawing ' + p)
    in_file.Get(p).Draw()
    c1.SaveAs('plots/track-reco-investigation/v1/' + filename + '/tauSV_' + p + ".png")

for p in params_2D:
    in_file.Get(p).Draw("col")
    c1.SaveAs('plots/track-reco-investigation/v1/' + filename + '/tauSV_' + p + ".png")


import logging
from sys import argv

logging.basicConfig(level=logging.DEBUG)

logging.debug("loading ROOT...")
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, gStyle, TCanvas
from ROOT.TMath import Sqrt
gROOT.SetBatch(True)
gStyle.SetOptStat(0)

filename = 'test_full_loop_on_taus2_matchQuality_runMillion_tt.root'
in_file = TFile(filename, "read")

h_refit_SV_cov00_tt = in_file.Get("refit_SV_cov00_tt")
h_refit_SV_cov01_tt = in_file.Get("refit_SV_cov01_tt")
h_refit_SV_cov02_tt = in_file.Get("refit_SV_cov02_tt")
h_refit_SV_cov10_tt = in_file.Get("refit_SV_cov10_tt")
h_refit_SV_cov11_tt = in_file.Get("refit_SV_cov11_tt")
h_refit_SV_cov12_tt = in_file.Get("refit_SV_cov12_tt")
h_refit_SV_cov20_tt = in_file.Get("refit_SV_cov20_tt")
h_refit_SV_cov21_tt = in_file.Get("refit_SV_cov21_tt")
h_refit_SV_cov22_tt = in_file.Get("refit_SV_cov22_tt")

h_refit_PV_cov00_tt = in_file.Get("refit_PV_cov00_tt")
h_refit_PV_cov01_tt = in_file.Get("refit_PV_cov01_tt")
h_refit_PV_cov02_tt = in_file.Get("refit_PV_cov02_tt")
h_refit_PV_cov10_tt = in_file.Get("refit_PV_cov10_tt")
h_refit_PV_cov11_tt = in_file.Get("refit_PV_cov11_tt")
h_refit_PV_cov12_tt = in_file.Get("refit_PV_cov12_tt")
h_refit_PV_cov20_tt = in_file.Get("refit_PV_cov20_tt")
h_refit_PV_cov21_tt = in_file.Get("refit_PV_cov21_tt")
h_refit_PV_cov22_tt = in_file.Get("refit_PV_cov22_tt")

logging.debug("plotting SV cov-s")

c1 = TCanvas("cSV", "SV canvas", 800, 800) #
c1.Clear()
c1.Divide(3,3)
c1.cd(0)
h_refit_SV_cov00_tt.Draw()
c1.cd(1)
h_refit_SV_cov01_tt.Draw()
c1.cd(2)
h_refit_SV_cov02_tt.Draw()
c1.cd(3)
h_refit_SV_cov10_tt.Draw()
c1.cd(4)
h_refit_SV_cov11_tt.Draw()
c1.cd(5)
h_refit_SV_cov12_tt.Draw()
c1.cd(6)
h_refit_SV_cov20_tt.Draw()
c1.cd(7)
h_refit_SV_cov21_tt.Draw()
c1.cd(8)
h_refit_SV_cov22_tt.Draw()

c1.SaveAs("test_full_loop_on_taus2_matchQuality_runMillion_tt_SV.png")

logging.debug("plotting PV cov-s")

c2 = TCanvas("cPV", "PV canvas", 800, 800) #
c2.Divide(3,3)
c2.cd(0)
h_refit_PV_cov00_tt.Draw()
c2.cd(1)
h_refit_PV_cov01_tt.Draw()
c2.cd(2)
h_refit_PV_cov02_tt.Draw()
c2.cd(3)
h_refit_PV_cov10_tt.Draw()
c2.cd(4)
h_refit_PV_cov11_tt.Draw()
c2.cd(5)
h_refit_PV_cov12_tt.Draw()
c2.cd(6)
h_refit_PV_cov20_tt.Draw()
c2.cd(7)
h_refit_PV_cov21_tt.Draw()
c2.cd(8)
h_refit_PV_cov22_tt.Draw()

c2.SaveAs("test_full_loop_on_taus2_matchQuality_runMillion_tt_PV.png")



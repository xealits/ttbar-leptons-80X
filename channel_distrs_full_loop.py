from array import array
from sys import argv
from os import environ
from collections import OrderedDict
import logging
import cProfile


logging.basicConfig(level=logging.DEBUG)

logging.info('importing ROOT')
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, TCanvas, TGraphAsymmErrors, TMath

pileup_ratio = array('d', [0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922,])
# 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

pileup_ratio_up = array('d', [0, 0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536,])
#0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

pileup_ratio_down = array('d', [0, 0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886,])
#0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])









logging.info("leptonic SFs")

muon_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/muon-effs/"
gSystem.ExpandPathName(muon_effs_dirname    )
#TString muon_effs_dirname = "analysis/muon-effs/";
logging.info("muon SFs from " + muon_effs_dirname)

muon_effs_tracking_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_tracking_more_BCDEF_fits.root" )
muon_effs_tracking_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_tracking_more_GH_fits.root" )
muon_effs_tracking_BCDEF_graph = muon_effs_tracking_BCDEF_file.Get("ratio_eff_aeta_dr030e030_corr")
muon_effs_tracking_GH_graph    = muon_effs_tracking_GH_file.Get("ratio_eff_aeta_dr030e030_corr")

print type(muon_effs_tracking_BCDEF_graph)

muon_effs_id_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_BCDEF.root" )
muon_effs_id_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_GH.root" )
muon_effs_id_BCDEF_histo = muon_effs_id_BCDEF_file.Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta").Get("pt_abseta_ratio")
muon_effs_id_GH_histo    = muon_effs_id_GH_file   .Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta").Get("pt_abseta_ratio")

print type(muon_effs_id_BCDEF_histo)

muon_effs_iso_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_BCDEF.root" )
muon_effs_iso_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_GH.root" )
muon_effs_iso_BCDEF_histo = muon_effs_iso_BCDEF_file.Get("TightISO_TightID_pt_eta").Get("pt_abseta_ratio")
muon_effs_iso_GH_histo    = muon_effs_iso_GH_file   .Get("TightISO_TightID_pt_eta").Get("pt_abseta_ratio")

print type(muon_effs_iso_BCDEF_histo)

# --- yep, everywhere here Tight ID and ISO is used, since that's the leptons I use

# ------------------------------------------------------------------
# Muon trigger SF
# 

muon_effs_trg_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_RunBtoF.root" )
muon_effs_trg_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_Period4.root" )
muon_effs_trg_BCDEF_histo = muon_effs_trg_BCDEF_file.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins").Get("pt_abseta_ratio")
muon_effs_trg_GH_histo    = muon_effs_trg_GH_file   .Get("IsoMu24_OR_IsoTkMu24_PtEtaBins").Get("pt_abseta_ratio")

print type(muon_effs_trg_BCDEF_file)

muon_effs_id_BCDEF_histo_max_x = muon_effs_id_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_id_BCDEF_histo_max_y = muon_effs_id_BCDEF_histo.GetYaxis().GetXmax()
muon_effs_iso_BCDEF_histo_max_x = muon_effs_iso_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_iso_BCDEF_histo_max_y = muon_effs_iso_BCDEF_histo.GetYaxis().GetXmax()

print muon_effs_id_BCDEF_histo_max_x, muon_effs_id_BCDEF_histo_max_y, muon_effs_iso_BCDEF_histo_max_x, muon_effs_iso_BCDEF_histo_max_y

muon_effs_id_GH_histo_max_x = muon_effs_id_GH_histo.GetXaxis().GetXmax()
muon_effs_id_GH_histo_max_y = muon_effs_id_GH_histo.GetYaxis().GetXmax()
muon_effs_iso_GH_histo_max_x = muon_effs_iso_GH_histo.GetXaxis().GetXmax()
muon_effs_iso_GH_histo_max_y = muon_effs_iso_GH_histo.GetYaxis().GetXmax()

print muon_effs_id_GH_histo_max_x, muon_effs_id_GH_histo_max_y, muon_effs_iso_GH_histo_max_x, muon_effs_iso_GH_histo_max_y

h_weight_mu_trk_bcdef_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200)
h_weight_mu_idd_bcdef_pt = TH1D("weight_mu_idd_bcdef_pt", "", 50, 0, 200)
h_weight_mu_iso_bcdef_pt = TH1D("weight_mu_iso_bcdef_pt", "", 50, 0, 200)
h_weight_mu_trg_bcdef_pt = TH1D("weight_mu_trg_bcdef_pt", "", 50, 0, 200)

h_weight_mu_trk_bcdef_eta = TH1D("weight_mu_trk_bcdef_eta", "", 50, 0, 3)
h_weight_mu_idd_bcdef_eta = TH1D("weight_mu_idd_bcdef_eta", "", 50, 0, 3)
h_weight_mu_iso_bcdef_eta = TH1D("weight_mu_iso_bcdef_eta", "", 50, 0, 3)
h_weight_mu_trg_bcdef_eta = TH1D("weight_mu_trg_bcdef_eta", "", 50, 0, 3)

#
#h_weight_mu_idd_bcdef_nbins = TH2D("weight_mu_idd_bcdef_nbins", "", 100, 0, muon_effs_id_BCDEF_histo.GetSize(), 100, 0, muon_effs_id_BCDEF_histo.GetSize())

def lepton_muon_SF(abs_eta, pt): #, SingleMuon_data_bcdef_fraction, SingleMuon_data_gh_fraction):
    #weight *= 0.98; // average mu trig SF
    #weight *= 0.97; // average mu ID
    weight_muon_sfs = 1

    # muon ID SFs:
    #double weight_muon_sfs = 1;
    bcdef_weight = 1
    gh_weight = 1
    #double SingleMuon_data_bcdef_fraction = 19716.274 / (19716.274 + 15931.028);
    #double SingleMuon_data_gh_fraction    = 15931.028 / (19716.274 + 15931.028);

    #double abs_eta = abs(lep0_eta);
    #double pt = lep0_pt;

    # for control (TODO: remove this later)
    #double bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso,
    #	gh_weight_trk, gh_weight_id, gh_weight_iso;

    # hopefully tracking won't overflow in eta range:
    bcdef_weight_trk = muon_effs_tracking_BCDEF_graph.Eval(abs_eta)
    #h_weight_mu_trk_bcdef_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200),
    h_weight_mu_trk_bcdef_eta.Fill(abs_eta)
    # the id-s totally can overflow:
    bin_x = pt      if pt < muon_effs_id_BCDEF_histo_max_x      else muon_effs_id_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_id_BCDEF_histo_max_y else muon_effs_id_BCDEF_histo_max_y - 0.01 # checked. the binnint is about 0.2 there
    bcdef_weight_id = muon_effs_id_BCDEF_histo.GetBinContent(muon_effs_id_BCDEF_histo.FindBin(bin_x, bin_y))
    h_weight_mu_idd_bcdef_pt .Fill(bin_x)
    h_weight_mu_idd_bcdef_eta.Fill(bin_y)

    # these too:
    bin_x = pt      if pt      < muon_effs_iso_BCDEF_histo_max_x else muon_effs_iso_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_iso_BCDEF_histo_max_y else muon_effs_iso_BCDEF_histo_max_y - 0.01
    bcdef_weight_iso = muon_effs_iso_BCDEF_histo.GetBinContent(muon_effs_iso_BCDEF_histo.FindBin(bin_x, bin_y))
    h_weight_mu_iso_bcdef_pt .Fill(bin_x)
    h_weight_mu_iso_bcdef_eta.Fill(bin_y)
    #bin_x = (pt < muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax() - 1);
    #bin_y = (abs_eta < muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() - 1);
    #bcdef_weight_iso = muon_effs_iso_BCDEF_histo->GetBinContent (muon_effs_iso_BCDEF_histo->FindBin(bin_x, bin_y));

    #fill_1d(string("weight_muon_effs_BCDEF_trk"),  200, 0., 1.1,   bcdef_weight_trk, 1);
    #fill_1d(string("weight_muon_effs_BCDEF_id"),   200, 0., 1.1,   bcdef_weight_id,  1);
    #fill_1d(string("weight_muon_effs_BCDEF_iso"),  200, 0., 1.1,   bcdef_weight_iso, 1);
    #bcdef_weight *= bcdef_weight_trk * bcdef_weight_id * bcdef_weight_iso;

    gh_weight_trk = muon_effs_tracking_GH_graph.Eval(abs_eta)
    #h_weight_mu_trk_gh_eta.Fill(abs_eta)
    # the id-s totally can overflow:
    bin_x = pt      if pt < muon_effs_id_GH_histo_max_x      else muon_effs_id_GH_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_id_GH_histo_max_y else muon_effs_id_GH_histo_max_y - 0.01
    gh_weight_id = muon_effs_id_GH_histo.GetBinContent(muon_effs_id_GH_histo.FindBin(bin_x, bin_y))
    #h_weight_mu_idd_gh_pt .Fill(bin_x)
    #h_weight_mu_idd_gh_eta.Fill(bin_y)
    # these too:
    bin_x = pt      if pt < muon_effs_iso_GH_histo_max_x      else muon_effs_iso_GH_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_iso_GH_histo_max_y else muon_effs_iso_GH_histo_max_y - 0.01
    gh_weight_iso = muon_effs_iso_GH_histo.GetBinContent(muon_effs_iso_GH_histo.FindBin(bin_x, bin_y))
    #h_weight_mu_iso_gh_pt .Fill(bin_x)
    #h_weight_mu_iso_gh_eta.Fill(bin_y)

    return bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso, gh_weight_trk, gh_weight_id, gh_weight_iso


muon_effs_trg_BCDEF_histo_max_x = muon_effs_trg_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_trg_BCDEF_histo_max_y = muon_effs_trg_BCDEF_histo.GetYaxis().GetXmax()

print muon_effs_trg_BCDEF_histo_max_x, muon_effs_trg_BCDEF_histo_max_y

# MUON Trigger
#double lepton_muon_trigger_SF ()
def lepton_muon_trigger_SF(abs_eta, pt): #, double SingleMuon_data_bcdef_fraction, double SingleMuon_data_gh_fraction)
    no_mu_trig = 1;
    #double mu_trig_weight = 1;
    # calculate it the inverse-probbility way
    #double abs_eta = abs(NT_lep_eta_0);
    #double pt = NT_lep_pt_0;
    bin_x = pt      if pt      < muon_effs_trg_BCDEF_histo_max_x else  muon_effs_trg_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_trg_BCDEF_histo_max_y else  muon_effs_trg_BCDEF_histo_max_y - 0.01
    #no_mu_trig *= SingleMuon_data_bcdef_fraction * (1 - muon_effs_trg_BCDEF_histo->GetBinContent( muon_effs_trg_BCDEF_histo->FindBin(bin_x, bin_y) )) +
    #		SingleMuon_data_gh_fraction * (1 - muon_effs_trg_GH_histo->GetBinContent( muon_effs_trg_GH_histo->FindBin(bin_x, bin_y) ));
    #mu_trig_weight = 1 - no_trig; // so for 1 muon it will = to the SF, for 2 there will be a mix
    #fill_1d(string("weight_trigger_no_muon"),  200, 0., 1.1,   no_mu_trig, 1);

    h_weight_mu_trg_bcdef_pt .Fill(bin_x)
    h_weight_mu_trg_bcdef_eta.Fill(bin_y)

    #weight_muon_trig = 1 - no_mu_trig;
    return muon_effs_trg_BCDEF_histo.GetBinContent(muon_effs_trg_BCDEF_histo.FindBin(bin_x, bin_y)), muon_effs_trg_GH_histo.GetBinContent(muon_effs_trg_GH_histo.FindBin(bin_x, bin_y))


'''
/* ---------------------------------------------------
 * now, electrons have
 *      track(reconstruction) efficiency, which is recommended per eta of muon now (however there should be something about N vertices too..
 *      and ID sf
 *      also trigger
 *
 * the trig eff for dilepton case is: apply negative of it for both leptons
 */
'''
logging.info("unpacking electron eff SFs")

electron_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/electron-effs"
gSystem.ExpandPathName(electron_effs_dirname)

electron_effs_tracking_all_file  = TFile(electron_effs_dirname + "/2016_Sept23_ElectronReconstructionSF_egammaEffi.txt_EGM2D.root")
electron_effs_tracking_all_histo = electron_effs_tracking_all_file.Get("EGamma_SF2D")
logging.info("Y tracking (reconstruction)")

# for the selected electrons, Tight ID
# not for Veto
electron_effs_id_all_file  = TFile(electron_effs_dirname + "/2016_Sept23_ElectronID_TightCutBased_egammaEffi.txt_EGM2D.root")
electron_effs_id_all_histo = electron_effs_id_all_file.Get("EGamma_SF2D")
logging.info("Y id")

#analysis/electron-effs/2016_03Feb_TriggerSF_Run2016All_v1.root
electron_effs_trg_all_file  = TFile(electron_effs_dirname + "/2016_03Feb_TriggerSF_Run2016All_v1.root")
electron_effs_trg_all_histo = electron_effs_trg_all_file.Get("Ele27_WPTight_Gsf")
logging.info("Y trigger")

# --- these SFs will be applied to the selected leptons independently

electron_effs_tracking_all_histo_max_y = electron_effs_tracking_all_histo.GetYaxis().GetXmax()
electron_effs_id_all_histo_max_y = electron_effs_id_all_histo.GetYaxis().GetXmax()

h_weight_el_trk_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200)
h_weight_el_idd_pt = TH1D("weight_mu_idd_bcdef_pt", "", 50, 0, 200)
h_weight_mu_iso_bcdef_pt = TH1D("weight_mu_iso_bcdef_pt", "", 50, 0, 200)
h_weight_mu_trg_bcdef_pt = TH1D("weight_mu_trg_bcdef_pt", "", 50, 0, 200)

def lepton_electron_SF(eta, pt):
    #double weight_reco, weight_id;

    # here X axis is eta, Y axis is pt
    # X is from -2.5 to 2.5 -- our eta is up to 2.4 (2.5 in wide ntuples), should be ok
    #double bin_x = (pt < electron_effs_tracking_all_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_id_BCDEF_histo->GetXaxis()->GetXmax() - 1);
    bin_x = eta;
    bin_y = pt if pt < electron_effs_tracking_all_histo_max_y else electron_effs_tracking_all_histo_max_y - 1
    sf_reco = electron_effs_tracking_all_histo.GetBinContent (electron_effs_tracking_all_histo.FindBin(bin_x, bin_y))

    #bin_x = eta;
    bin_y = pt if pt < electron_effs_id_all_histo_max_y else electron_effs_id_all_histo_max_y - 1
    sf_id = electron_effs_id_all_histo.GetBinContent (electron_effs_id_all_histo.FindBin(bin_x, bin_y))

    return sf_reco, sf_id

electron_effs_trg_all_histo_max_x = electron_effs_trg_all_histo.GetXaxis().GetXmax()
electron_effs_trg_all_histo_max_y = electron_effs_trg_all_histo.GetYaxis().GetXmax()

def lepton_electron_trigger_SF(eta, pt):
    #double no_ele_trig = 1;
    # (calculate it the inverse-probbility way)
    # no, just return the SF, assume 1-lepton case
    #pat::Electron& el = selElectrons[i];
    #double eta = el.superCluster()->position().eta();
    #double pt = el.pt();
    # here X axis is pt, Y axis is eta (from -2.5 to 2.5)
    bin_x = pt  if pt  < electron_effs_trg_all_histo_max_x else electron_effs_trg_all_histo_max_x - 1

    bin_y = eta
    if eta > electron_effs_trg_all_histo_max_y:
        bin_y = electron_effs_trg_all_histo_max_y - 0.01
    elif eta < - electron_effs_trg_all_histo_max_y:
        bin_y = - electron_effs_trg_all_histo_max_y + 0.01

    return electron_effs_trg_all_histo.GetBinContent(electron_effs_trg_all_histo.FindBin(bin_x, bin_y))
    #el_trig_weight = 1 - no_trig; // so for 1 lepton it will = to the SF, for 2 there will be a mix
    #fill_1d(string("weight_trigger_no_electron"),  200, 0., 1.1,   no_ele_trig, 1);





'''
b-tagging SF

itak

standalone even if you get to it and  run .L
colapses with
error: constructor for 'BTagCalibration' must explicitly
      initialize the member 'data_'

-- same error as when the wrapper was made...

it should work like:
import ROOT

ROOT.gROOT.ProcessLine('.L BTagCalibrationStandalone.cc+') 

calib = ROOT.BTagCalibration("csv", "CSV.csv")
reader = ROOT.BTagCalibrationReader(0, "central")  # 0 is for loose op
reader.load(calib, 0, "comb")  # 0 is for b flavour, "comb" is the measurement type

# in your event loop
reader.eval(0, 1.2, 30.)  # jet flavor, eta, pt

-- maybe loading ttbar lib the calibrator class should be there?

yup:
>>> import ROOT
>>> ROOT.gROOT.Reset()
>>> ROOT.gROOT.ProcessLine(".L pu_weight.C+") # this is also needed for stuf to run
0L
>>> ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")
0
>>> 
>>> ROOT.BTagCalibration
<class 'ROOT.BTagCalibration'>

-- test tomorrow
'''

with_bSF = False

if with_bSF:
    logging.info("loading b-tagging SF stuff")

    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine(".L pu_weight.C+") # this is also needed for stuf to run
    ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")
    #ROOT.gROOT.ProcessLine("set_bSF_calibrators()")


    #calib = ROOT.BTagCalibration("csv", "CSV.csv")
    #reader = ROOT.BTagCalibrationReader(0, "central")  # 0 is for loose op
    #reader.load(calib, 0, "comb")  # 0 is for b flavour, "comb" is the measurement type

    logging.info("loading b-tagging SF callibration")
    #bCalib_filename = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv"
    bCalib_filename = environ['CMSSW_BASE'] + '/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv'
    #gSystem.ExpandPathName(bCalib_filename)

    logging.info("btag SFs from " + bCalib_filename)
    btagCalib = ROOT.BTagCalibration("CSVv2", bCalib_filename)
    #BTagCalibration* btagCalib; // ("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");

    logging.info("loading Reader")
    btagCal_sys = ROOT.vector('string')()
    btagCal_sys.push_back("up")
    btagCal_sys.push_back("down")
    btagCal = ROOT.BTagCalibrationReader(ROOT.BTagEntry.OP_MEDIUM,  # operating point
    # BTagCalibrationReader btagCal(BTagEntry::OP_TIGHT,  // operating point
    			     "central",            # central sys type
    			     btagCal_sys)      # other sys types
    # the type is:
    # const std::vector<std::string> & otherSysTypes={}


    logging.info("...flavours")
    btagCal.load(btagCalib,        # calibration instance
    	ROOT.BTagEntry.FLAV_B,      # btag flavour
    #            "comb");          # they say "comb" is better precision, but mujets are independent from ttbar dilepton channels
    	"mujets")              #
    btagCal.load(btagCalib,        # calibration instance
    	ROOT.BTagEntry.FLAV_C,      # btag flavour
    	"mujets")              # measurement type
    btagCal.load(btagCalib,        # calibration instance
    	ROOT.BTagEntry.FLAV_UDSG,   # btag flavour
    	"incl")                # measurement type

    logging.info("loaded b-tagging callibration")

    logging.info("loading b-tagging efficiencies")
    #bTaggingEfficiencies_filename = std::string(std::getenv("CMSSW_BASE")) + "/src/UserCode/ttbar-leptons-80X/jobing/configs/b-tagging-efficiencies.root";
    bTaggingEfficiencies_filename = '${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/b-tagging/v9.38-for-b-effs/beff_histos.root'
    gSystem.ExpandPathName(bTaggingEfficiencies_filename)
    bTaggingEfficiencies_file = TFile(bTaggingEfficiencies_filename)

    bEff_histo_b = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_b_hadronFlavour_candidates_tagged')
    bEff_histo_c = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_c_hadronFlavour_candidates_tagged')
    bEff_histo_udsg = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_udsg_hadronFlavour_candidates_tagged')

    logging.info("loaded b-tagging efficiencies")

    h_control_btag_eff_b    = TH1D("control_btag_eff_b",    "", 150, 0, 2)
    h_control_btag_eff_c    = TH1D("control_btag_eff_c",    "", 150, 0, 2)
    h_control_btag_eff_udsg = TH1D("control_btag_eff_udsg", "", 150, 0, 2)

    h_control_btag_weight_b    = TH1D("control_btag_weight_b",    "", 150, 0, 2)
    h_control_btag_weight_c    = TH1D("control_btag_weight_c",    "", 150, 0, 2)
    h_control_btag_weight_udsg = TH1D("control_btag_weight_udsg", "", 150, 0, 2)

    h_control_btag_weight_notag_b    = TH1D("control_btag_weight_notag_b",    "", 150, 0, 2)
    h_control_btag_weight_notag_c    = TH1D("control_btag_weight_notag_c",    "", 150, 0, 2)
    h_control_btag_weight_notag_udsg = TH1D("control_btag_weight_notag_udsg", "", 150, 0, 2)


#def calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
def calc_btag_sf_weight(hasCSVtag, flavId, pt, eta, sys="central"):
    # int flavId=jet.partonFlavour();
    #int flavId=jet.hadronFlavour();
    # also: patJet->genParton().pdgId()
    # fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

    sf, eff = 1.0, 1.0
    # If the jet is tagged -- weight *= SF of the jet
    # if not weight *= (1 - eff*SF)/(1 - eff)
    # 

    if abs(flavId) == 5:
        # get SF for the jet
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_B, eta, pt, 0.)
        # get eff for the jet
        #eff = bTagging_b_jet_efficiency(pt, eta)
        eff = bEff_histo_b.GetBinContent(bEff_histo_b.FindBin(pt, eta))
        h_control_btag_eff_b.Fill(eff)
    elif abs(flavId) == 4:
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_C, eta, pt, 0.)
        #eff = bTagging_c_jet_efficiency(pt, eta)
        eff = bEff_histo_c.GetBinContent(bEff_histo_c.FindBin(pt, eta))
        h_control_btag_eff_c.Fill(eff)
    else:
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_UDSG, eta, pt, 0.)
        #eff = bTagging_udsg_jet_efficiency(pt, eta)
        eff = bEff_histo_udsg.GetBinContent(bEff_histo_udsg.FindBin(pt, eta))
        h_control_btag_eff_udsg.Fill(eff)

    jet_weight_factor = 1.
    if hasCSVtag: # a tagged jet
        jet_weight_factor = sf
        if abs(flavId) == 5:
            h_control_btag_weight_b.Fill(jet_weight_factor)
        elif abs(flavId) == 4:
            h_control_btag_weight_c.Fill(jet_weight_factor)
        else:
            h_control_btag_weight_udsg.Fill(jet_weight_factor)
    else: # not tagged
        # truncate efficiency to 0 and 0.99
        #eff = 0. if eff < 0. else (0.999 if eff > 0.999 else eff)
        jet_weight_factor = (1 - sf*eff) / (1 - eff)
        if abs(flavId) == 5:
            h_control_btag_weight_notag_b.Fill(jet_weight_factor)
        elif abs(flavId) == 4:
            h_control_btag_weight_notag_c.Fill(jet_weight_factor)
        else:
            h_control_btag_weight_notag_udsg.Fill(jet_weight_factor)

    return jet_weight_factor



def top_pT_SF(x):
    # the SF function is SF(x)=exp(a+bx)
    # where x is pT of the top quark (at generation?)
    # sqrt(s) 	channel     	a     	b
    # 7 TeV 	all combined 	0.199 	-0.00166
    # 7 TeV 	l+jets      	0.174 	-0.00137
    # 7 TeV 	dilepton    	0.222 	-0.00197
    # 8 TeV 	all combined 	0.156 	-0.00137
    # 8 TeV 	l+jets       	0.159 	-0.00141
    # 8 TeV 	dilepton     	0.148 	-0.00129
    # 13 TeV	all combined	0.0615	-0.0005
    # -- taking all combined 13 TeV
    a = 0.0615
    b = -0.0005
    return TMath.Exp(a + b*x)

def ttbar_pT_SF(t_pt, tbar_pt):
    return TMath.Sqrt(top_pT_SF(t_pt) * top_pT_SF(tbar_pt))

def transverse_mass(v1, v2):
    return TMath.Sqrt(2*v1.pt()*v2.pt()*(1 - TMath.Cos(v1.phi() - v2.phi())))



def full_loop(t, dtag):
    '''full_loop(t, dtag)

    TTree t
    dtag string
    '''

    save_control = False
    isMC = 'MC' in dtag
    aMCatNLO = 'amcatnlo' in dtag
    isTT = 'TT' in dtag

    #set_bSF_effs_for_dtag(dtag)
    if with_bSF: print bEff_histo_b, bEff_histo_c, bEff_histo_udsg
    #print b_alljet, b_tagged, c_alljet, c_tagged, udsg_alljet, udsg_tagged
    #global bTagging_b_jet_efficiency, bTagging_c_jet_efficiency, bTagging_udsg_jet_efficiency
    #bTagging_b_jet_efficiency = bTagging_X_jet_efficiency(b_alljet, b_tagged)
    #bTagging_c_jet_efficiency = bTagging_X_jet_efficiency(c_alljet, c_tagged)
    #bTagging_udsg_jet_efficiency = bTagging_X_jet_efficiency(udsg_alljet, udsg_tagged)

    control_hs = OrderedDict([
    ('weight_pu', TH1D("weight_pu", "", 50, 0, 2)),
    ('weight_pu_up', TH1D("weight_pu_up", "", 50, 0, 2)),
    ('weight_pu_dn', TH1D("weight_pu_dn", "", 50, 0, 2)),
    ('weight_top_pt', TH1D("weight_top_pt", "", 50, 0, 2)),

    ('weight_mu_trk_bcdef', TH1D("weight_mu_trk_bcdef", "", 50, 0, 2)),
    ('weight_mu_id_bcdef' , TH1D("weight_mu_id_bcdef", "", 50, 0, 2)),
    ('weight_mu_iso_bcdef', TH1D("weight_mu_iso_bcdef", "", 50, 0, 2)),
    ('weight_mu_trg_bcdef', TH1D("weight_mu_trg_bcdef", "", 50, 0, 2)),
    ('weight_mu_all_bcdef', TH1D("weight_mu_all_bcdef", "", 50, 0, 2)),

    ('weight_mu_trk_gh', TH1D("weight_mu_trk_gh", "", 50, 0, 2)),
    ('weight_mu_id_gh' , TH1D("weight_mu_id_gh", "", 50, 0, 2)),
    ('weight_mu_iso_gh', TH1D("weight_mu_iso_gh", "", 50, 0, 2)),
    ('weight_mu_trg_gh', TH1D("weight_mu_trg_gh", "", 50, 0, 2)),
    ('weight_mu_all_gh', TH1D("weight_mu_all_gh", "", 50, 0, 2)),

    ('weight_mu_bSF', TH1D("weight_mu_bSF", "", 50, 0, 2)),
    ('weight_mu_bSF_up', TH1D("weight_mu_bSF_up", "", 50, 0, 2)),
    ('weight_mu_bSF_down', TH1D("weight_mu_bSF_down", "", 50, 0, 2)),

    ('weight_el_trk', TH1D("weight_el_trk", "", 50, 0, 2)),
    ('weight_el_idd', TH1D("weight_el_idd", "", 50, 0, 2)),
    ('weight_el_trg', TH1D("weight_el_trg", "", 50, 0, 2)),
    ('weight_el_all', TH1D("weight_el_all", "", 50, 0, 2)),

    ('weight_el_bSF', TH1D("weight_el_bSF", "", 50, 0, 2)),
    ('weight_el_bSF_up', TH1D("weight_el_bSF_up", "", 50, 0, 2)),
    ('weight_el_bSF_down', TH1D("weight_el_bSF_down", "", 50, 0, 2)),
    ])

    channels = {'el': ['foo'], 'mu': ['foo'], 'mujets': ['foo'],
    'mujets_b': ['foo'], 'taumutauh': ['foo'], 'taumutauh_antiMt_pretau_allb': ['foo'], 'taumutauh_antiMt_pretau': ['foo'], 'taumutauh_antiMt': ['foo'], 'taumutauh_antiMt_OS': ['foo']}

    systematic_names = ('NOMINAL',
                'JESUp'     ,
                'JESDown'   ,
                'JERUp'     ,
                'JERDown'   ,
                'TauESUp'   ,
                'TauESDown' ,
                'bSFUp'     ,
                'bSFDown'   ,
                'PUUp'      ,
                'PUDown'    ,
                )

    out_hs = OrderedDict([((chan, proc, sys), {'met': TH1D('met'+chan+proc+sys, '', 40, 0, 400),
                                               'Mt_lep_met': TH1D('Mt_lep_met'+chan+proc+sys, '', 20, 0, 200),
                                               'Mt_lep_met_d': TH1D('Mt_lep_met_d'+chan+proc+sys, '', 20, 0, 200),
                                               'dijet_trijet_mass': TH1D('dijet_trijet_mass'+chan+proc+sys, '', 20, 0, 400) })
                for chan, procs in channels.items() for proc in procs for sys in systematic_names])

    '''
    for d, histos in out_hs.items():
        for name, histo in histos.items():
            print(d, name)
            histo.Print()
    '''

    print "N entries:", t.GetEntries()
    profile = cProfile.Profile()
    profile.enable()
    for i, ev in enumerate(t):
        '''
        HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0
        '''

        '''
        scheme:

        dtag -> process [subprocesses]
        for event:
            for each SYS -> jet_pts, tau_pts, weight
                which (reco) [channels] it passes
                find (gen) subprocess <------! different subprocess def-s for diff channels
                record distr-s for each
        '''

        # the lepton requirements for all 1-lepton channels:
        # TODO: is there relIso cut now?
        pass_mu = abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0] and ev.lep_p4[0].pt() > 26 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02 # lep_relIso[0]
        pass_el = abs(ev.leps_ID) == 11 and ev.HLT_el and ev.lep_matched_HLT[0] and ev.lep_p4[0].pt() > 29 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02 # lep_relIso[0]

        # only 1-lep channels
        if not (pass_mu or pass_el): continue

        # expensive calls and they don't depend on systematics now
        Mt_lep_met = transverse_mass(ev.lep_p4[0], ev.met_corrected)
        # also
        Mt_lep_met_d = (ev.lep_p4[0] + ev.met_corrected).Mt()

        weight = 1. # common weight of event (1. for data)
        if isMC:
            try:
                weight_pu    = pileup_ratio[ev.nvtx]
                weight_pu_up = pileup_ratio_up[ev.nvtx]
                weight_pu_dn = pileup_ratio_down[ev.nvtx]
                control_hs['weight_pu']   .Fill(weight_pu)
                control_hs['weight_pu_up'].Fill(weight_pu_up)
                control_hs['weight_pu_dn'].Fill(weight_pu_dn)
            except:
                #print i, ev.nvtx
                continue

            if aMCatNLO and ev.aMCatNLO_weight < 0:
                weight *= -1
            #weight_top_pt = 1.
            if isTT:
                weight_top_pt = ttbar_pT_SF(ev.gen_t_pt, ev.gen_tb_pt)
                weight *= weight_top_pt
                control_hs['weight_top_pt']   .Fill(weight_top_pt)

            if pass_mu and isMC:
                mu_sfs = lepton_muon_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                mu_trg_sf = lepton_muon_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                # bcdef gh eras
                weight *= 0.6 * mu_trg_sf[0] * mu_sfs[0] * mu_sfs[1] * mu_sfs[2] + 0.4 * mu_trg_sf[1] * mu_sfs[3] * mu_sfs[4] * mu_sfs[5]

            if pass_el and isMC:
                el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                weight *= el_trg_sf * el_sfs[0] * el_sfs[1]

            #JES corrections
            #[((1-f)*j.pt(), (1+f)*j.pt()) for f,j in zip(ev.jet_jes_correction_relShift, ev.jet_p4)]
            jet_pts_jes_up, jet_pts_jes_down = [], []
            #for f,j in zip(ev.jet_jes_correction_relShift, ev.jet_p4):
            #    jet_pts_jes_up  .append(j.pt()*(1-f))
            #    jet_pts_jes_down.append(j.pt()*(1+f))
            #JER
            #[(j.pt()*(down/factor), j.pt()*(up/factor)) for j, factor, up, down in zip(ev.jet_p4, ev.jet_jer_factor, ev.jet_jer_factor_up, ev.jet_jer_factor_down)]
            jet_pts_jer_up, jet_pts_jer_down = [], []
            #for j, factor, up, down in zip(ev.jet_p4, ev.jet_jer_factor, ev.jet_jer_factor_up, ev.jet_jer_factor_down):
            for i in range(ev.jet_p4.size()):
                j, factor, up, down = ev.jet_p4[i], ev.jet_jer_factor[i], ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
                jet_pts_jer_up  .append(j.pt()*(up/factor))
                jet_pts_jer_down.append(j.pt()*(down/factor))
                jes_shift = ev.jet_jes_correction_relShift[i]
                jet_pts_jes_up  .append(j.pt()*(1-jes_shift))
                jet_pts_jes_down.append(j.pt()*(1+jes_shift))

        jet_pts = [] # nominal jet pts
        for j in ev.jet_p4:
            jet_pts.append(j.pt())

        # number of b-tagged jets is not varied, it's the same for given jets
        # TODO: however, the N b jets should be calculated for JER/JES varied jets
        #       also the variation should be propagated to MET
        N_bjets = 0
        weight_bSF, weight_bSF_up, weight_bSF_down = 1., 1., 1.
        #for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
        for i in range(ev.jet_p4.size()):
            p4, flavId, b_discr = ev.jet_p4[i], ev.jet_hadronFlavour[i], ev.jet_b_discr[i]
            # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
            if p4.pt() < 30: continue
            N_bjets += b_discr > 0.8484
            if isMC and with_bSF:
                weight_bSF      *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta())
                weight_bSF_up   *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up")
                weight_bSF_down *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down")

        # tau pt-s
        # ES correction
        # modes have different correction but the same uncertainty = +- 1.2% = 0.012
        # uncertainties are not correlated, but I'll do correlated UP/DOWN -- all modes UP or all modes DOWN
        tau_pts_corrected = []
        tau_pts_corrected_up = []
        tau_pts_corrected_down = []
        #for i, (p4, DM, IDlev) in enumerate(zip(ev.tau_p4, ev.tau_decayMode, ev.tau_IDlev)):
        for i in range(ev.tau_p4.size()):
            p4, DM, IDlev = ev.tau_p4[i], ev.tau_decayMode[i], ev.tau_IDlev[i]
            if IDlev < 2: continue # only Medium taus
            if DM == 0:
              tau_pts_corrected.append(p4.pt() * 0.995)
              if isMC:
                tau_pts_corrected_up.append(p4.pt() * (0.995 + 0.012))
                tau_pts_corrected_down.append(p4.pt() * (0.995 - 0.012))
            elif DM < 10:
              tau_pts_corrected.append(p4.pt() * 1.011)
              if isMC:
                tau_pts_corrected_up.append(p4.pt()   * (1.011 + 0.012))
                tau_pts_corrected_down.append(p4.pt() * (1.011 - 0.012))
            else:
              tau_pts_corrected.append(p4.pt() * 1.006)
              if isMC:
                tau_pts_corrected_up.append(p4.pt()   * (1.006 + 0.012))
                tau_pts_corrected_down.append(p4.pt() * (1.006 - 0.012))
        #has_medium_tau = any(IDlev > 2 and p4.pt() > 30 for IDlev, p4 in zip(ev.tau_IDlev, ev.tau_p4))
        #has_medium_tau = ev.tau_IDlev.size() > 0 and ev.tau_IDlev[0] > 2 and ev.tau_p4[0].pt() > 30
        #has_medium_tau = bool(tau_pts_corrected)
        #TODO: propagate TES to MET?

        # shape systematics are:
        # corrected jet pts and tau pts
        # weight variations
        #nominal systematics
        # jet pts, tau pts, b weight (=1 for data), pu weight (=1 for data)
        if isMC:
            systematics = {'NOMINAL': [jet_pts, tau_pts_corrected, weight_bSF, weight_pu],
                'JESUp'     : [jet_pts_jes_up,   tau_pts_corrected, weight_bSF, weight_pu],
                'JESDown'   : [jet_pts_jes_down, tau_pts_corrected, weight_bSF, weight_pu],
                'JERUp'     : [jet_pts_jer_up,   tau_pts_corrected, weight_bSF, weight_pu],
                'JERDown'   : [jet_pts_jer_down, tau_pts_corrected, weight_bSF, weight_pu],
                'TauESUp'   : [jet_pts, tau_pts_corrected_up  , weight_bSF, weight_pu],
                'TauESDown' : [jet_pts, tau_pts_corrected_down, weight_bSF, weight_pu],
                'bSFUp'     : [jet_pts, tau_pts_corrected, weight_bSF_up  , weight_pu],
                'bSFDown'   : [jet_pts, tau_pts_corrected, weight_bSF_down, weight_pu],
                'PUUp'      : [jet_pts, tau_pts_corrected, weight_bSF, weight_pu_up],
                'PUDown'    : [jet_pts, tau_pts_corrected, weight_bSF, weight_pu_dn]
                }
        else:
            systematics = {'NOMINAL': [jet_pts, tau_pts_corrected, 1., 1.]}

        # for each systematic
        # pass one of the reco selections
        # check the subprocess
        # store distr

        for sys_name, (jet_pts, tau_pts, weight_bSF, weight_PU) in systematics.items():
            # pass reco selections

            # from channel_distrs
            # --- 'mu': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0]  && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4',
            # met_corrected.pt() > 40 && nbjets > 0 && tau_IDlev[0] > 2 && tau_id[0]*lep_id[0] < 0

            # --- 'mujets':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4
            #&& tau_IDlev[0] > 0 && transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',
            # --- 'mujets_b':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4
            # nbjets > 0 
            # transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',

            # --- 'taumutauh':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] 
            # nbjets == 0 
            # tau_IDlev[0] > 2 && tau_p4[0].pt() > 30 
            # lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 
            # sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) < 40 
            # divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 45 
            # divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 85',
            # --- taumutauh_antiMt <-- anti Mt and ant DY dilep mass
            # --- taumutauh_antiMt_pretau
            # --- taumutauh_antiMt_pretau_allb

            # piecemeal
            # njets > 2 && met_corrected.pt() > 40 && nbjets > 0
            met_pt = ev.met_corrected.pt()
            large_met = met_pt > 40
            large_Mt_lep = Mt_lep_met > 40
            has_bjets = N_bjets > 0
            has_pre_tau = len(tau_pts) > 0
            has_medium_tau = has_pre_tau and tau_pts[0] > 30 and ev.tau_IDlev[0] > 2
            os_lep_med_tau = has_medium_tau and ev.tau_id[0]*ev.lep_id[0] < 0
            # dilep_mass is not TES corrected
            # it's used only in control regions and shouldn't be a big deal
            dy_dilep_mass = has_pre_tau and (45 < (ev.lep_p4[0] + ev.tau_p4[0]).mass() < 85)

            pass_single_lep_presel = large_met and has_bjets and os_lep_med_tau
            pass_mujets    = pass_mu and has_pre_tau and large_Mt_lep
            pass_mujets_b  = pass_mu and has_pre_tau and large_Mt_lep and has_bjets
            pass_taumutauh = pass_mu and os_lep_med_tau and dy_dilep_mass and not large_Mt_lep and not has_bjets
            pass_antiMt_taumutauh_pretau_allb = pass_mu and not dy_dilep_mass and large_Mt_lep
            pass_antiMt_taumutauh_pretau      = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets
            pass_antiMt_taumutauh             = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets and has_medium_tau
            pass_antiMt_taumutauh_OS          = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets and os_lep_med_tau

            # figure out subprocess name and save distrs
            chan_subproc_pairs = []
            if pass_single_lep_presel:
                chan_subproc_pairs.append(('mu' if pass_mu else 'el', 'foo'))
            if pass_mujets:
                chan_subproc_pairs.append(('mujets', 'foo'))
            if pass_mujets_b:
                chan_subproc_pairs.append(('mujets_b', 'foo'))
            if pass_taumutauh:
                chan_subproc_pairs.append(('taumutauh', 'foo'))
            if pass_antiMt_taumutauh_pretau_allb:
                chan_subproc_pairs.append(('taumutauh_antiMt_pretau_allb', 'foo'))
            if pass_antiMt_taumutauh_pretau:
                chan_subproc_pairs.append(('taumutauh_antiMt_pretau', 'foo'))
            if pass_antiMt_taumutauh:
                chan_subproc_pairs.append(('taumutauh_antiMt', 'foo'))
            if pass_antiMt_taumutauh_OS:
                chan_subproc_pairs.append(('taumutauh_antiMt_OS', 'foo'))

            # save distrs
            for chan, proc in chan_subproc_pairs:
                out_hs[(chan, proc, sys_name)]['met'].Fill(met_pt, weight)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met'].Fill(Mt_lep_met, weight)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_d'].Fill(Mt_lep_met_d, weight)
                out_hs[(chan, proc, sys_name)]['dijet_trijet_mass'].Fill(25, weight)

        if save_control:
          if pass_mu:
            # bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso, gh_weight_trk, gh_weight_id, gh_weight_iso
            #mu_sfs = lepton_muon_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #mu_trg_sf = lepton_muon_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_mu_trk_bcdef'].Fill(mu_sfs[0])
            control_hs['weight_mu_id_bcdef'] .Fill(mu_sfs[1])
            control_hs['weight_mu_iso_bcdef'].Fill(mu_sfs[2])
            control_hs['weight_mu_trg_bcdef'].Fill(mu_trg_sf[0])
            control_hs['weight_mu_all_bcdef'].Fill(mu_trg_sf[0] * mu_sfs[0] * mu_sfs[1] * mu_sfs[2])

            control_hs['weight_mu_trk_gh'].Fill(mu_sfs[3])
            control_hs['weight_mu_id_gh'] .Fill(mu_sfs[4])
            control_hs['weight_mu_iso_gh'].Fill(mu_sfs[5])
            control_hs['weight_mu_trg_gh'].Fill(mu_trg_sf[1])
            control_hs['weight_mu_all_gh'].Fill(mu_trg_sf[1] * mu_sfs[3] * mu_sfs[4] * mu_sfs[5])

            for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
                # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
                control_hs['weight_mu_bSF'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta()))
                control_hs['weight_mu_bSF_up']  .Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up"))
                control_hs['weight_mu_bSF_down'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down"))
            # I do need a working standalone b-tag calibrator instead of this jogling
            #b_taggin_SF(jet0_p4.pt(), jet0_p4.eta(), jet0_b_discr, jet0_hadronFlavour, 0.8484)

          elif pass_el:
            #el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_el_trk'].Fill(el_sfs[0])
            control_hs['weight_el_idd'].Fill(el_sfs[1])
            control_hs['weight_el_trg'].Fill(el_trg_sf)
            control_hs['weight_el_all'].Fill(el_trg_sf * el_sfs[0] * el_sfs[1])

            for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
                # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
                control_hs['weight_el_bSF'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta()))
                control_hs['weight_el_bSF_up']  .Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up"))
                control_hs['weight_el_bSF_down'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down"))
    profile.disable()

    profile.print_stats()

    return control_hs, out_hs


def main(argv):
    if '-w' in argv:
        input_filename, nick = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_WJets_madgraph.root", 'wjets'
        dtag = "MC2016_Summer16_WJets_madgraph"
        #init_bSF_call = 'set_bSF_effs_for_dtag("' + dtag + '");'
        #logging.info("init b SFs with: " + init_bSF_call)
        #gROOT.ProcessLine(init_bSF_call)

    else:
        input_filename, nick = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_TTJets_powheg.root" , 'tt'
        dtag = "MC2016_Summer16_TTJets_powheg"
        #init_bSF_call = 'set_bSF_effs_for_dtag("' + dtag + '");'
        #logging.info("init b SFs with: " + init_bSF_call)
        #gROOT.ProcessLine(init_bSF_call)

    #dtag = input_filename.split('/')[-1].split('.')[0]
    logging.info("dtag = " + dtag)
    f = TFile(input_filename)
    #f = TFile('outdir/v12.3/merged-sets/MC2016_Summer16_TTJets_powheg.root')
    t = f.Get('ntupler/reduced_ttree')
    logging.info("N entries = %s" % t.GetEntries())
    c_hs, out_hs = full_loop(t, dtag)
    for name, h in c_hs.items():
        print "%20s %9.5f" % (name, h.GetMean())

    if with_bSF:
        print "eff_b             ", h_control_btag_eff_b             .GetMean()
        print "eff_c             ", h_control_btag_eff_c             .GetMean()
        print "eff_udsg          ", h_control_btag_eff_udsg          .GetMean()
        print "weight_b          ", h_control_btag_weight_b          .GetMean()
        print "weight_c          ", h_control_btag_weight_c          .GetMean()
        print "weight_udsg       ", h_control_btag_weight_udsg       .GetMean()
        print "weight_notag_b    ", h_control_btag_weight_notag_b    .GetMean()
        print "weight_notag_c    ", h_control_btag_weight_notag_c    .GetMean()
        print "weight_notag_udsg ", h_control_btag_weight_notag_udsg .GetMean()

    '''
    for d, histos in out_hs.items():
        for name, histo in histos.items():
            print(d, name)
            histo.Print()
    '''

if __name__ == '__main__':
    main(argv)


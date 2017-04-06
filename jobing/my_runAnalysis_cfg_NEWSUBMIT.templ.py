import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

import PhysicsTools.PythonAnalysis.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('') )

# the PUJetID
from RecoJets.JetProducers.PileupJetIDParams_cfi import cutbased as pu_jetid # this is probably outdated

from RecoJets.JetProducers.PileupJetIDParams_cfi import *

#_stdalgos_4x = cms.VPSet(full,   cutbased,PhilV1)
_stdalgos_5x = cms.VPSet(full_5x,cutbased,PhilV1)

#_chsalgos_4x = cms.VPSet(full,   cutbased) 
_chsalgos_5x = cms.VPSet(full_5x_chs,cutbased)
_chsalgos_74x = cms.VPSet(full_74x_chs,cutbased)
_chsalgos_76x = cms.VPSet(full_76x_chs,cutbased)
_chsalgos_80x = cms.VPSet(full_80x_chs,cutbased)
_chsalgos_81x = cms.VPSet(full_81x_chs,cutbased)

_std_PUJetID_algos = _chsalgos_80x

###### Electron VID
from RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff import *

if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_loose,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_loose.isPOGApproved
if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_medium,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_medium.isPOGApproved
if hasattr(cutBasedElectronID_Spring15_25ns_V1_standalone_tight,'isPOGApproved'):
    del cutBasedElectronID_Spring15_25ns_V1_standalone_tight.isPOGApproved

myVidElectronId = cms.PSet(
    loose = cutBasedElectronID_Spring15_25ns_V1_standalone_loose,
    medium = cutBasedElectronID_Spring15_25ns_V1_standalone_medium,
    tight = cutBasedElectronID_Spring15_25ns_V1_standalone_tight
)
#######

#from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jFullNoQG as mySignalMVA
from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jFull as mySignalMVA
#from UserCode.llvv_fwk.mvaConfig_cfi import ewkzp2jBase as mySignalMVA


pileup2016_direct_old = cms.vdouble(
0.0001963997, 0.0099803283, 0.0140393847, 0.0236650192, 0.0393018959, 0.0309978672, 0.0335953171, 0.0693976520, 0.1645644073, 0.3172568412, 0.4958772343,
0.7212457135, 0.9722619718, 1.1588673286, 1.2629741307, 1.2206370438, 1.1412999298, 1.2127949022, 1.1814431982, 1.2755529410, 1.1356762690, 1.0576495063,
1.0698530227, 1.1106773694, 1.0700410962, 1.1490349905, 1.0178461013, 1.0042396010, 0.8374938723, 0.7532931817, 0.6394506673, 0.5333444236, 0.5741482254,
0.5579844927, 0.7467029777, 1.1796475320, 2.4277076607, 3.2817940991, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
0.0000000000)


# corrected to start from 0
pileup2016_direct2      = cms.vdouble({pileup_reweight_direct})
pileup2016_direct2_up   = cms.vdouble({pileup_reweight_direct_up})
pileup2016_direct2_down = cms.vdouble({pileup_reweight_direct_down})


from os import path as path

theLumiMask = path.expandvars("{lumiMask}")

runProcess = cms.PSet(
    dtag  = cms.string("{dtag}"),
    job_num  = cms.string("{job_num}"),
    input = cms.untracked.vstring({input}),
    outfile = cms.string("{outfile}"),
    outdir = cms.string("{outdir}"),
    isMC = cms.bool({isMC}),
    elHLT_MC     = cms.string("{elHLT_MC}"), # should I just pass v8 and in matching always add * at the end?
    elHLT_Data   = cms.string("{elHLT_Data}"),
    muHLT_MC1    = cms.string("{muHLT_MC1}"),
    muHLT_MC2    = cms.string("{muHLT_MC2}"),
    muHLT_Data1  = cms.string("{muHLT_Data1}"),
    muHLT_Data2  = cms.string("{muHLT_Data2}"),
    jetHLT       = cms.string("{jetHLT}"),
    tau_decayMode = cms.string("{tau_decayMode}"),
    tau_ID        = cms.string("{tau_ID}"),
    tau_againstMuon     = cms.string("{tau_againstMuon}"),
    tau_againstElectron = cms.string("{tau_againstElectron}"),
    jetID   = cms.string("{jetID}"),
    jetPUID = cms.string("{jetPUID}"),
    with_PU  = cms.bool({with_PU}),

    jet_kino_cuts_pt    = cms.double({jet_kino_cuts_pt}),
    jet_kino_cuts_eta   = cms.double({jet_kino_cuts_eta}),
    tau_kino_cuts_pt    = cms.double({tau_kino_cuts_pt}),
    tau_kino_cuts_eta   = cms.double({tau_kino_cuts_eta}),
    jettaufr_jet_kino_cuts_pt    = cms.double({jettaufr_jet_kino_cuts_pt}),
    jettaufr_jet_kino_cuts_eta   = cms.double({jettaufr_jet_kino_cuts_eta}),
    jettaufr_tau_kino_cuts_pt    = cms.double({jettaufr_tau_kino_cuts_pt}),
    jettaufr_tau_kino_cuts_eta   = cms.double({jettaufr_tau_kino_cuts_eta}),
    conf_record_electrons  = cms.bool({conf_record_electrons}),
    conf_record_muons      = cms.bool({conf_record_muons}),
    conf_record_taus_ID    = cms.bool({conf_record_taus_ID}),
    conf_record_taus_kino  = cms.bool({conf_record_taus_kino}),
    withTauIDSFs = cms.bool({withTauIDSFs}),
    PUJetID_algos = cms.VPSet(_std_PUJetID_algos),
    #triggerstudy = cms.bool(@trig),
    triggerstudy = cms.bool(False),
    #xsec = cms.double(@xsec),
    xsec = cms.double(1.0),
    #suffix = cms.string("@suffix"), 
    suffix = cms.string(""), 
    #cprime = cms.double(@cprime),	
    #brnew = cms.double(@brnew),	
    cprime = cms.double(-1),      # 111
    brnew = cms.double(-1),	  # 111
    #mctruthmode = cms.int32(@mctruthmode),
    mctruthmode = cms.int32(0),
    #jacknife = cms.vint32(@jacknife,@jacks),
    #saveSummaryTree = cms.bool(@saveSummaryTree),
    #runSystematics = cms.bool(@runSystematics),
    jacknife = cms.vint32(-1,-1), # 111
    saveSummaryTree = cms.bool(False),
    runSystematics = cms.bool(False),
    #weightsFile = cms.vstring("@weightsFile"),
    weightsFile = cms.vstring(""),
    dirName = cms.string("dataAnalyzer"),
    #useMVA = cms.bool(@useMVA),
    useMVA = cms.bool(False),
    tmvaInput = mySignalMVA,
    muscleDir                         = cms.string('{muscleDir}'),
    muon_effs                      = cms.string('{muon_effs}'),
    electron_effs                  = cms.string('{electron_effs}'),

    jecDir                            = cms.string('{jecDir}'),
    resolutionFile                    = cms.string('{resolutionFile}'),
    scaleFactorFile                   = cms.string('{scaleFactorFile}'),
    dataDriven_tauFakeRates1          = cms.string('{dataDriven_tauFakeRates1}'),
    dataDriven_tauFakeRates2          = cms.string('{dataDriven_tauFakeRates2}'),
    dataDriven_tauFakeRates_dileptons = cms.string('{dataDriven_tauFakeRates_dileptons}'),
    bTaggingEfficiencies              = cms.string('{bTaggingEfficiencies}'),
    tau_fake_rate_histo1_fraction     = cms.double({tau_fake_rate_histo1_fraction}),
    #muscleDir =  cms.string(''),
    #jecDir = cms.string('{{project_dir}}/src/UserCode/llvv_fwk/data/jec/25ns/'),
    #pileup_reweight_direct = pileup_reweight_direct04,
    #pileup_reweight_direct = pileup_reweight_direct05,
    pileup_reweight_direct      = pileup2016_direct2,
    pileup_reweight_direct_down = pileup2016_direct2_down,
    pileup_reweight_direct_up   = pileup2016_direct2_up,
    debug = cms.bool(False),
    debug_len = cms.int32(100),
    lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange(),
    pujetidparas = cms.PSet(pu_jetid),
    electronidparas = cms.PSet(myVidElectronId),
    maxevents = cms.int32(-1) # set to -1 when running on grid. 
)


#STUFF BELLOW IS RUN ONLY IF CRAB IS USED
def lfn_to_pfn(f):
    import subprocess
    proc = subprocess.Popen(["edmFileUtil -d %s" %f],
                            stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    pfn = out.strip()
    return pfn


try:
    import PSet
    fnames = [ lfn_to_pfn(f) for f in list(PSet.process.source.fileNames)]    
    inputFilesPfn =  cms.untracked.vstring(fnames)                
    print inputFilesPfn
    runProcess.input = inputFilesPfn
    runProcess.outfile = cms.string("output.root")
    process.source.fileNames = PSet.process.source.fileNames
except:
    pass

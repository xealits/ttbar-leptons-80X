import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

import PhysicsTools.PythonAnalysis.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('') )
from RecoJets.JetProducers.PileupJetIDParams_cfi import cutbased as pu_jetid


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


pileup2016_direct2 = cms.vdouble(
0.362487235656909, 0.915590545812284, 1.21132246509444, 0.947925601013885, 1.08017689821023, 1.11710415838814, 0.767987830322152, 0.481845466639102, 0.737226697541913, 0.88542246406117,
0.968234885998274, 1.07515806511256, 1.13011024238979, 1.17914792355314, 1.2081221075371, 1.21396216274561, 1.20627443861416, 1.18852262623873, 1.14957508795551, 1.10152205723583,
1.07069598583705, 1.05669944544375, 1.05681198316238, 1.05584403601945, 1.05432721315647, 1.06286485823505, 1.07791370920327, 1.08874014134848, 1.10005952131492, 1.11346269075817,
1.09673780288763, 1.0837667430807, 1.04106716957985, 0.978950963457705, 0.897772213183913, 0.801528577203105, 0.688788869236312, 0.573188095041245, 0.459464269533227, 0.357497209983643,
0.261420202686356, 0.183111774431868, 0.123777241466712, 0.0805931851814276, 0.0512536704394901, 0.0313586759601685, 0.0182617446559743, 0.0106559801209177, 0.00604011048146159, 0.00339109943597573,
0.00188612750875933, 0.00108135060186278, 0.000659922139719069, 0.000451314775712517, 0.0004191474663296, 0.00048513875353278, 0.000627355131909455, 0.000881284041932056, 0.00131150610243989, 0.00185916799436497,
0.00315009081785925, 0.00408586938467248, 0.00470129211639372, 0.00507783574951596, 0.00565732590241663, 0.00561798767875606, 0.0052558923056507, 0.0045310067042737, 0.00397360705412749, 0.00332299989704563,
0.00309960042285951, 0.00279294518017171, 0.00225110283573207, 0.00197674921968105, 0.00185705531699362, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


from os import path as path

theLumiMask = path.expandvars("{lumiMask}")

runProcess = cms.PSet(
    dtag  = cms.string("{dtag}"),
    job_num  = cms.string("{job_num}"),
    input = cms.untracked.vstring({input}),
    outfile = cms.string("{outfile}"),
    outdir = cms.string("{outdir}"),
    isMC = cms.bool({isMC}),
    withTauIDSFs = cms.bool({withTauIDSFs}),
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
    muscleDir =  cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/data/jec/'),
    jecDir = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/data/jec/25ns/'),
    resolutionFile =  cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/data/jec/25ns/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
    scaleFactorFile = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/data/jec/25ns/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
    dataDriven_tauFakeRates1 = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/bin/ttbar-leptons-80X/jet_to_tau_fakerates1.root'),
    dataDriven_tauFakeRates2 = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/bin/ttbar-leptons-80X/jet_to_tau_fakerates2.root'),
    dataDriven_tauFakeRates_dileptons = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/bin/ttbar-leptons-80X/jet_to_tau_fakerates_dileptons.root'),
    bTaggingEfficiencies = cms.string('${{CMSSW_BASE}}/src/UserCode/llvv_fwk/bin/ttbar-leptons-80X/b-tagging-efficiencies.root'),
    tau_fake_rate_histo1_fraction = cms.double(0.5),
    #muscleDir =  cms.string('{{project_dir}}/src/UserCode/llvv_fwk/data/jec/'),
    #jecDir = cms.string('{{project_dir}}/src/UserCode/llvv_fwk/data/jec/25ns/'),
    #pileup_reweight_direct = pileup_reweight_direct04,
    #pileup_reweight_direct = pileup_reweight_direct05,
    pileup_reweight_direct      = pileup2016_direct2,
    pileup_reweight_direct_down = pileup2016_direct2,
    pileup_reweight_direct_up   = pileup2016_direct2,
    datapileup = datapileup_official, 
    datapileupSingleLep = datapileup_official,
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

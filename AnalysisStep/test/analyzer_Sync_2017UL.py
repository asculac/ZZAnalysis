#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2018  # current default = 2017 = Moriond2017
#ELECORRTYPE = "None" # "None" to switch off
#ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = True     #
APPLYJER = True     #
RECORRECTMET = True #
#KINREFIT = True    # control KinZFitter (very slow)
PROCESS_CR = True   # Uncomment to run CR paths and trees
ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # hack for samples missing correct triggers - use with caution
KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = True      # Add tree for Z analysis
#SAMPLENAME = "THW" # For running locally, some samples needs this to be specified (TTZ, THW, WWZ,...) See MCHistoryTools for all samples
FAILED_TREE_LEVEL = True # To print candTree_failed, if you don't want to save it comment this line

PD = ""
MCFILTER = ""

#For DATA:
#IsMC = False
#PD = "DoubleMu"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
### ULTRA-LEGACY PAPER - 2017 sync files
#'/store/mc/RunIISummer20UL17MiniAOD/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/794448BF-6D5B-7149-90C7-2F7D0F3E1DA6.root'
#"/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/40000/DD0D06EF-5C77-8147-84F2-779DEB3E0F52.root"
#"/store/data/Run2018A/EGamma/MINIAOD/UL2018_MiniAODv2-v1/230000/1DC29AF8-7091-4245-A0D8-CFDF650310CC.root"
"/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/B5B6C422-F0C1-844E-974E-7D5CD741C10F.root"
)

#process.calibratedPatElectrons.isSynchronization = cms.bool(True) #process.calibratedPatElectrons.isSynchronization = cms.bool(True) # Not needed anymore since new EGamma smearing is event deterministic
#process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = 1000
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = -1


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


# process.source.eventsToProcess = cms.untracked.VEventRange("1:752:751971")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
#       slimmedElectron = cms.InputTag("slimmedElectrons"),
        electrons = cms.InputTag("appendPhotons:electrons"),
#        RSE = cms.InputTag("appendPhotons:looseElectrons"),
#        TLE = cms.InputTag("appendPhotons:electronstle"), #These are actually photons, should add a photonSrcs section for them.
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZRSE     = cms.InputTag("ZCandlooseEle"),
#        ZTLE     = cms.InputTag("ZCandtle"),
        ZZ  = cms.InputTag("ZZCand"),
#        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
#        ZZTLE     = cms.InputTag("ZZCandtle"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

#process.source.eventsToProcess = cms.untracked.VEventRange("1:6991:641744")

# Create lepton sync file
# process.PlotsZZ.dumpForSync = True;
# process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False
#process.ZZTree.failedTreeLevel = 1

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
# process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

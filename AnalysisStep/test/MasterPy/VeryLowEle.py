#
# Define full paths for ZZ candidates with one loose electron (SR and CRs)
#

SIP_LOOSE = "userFloat('SIP')<100000"
GOODLEPTON_LOOSE = SIP_LOOSE



process.selectedLowElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedLowPtElectrons"),
    cut = cms.string("pt>1.&& electronID('ID')>1.5"),
)

process.bareSoftLowElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("selectedLowElectrons"),
    cut = cms.string("") 
)

process.selectedSlimmedElectrons1 = cms.EDFilter("PATElectronSelector",  
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)<2.5")
)
process.bareSoftElectrons1 = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("selectedSlimmedElectrons1"),
   cut = cms.string("") 
   )
    

process.softLooseElectrons = cms.EDProducer("LowptEleFiller",
   src    = cms.InputTag("bareSoftLowElectrons"),
   src_pf    = cms.InputTag("bareSoftElectrons"),
   src_vertex= cms.InputTag("offlineSlimmedPrimaryVertices"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("pt>4 && abs(eta) < 2.5 && userFloat('dxy')<0.5 && userFloat('dz')<1"),# removed cut because variable is in Spring15 ID  && userFloat('missingHit')<=1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODELECTRON),#GOODLEPTON_LOOSE),#(GOODLEPTON_LOOSE),
        isGoodRegular = cms.string(GOODELECTRON),#GOODELECTRON), # the "regular" (tight) selection
        isIsoFSRUncorr  = cms.string("abs(1)"),#"userFloat('combRelIsoPF')<"+str(ELEISOCUT)), #FIXME add iso and fsr selection for lowpt electrons!!
        isLoose = cms.string("userFloat('isLOOSE')"), #FIXME: I'd set this to  (isGood&&!isGoodTight), that would be clearer I think.
        isPFoverlap = cms.string("userFloat('isPFoverlap')"),
        trackIso = cms.string("userFloat('trackIso')"),
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
   mvaValuesMap = cms.InputTag(""), 
   )
### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftLooseElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softLooseElectrons"),
    preselection = cms.string(''),
    # finalCut (any string-based cut for pat::Electron)
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("userFloat('isGood')"),
           deltaR              = cms.double(0.05),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    finalCut = cms.string(''),
)



process.loose_electrons = cms.Sequence(process.selectedLowElectrons + process.bareSoftLowElectrons 
+ process.selectedSlimmedElectrons1 + process.bareSoftElectrons1 ##this part should be remved and just use the selectedSlimmedElectrons and bareSoftElectrons
+ process.softLooseElectrons + process.cleanSoftLooseElectrons )
process.electrons += process.loose_electrons


 
SKIPPERMUTATIONS = ""#"((daughter(0).pt>daughter(1).pt)||(!daughter(1).masterClone.userFloat('isGoodRegular')))"  #----- ?


#-------------MERGING------- meging low pt collections with electrons and muons, with cleaning duplicates
process.softlooseLeptons = cms.EDProducer("CandViewMerger",
 #   src = cms.VInputTag(cms.InputTag("cleanSoftLooseElectrons"), cms.InputTag("cleanSoftElectrons"))
    src = cms.VInputTag(cms.InputTag("softLeptons") , cms.InputTag("appendPhotons:looseElectrons"))
)
#-----------------------


#KEEPLOOSECOMB_CUT = 'mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())' --keep the one from ZZ4lAnalysis where we keep combinations of tight leptons (passing ID, SIP and ISO)
#"REMOVE2LOWPT="!(daughter(0).masterClone.userFloat('isLoose') && daughter(1).masterClone.userFloat('isLoose'))"

process.bareZCandlooseEle = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('softlooseLeptons@+ softlooseLeptons@-'),
    #cut=cms.string(""),
    cut = cms.string(KEEPLOOSECOMB_CUT),
    checkCharge = cms.bool(True)
)

# process.bareZCandlooseEle = cms.EDProducer("CandViewShallowCloneCombiner",
#     decay = cms.string('appendPhotons:electrons appendPhotons:looseElectrons'),
#     cut = cms.string(SKIPPERMUTATIONS),
#     checkCharge = cms.bool(False)
# )



process.ZCandlooseEle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
        #isLoose = cms.string("abs(1)"),
        isTrueTightRSEZ = cms.string(""),
    )
)

# used to be CandViewShallowCloneCombiner
process.bareZZCandlooseEle= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCandlooseEle ZCandlooseEle'),
    cut = cms.string(LLLLPRESEL),
    checkCharge = cms.bool(True),
    flags = cms.PSet(
        #isLoose = cms.string("abs(1)"),
    )
)


process.ZZCandlooseEle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZZCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandAmong = cms.PSet(isBestCand = cms.string(BESTCAND_AMONG)),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    ZRolesByMass = cms.bool(True),
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
        GoodLeptons =  cms.string(FOURGOODLEPTONS),
        Z2Mass  = cms.string(Z2MASS),
        MAllComb = cms.string(MLLALLCOMB),
        SR = cms.string(SR),
        FullSel70 = cms.string(SR), #Obsolete, use "SR"
        FullSel = cms.string(FULLSEL),
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)

###
#
# CR Need to be implemnted properly! 
#
###
### Trees for control regions only
# what to do in CR where charge is needed?

# Z (OSSF,both e/mu) + LL (any F/C, with no ID/iso); this is the starting point for control regions

# Need to add the CR where the (tight) Z1  and the LL-pair is from the e/mu colections with low pt electrons



#Z2LL_looseEle = "abs(daughter(1).daughter(0).pdgId * daughter(1).daughter(1).pdgId) == 242"   #Z2 = e * e w/o track
Z2LL_SS_looseEle = "daughter(1).daughter(0).pdgId()==daughter(1).daughter(1).pdgId()"       #Z2 = same-sign, same-flavour #= Z2LL_looseEle 
Z2LL_OS_looseEle = "abs(1)"

# Need to drop the SIP cut for teh Z2 candidate
Z2SIP_looseEle = "userFloat('d1.d0.isSIP')< 4 && userFloat('d1.d1.isSIP')"  
CR_BESTCANDBASE_AA_looseEle = ("userFloat('d0.Z1Presel') && userFloat('d0.worstEleIso') <" + str(ELEISOCUT) +
                               "&& userFloat('d0.worstMuIso') <" + str(MUISOCUT)  )
                     #removing SIP cut for bkg estimetion
                     # + "&&"  + Z2SIP_looseEle) # base for AA CR: # Z1 with tight leptons passing SIP and ISO, mass cuts; SIP on Z2

if SELSETUP == "allCutsAtOncePlusSmart" :
    CR_BESTZLLss_looseEle = CR_BESTCANDBASE_AA_looseEle + "&&" + Z2LL_SS + "&&" + CR_Z2MASS + "&&" + MLLALLCOMB + "&&" + "mass>70" + "&&" + "daughter(1).mass>12" + "&&" + SMARTMALLCOMB + "&&" + Z2LL_SS_looseEle 


CR_BASESEL_NOPT20_10 = (CR_Z2MASS + "&&" +              # mass cuts on LL
              MLLALLCOMB + "&&" +             # mass cut on all lepton pairs
              PT20_10    + "&&" +             # pT> 20/10 over all 4 l   ?? we keep it also for low pts!
              "daughter(1).mass>12 &&" +      # mZ2 >12
              "mass>70" )  

# ll, any combination of flavour/charge, for control regions only
process.bareLLCandlooseEle = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('softlooseLeptons softlooseLeptons'),
    cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())'), # protect against ghosts && same flavour
    #cms.string(KEEPLOOSECOMB_CUT),#deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts && skip permuations of the same 2 electrons (See above)
    checkCharge = cms.bool(False)
)

process.LLCandlooseEle = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareLLCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),                     
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        Z1Presel = cms.string(Z1PRESEL),
    )
)
process.bareZLLCandlooseEle= cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('ZCandlooseEle LLCandlooseEle'),
    cut = cms.string(LLLLPRESEL),#'deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())'),
    checkCharge = cms.bool(False)
)
process.ZLLCandlooseEle = cms.EDProducer("ZZCandidateFiller",
    src = cms.InputTag("bareZLLCandlooseEle"),
    sampleType = cms.int32(SAMPLE_TYPE),                    
    setup = cms.int32(LEPTON_SETUP),
    superMelaMass = cms.double(SUPERMELA_MASS),
    isMC = cms.bool(IsMC),
    bestCandComparator = cms.string(BESTCANDCOMPARATOR),
    bestCandAmong = cms.PSet(
      isBestCand    = cms.string("0"), #do not set SR best cand flag
      isBestCRZLLss = cms.string(CR_BESTZLLss_looseEle),
      isBestCRZLLos_2P2F = cms.string("0"),
      isBestCRZLLos_3P1F = cms.string("0")

    ),
    ZRolesByMass = cms.bool(False),  # daughter('Z1') = daughter(0)
    doKinFit = cms.bool(KINREFIT),
    flags = cms.PSet(
      SR = cms.string(SR),
      CRZLLss = cms.string(CR_BASESEL_NOPT20_10  ),        #combine with proper isBestCRZLLss for AA ss/os CRss    
      CRZLLos_2P2F = cms.string("0"), 
      CRZLLos_3P1F = cms.string("0"),        
      number_trackless_electrons = cms.string("abs(1)"),
    ),
    recoProbabilities = cms.vstring(),
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
)





# process.ZlCandlooseEle = cms.EDProducer("PATCandViewShallowCloneCombiner",
#     decay = cms.string('ZCand appendPhotons:looseElectrons'),
#     cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
#                      "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
#                      ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) + # mLL>4 for any pair cause no charge for looseEle OS pair (Giovanni's impl)
#                      ("( %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) + # for trackless electrons we check all the combinations (to be validated)
#                      "daughter(0).masterClone.userFloat('isBestZ') &&" +
#                      "daughter(0).masterClone.userFloat('Z1Presel')"
#                      ),
#     checkCharge = cms.bool(False)
#)



# # Prepare lepton collections
# process.Candidates_loose = cms.Path(
# #       process.muons             +
# #       process.electrons         + process.cleanSoftElectrons +
#        process.appendPhotons      + 
#        process.ZZCandSR           + ~process.ZZCandFilter +
# #       process.fsrPhotons        + process.boostedFsrPhotons +
# #       process.appendPhotons     +
# #       process.softLeptons       +
# #       process.cleanJets         +
# # Build 4-lepton candidates
# #       process.bareZCand         + process.ZCand     +
#        process.bareZCandlooseEle  + process.ZCandlooseEle +  
#        process.bareZZCandlooseEle + process.ZZCandlooseEle
#     )
# Prepare lepton collections
process.Candidates_loose = cms.Path(
#       process.muons             +
#       process.electrons         + process.cleanSoftElectrons +
       process.appendPhotons      + 
       process.ZZCandSR           + ~process.ZZCandFilter +
#       process.fsrPhotons        + process.boostedFsrPhotons +
#       process.appendPhotons     +
#       process.softLeptons       +
#       process.cleanJets         +
# Build 4-lepton candidates
#       process.bareZCand         + process.ZCand     +
       process.softlooseLeptons +
       process.bareZCandlooseEle  + process.ZCandlooseEle +  
       process.bareZZCandlooseEle + process.ZZCandlooseEle
    )

process.CRlooseEle = cms.Sequence(
    #    process.trackless_electrons +
       process.bareZCand + 
       process.softlooseLeptons +
       process.bareZCandlooseEle + process.ZCandlooseEle +
       process.bareLLCandlooseEle       + process.LLCandlooseEle    +
       process.bareZLLCandlooseEle       + process.ZLLCandlooseEle   #+
    #    process.bareZLLCandlooseEleZ1RSE + process.ZLLCandZ1RSE #+ 

    #   process.ZlCandlooseEle 
   )


# process.CRZllooseEle = cms.Sequence(
# #       process.bareZCand         + process.ZCand     +  
#        process.ZlCandlooseEle            
#    )





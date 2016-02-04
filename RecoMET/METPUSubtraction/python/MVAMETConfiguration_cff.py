import FWCore.ParameterSet.Config as cms
import sys


import FWCore.ParameterSet.Config as cms
import sys

from RecoMET.METPUSubtraction.LeptonSelectionTools_cff import applyMuonID
from RecoMET.METPUSubtraction.LeptonSelectionTools_cff import applyElectronID
from RecoMET.METPUSubtraction.LeptonSelectionTools_cff import applyTauID
from RecoMET.METPUSubtraction.LeptonSelectionTools_cff import cleanJetsFromLeptons
from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

def runMVAMET(process,
                 srcMuons =  "slimmedMuons", ## inputMuonCollection
                 muonTypeID    = "Medium", ## type of muon ID to be applied                                                                                                   
                 typeIsoMuons  = "dBeta", ## isolation type to be used for muons                                                                                               
                 relativeIsoCutMuons = 0.25,
                 srcElectrons = "slimmedElectrons", 
                 #electronTypeID= "Tight", 
                 #electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight',
                 typeIsoElectrons = "rhoCorr",
                 relativeIsoCutEletrons = 0.12,
                 srcTaus = "slimmedTaus",
                 tauTypeID = "Loose",
                 jetCollectionPF    = "slimmedJets",
                 dRCleaning= 0.3, 
                 jetPtCut = 10, 
                 jetEtaCut =5.,
                 saveMapForTraining = False
                 ):

    # additional contribution from hadronically decaying taus
    from RecoMET.METPUSubtraction.tausSignificance import tausSignificance, tauMET, tauPFMET, tauDecayProducts
    process.tausSignificance = tausSignificance
    process.tauMET = tauMET
    process.tauPFMET = tauPFMET
    process.tauDecayProducts = tauDecayProducts
   
    relativeIsoCutMuonsLoose = relativeIsoCutMuons+0.05;
    relativeIsoCutEletronsLoose = relativeIsoCutEletrons+0.05;    

    ### run Muon ID
    applyMuonID(process, 
                src   = srcMuons,
                label = muonTypeID, 
                iso_map_charged_hadron  = '',
                iso_map_neutral_hadron  = '',
                iso_map_photon          = '',
                typeIso                 = typeIsoMuons,
                relativeIsolationCutVal = relativeIsoCutMuons
                )


    ### run Electron ID
    """
    applyElectronID(process, 
                    label = electronTypeID, 
                    src   = srcElectrons,
                    iso_map_charged_hadron  = '',
                    iso_map_neutral_hadron  = '',
                    iso_map_photon          = '',
                    typeIso = typeIsoElectrons,
                    electron_id_map = electronID_map,
                    relativeIsolationCutVal = relativeIsoCutEletrons
                    )
    """

    ### run tau ID                                        
    applyTauID( process, 
                label = tauTypeID, 
                src = srcTaus,
                muonCollection     = srcMuons+muonTypeID,
                electronCollection = "slimmedElectrons")

    ## jet lepton cleaning

    cleanJetsFromLeptons(process,
                         label = "Cleaned",
                         jetCollection      = jetCollectionPF,
                         muonCollection     = srcMuons+muonTypeID,
                         electronCollection = "slimmedElectrons",
                         tauCollection      = srcTaus+tauTypeID+"Cleaned",
                         jetPtCut   = jetPtCut,
                         jetEtaCut  = jetEtaCut,
                         dRCleaning = dRCleaning)



    #### Input definitions like in classic MVA MET
    #### tracks from PV
    process.pfChargedPV = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('fromPV && charge !=0'),
                                       src = cms.InputTag("packedPFCandidates")
                                       )
    #### tracks not from PV
    process.pfChargedPU = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('!fromPV && charge !=0'),
                                       src = cms.InputTag("packedPFCandidates")
                                       )
    process.pfNeutrals  = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('charge ==0'),
                                       src = cms.InputTag("packedPFCandidates")
                                       )
    #### Neutrals in Jets passing PU Jet ID
    #### and Neutrals in Jets not passing PU Jet ID
    process.neutralInJets = cms.EDProducer("neutralCandidatePUIDJets",
                                           srcJets = cms.InputTag(jetCollectionPF+"Cleaned"),
                                           srcCandidates = cms.InputTag("pfNeutrals"),
                                           neutralParticlesPVJetsLabel = cms.string("neutralPassingPUIDJets"),
                                           neutralParticlesPUJetsLabel = cms.string("neutralFailingPUIDJets"),
                                           neutralParticlesUnclustered = cms.string("neutralParticlesUnclustered"),
                                           jetPUDIWP = cms.string("user"),
                                           jetPUIDMapLabel = cms.string("fullDiscriminant"))
  

    #### Merge collections to produce corresponding METs
    process.pfMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(          
                                                                                          cms.InputTag("pfChargedPV"),
                                                                                          cms.InputTag("pfChargedPU"),
                                                                                          cms.InputTag("neutralInJets", "neutralPassingPUIDJets"),
                                                                                          cms.InputTag("neutralInJets", "neutralFailingPUIDJets"),
                                                                                          cms.InputTag("neutralInJets", "neutralParticlesUnclustered")
    ))
    #### Track MET
    process.pfTrackMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("pfChargedPV")))
    ## No-PU MET
    process.pfNoPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(        cms.InputTag("pfChargedPV"),
                                                                                          cms.InputTag("neutralInJets", "neutralPassingPUIDJets")))
    ## PU corrected MET
    process.pfPUCorrectedMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag( cms.InputTag("pfChargedPV"), 
                                                                                          cms.InputTag("neutralInJets", "neutralPassingPUIDJets"),
                                                                                          cms.InputTag("neutralInJets", "neutralParticlesUnclustered")))
    ## PU MET
    process.pfPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(          cms.InputTag("pfChargedPU"),
                                                                                          cms.InputTag("neutralInJets", "neutralFailingPUIDJets")))
    
    ##Experimental
    process.pfChargedPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("pfChargedPU")))
    process.pfNeutralPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("neutralInJets", "neutralFailingPUIDJets")))
    process.pfNeutralPVMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("neutralInJets", "neutralPassingPUIDJets")))
    process.pfNeutralUnclusteredMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("neutralInJets", "neutralParticlesUnclustered")))

                                                          
    from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
    patMETsForMVA = patMETs.clone()
    patMETsForMVA.computeMETSignificance = cms.bool(True)
    patMETsForMVA.addGenMET = cms.bool(False)
    patMETsForMVA.srcJets = cms.InputTag(jetCollectionPF)
    #patMETsForMVA.srcLeptons = cms.InputTag("selectedPatJetsAK4PF")
    setattr(patMETsForMVA,"srcLeptons", cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedTaus"))

    for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET", "pfChargedPUMET", "pfNeutralPUMET", "pfNeutralPVMET", "pfNeutralUnclusteredMET"]:
        setattr(process, met, pfMet.clone())
        setattr(getattr(process, met), "src", cms.InputTag(met+"Cands"))
        setattr(getattr(process, met), "alias", cms.string(met))
        setattr(process, "pat"+met, patMETsForMVA.clone())
        setattr(getattr(process, "pat"+met), "metSource", cms.InputTag(met))


    ### MVA MET
    setattr(process,"MVAMET", cms.EDProducer("MVAMET",                                                
                                                referenceMET = cms.InputTag("slimmedMETs"),
                                                debug = cms.bool(False),
                                                requireOS = cms.bool(True),
                                                combineNLeptons = cms.int32(2),
                                                MVAMETLabel = cms.string("MVAMET"),
                                                srcMETs      = cms.VInputTag(
                                                                             cms.InputTag("slimmedMETs"),
                                                                             cms.InputTag("patpfMET"),
                                                                             cms.InputTag("patpfTrackMET"),
                                                                             cms.InputTag("patpfNoPUMET"),
                                                                             cms.InputTag("patpfPUCorrectedMET"),
                                                                             cms.InputTag("patpfPUMET"),
                                                                             cms.InputTag("slimmedMETsPuppi"),
                                                                             cms.InputTag("patpfChargedPUMET"),
                                                                             cms.InputTag("patpfNeutralPUMET"),
                                                                             cms.InputTag("patpfNeutralPVMET"),
                                                                             cms.InputTag("patpfNeutralUnclusteredMET")
                                                                            ),
                                                inputMETFlags = cms.vint32(0,0,0,1,0,0,2,0,2,1,1,1),
                                                srcJets        = cms.InputTag(jetCollectionPF+"Cleaned"),
                                                srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                srcTaus        = cms.InputTag(srcTaus+tauTypeID+"Cleaned"),
                                                srcMuons       = cms.InputTag(srcMuons+muonTypeID),
                                                srcElectrons   = cms.InputTag("slimmedElectrons"),
                                                weightFile     = cms.FileInPath('RecoMET/METPUSubtraction/data/weightfile.root'),
                                                srcLeptons  = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus"), # to produce all possible combinations
                                                tausSignificance = cms.InputTag('tausSignificance', 'METCovariance'),
                                                produceRecoils = cms.bool(False),
                                                saveMap = cms.bool(True)
                                                ))

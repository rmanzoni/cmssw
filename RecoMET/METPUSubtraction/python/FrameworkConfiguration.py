import FWCore.ParameterSet.Config as cms
import sys


import FWCore.ParameterSet.Config as cms
import sys

from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyMuonID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyElectronID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyTauID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanJetsFromLeptons
from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

def runMVAMET(process,
                 processName,
                 isMC,
                 srcMuons =  "slimmedMuons", ## inputMuonCollection
                 muonTypeID    = "Tight", ## type of muon ID to be applied                                                                                                   
                 iso_map_muons = [], ## isolation maps in case they have been re-computed (charged, neutral, photon)                                                         
                 typeIsoMuons  = "dBeta", ## isolation type to be used for muons                                                                                               
                 relativeIsoCutMuons = 0.12,
                 srcElectrons = "slimmedElectrons", 
                 electronTypeID= "Tight", 
                 electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium',
                 electronID_map_loose = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose',
                 iso_map_electrons = [], 
                 typeIsoElectrons = "rhoCorr",
                 relativeIsoCutEletrons = 0.12,
                 srcTaus = "slimmedTaus",
                 tauTypeID = "Loose",
                 doTauCleaning = True,
                 jetCollectionPF    = "slimmedJets",
                 dRCleaning= 0.3, 
                 jetPtCut = 10, 
                 jetEtaCut =5.,
                 genJetCollection = "",
                 cleanGenJets = False,
                 etaCutForMetDiagnostic = 10.,
                 applyTypeICorrection = True, 
                 applyZSelections = True,
                 putRecoilInsideEvent = True
                 ):

    relativeIsoCutMuonsLoose = relativeIsoCutMuons+0.05;
    relativeIsoCutEletronsLoose = relativeIsoCutEletrons+0.05;    

    ### run Muon ID
    if len(iso_map_muons) < 3 :
        
        applyMuonID(process, 
                    src   = srcMuons,
                    label = muonTypeID, 
                    iso_map_charged_hadron  = '',
                    iso_map_neutral_hadron  = '',
                    iso_map_photon          = '',
                    typeIso                 = typeIsoMuons,
                    relativeIsolationCutVal = relativeIsoCutMuons
                    )
    else:

        applyMuonID(process, 
                    src   = srcMuons, 
                    label = muonTypeID, 
                    iso_map_charged_hadron  = iso_map_muons[0],
                    iso_map_neutral_hadron  = iso_map_muons[1],
                    iso_map_photon          = iso_map_muons[2],
                    rho = 'fixedGridRhoFastjetAll',
                    typeIso                 = typeIsoMuons,
                    relativeIsolationCutVal = relativeIsoCutMuons
                )


    ### run Electron ID
    if len(iso_map_electrons) < 3 :
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
    else:
        applyElectronID(process, 
                        label = electronTypeID, 
                        src   = srcElectrons,
                        iso_map_charged_hadron  = iso_map_electrons[0],
                        iso_map_neutral_hadron  = iso_map_electrons[1],
                        iso_map_photon          = iso_map_electrons[2],
                        typeIso = typeIsoElectrons,
                        electron_id_map = electronID_map,
                        relativeIsolationCutVal = relativeIsoCutEletrons
                        )

    ### run tau ID                                        
    if doTauCleaning :
        applyTauID( process, 
                    label = tauTypeID, 
                    src = srcTaus,
                    muonCollection     = srcMuons+muonTypeID,
                    electronCollection = srcElectrons+electronTypeID)

    else:
        applyTauID( process, label = tauTypeID, 
                    src = srcTaus,
                    muonCollection     = "",
                    electronCollection = "")

    ## jet lepton cleaning

    cleanJetsFromLeptons(process,
                         label = "Cleaned",
                         jetCollection      = jetCollectionPF,
                         muonCollection     = srcMuons+muonTypeID,
                         electronCollection = srcElectrons+electronTypeID,
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
                                                              
    from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
    patMETsForMVA = patMETs.clone()
    patMETsForMVA.computeMETSignificance = cms.bool(True)
    patMETsForMVA.addGenMET = cms.bool(False)
    patMETsForMVA.srcJets = cms.InputTag(jetCollectionPF)
    #patMETsForMVA.srcLeptons = cms.InputTag("selectedPatJetsAK4PF")
    setattr(patMETsForMVA,"srcLeptons", cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned"))

    for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET"]:
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
                                                srcMETs      = cms.VInputTag(
                                                                             cms.InputTag("slimmedMETs"),
                                                                             cms.InputTag("patpfMET"),
                                                                             cms.InputTag("patpfTrackMET"),
                                                                             cms.InputTag("patpfNoPUMET"),
                                                                             cms.InputTag("patpfPUCorrectedMET"),
                                                                             cms.InputTag("patpfPUMET"),
                                                                             cms.InputTag("slimmedMETsPuppi"),
                                                                            ),
                                                inputMETFlags = cms.vint32(0,0,0,1,0,0,2,0),
                                                srcJets        = cms.InputTag(jetCollectionPF+"Cleaned"),
                                                srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                srcTaus        = cms.InputTag(srcTaus+tauTypeID+"Cleaned"),
                                                srcMuons       = cms.InputTag(srcMuons+muonTypeID),
                                                srcElectrons   = cms.InputTag("slimmedElectrons"),
                                                weightFile     = cms.FileInPath('JMEAnalysis/JMEValidator/data/weightfile.root'),
                                                #srcLeptons  = cms.VInputTag("LeptonMerge"),
                                                #srcLeptons  = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus"),
                                                srcLeptons  = cms.VInputTag("slimmedMuons"),
                                                #srcLeptons  = cms.VInputTag(cms.InputTag(srcMuons+muonTypeID)),
                                                saveMap = cms.bool(True)
                                                ))

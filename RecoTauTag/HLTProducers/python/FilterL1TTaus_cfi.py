import FWCore.ParameterSet.Config as cms

hltFilterL1TTaus = cms.EDProducer( 
    "FilterL1TTaus",
    L1TauCollection = cms.InputTag("hltCaloStage2Digis", "Tau"),
    minPt           = cms.double(28.),
    maxEta          = cms.double(2.1315),
    iso             = cms.int32(1),
    maxNumber       = cms.int32(99),
    verbose         = cms.bool(False),
)


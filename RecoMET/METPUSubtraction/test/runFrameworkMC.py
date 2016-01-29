import sys

def createProcess(isMC, ## isMC flag
                  processName,
                  globalTag, ## global tag
                  muonTypeID, muonIsoCone, ## muons
                  electronTypeID, ## electrons
                  tauTypeID, ## taus
                  applyZSelections, 
                  jetPtCut,
                  useJECFromLocalDB
                  ):

    process = cms.Process(processName)

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')

    process.GlobalTag.globaltag = globalTag

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ###########################
    ## Electrons and photons ##
    ###########################
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, switchOnVIDPhotonIdProducer, DataFormat, setupAllVIDIdsInModule, setupVIDElectronSelection, setupVIDPhotonSelection

    if not hasattr(process,"egmGsfElectronIDs"):
        electronIdModules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                             'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

        for idMod in electronIdModules:
            setupAllVIDIdsInModule(process, idMod, setupVIDElectronSelection)

    if not hasattr(process,"VersionedPhotonIdProducer"):
        switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

        photonIdModules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

        for idMod in photonIdModules:
            setupAllVIDIdsInModule(process, idMod, setupVIDPhotonSelection)
     ## create the Path
    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)

    return process 


### sample configuration

if not hasattr(sys, 'argv'):
  sys.argv = ["cmsRun", "runFrameworkMC.py"]

import FWCore.ParameterSet.Config as cms
## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC
options.register ('isMC',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to indicate data or MC');
## conditions
options.register ('globalTag',"MCRUN2_74_V9",VarParsing.multiplicity.singleton,VarParsing.varType.string,'input global tag to be used');

## Lepton ID
options.register ('muonTypeID',    "Tight",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'muon ID to be considered for MVA MET analysis ');
options.register ('electronTypeID',"Medium", VarParsing.multiplicity.singleton, VarParsing.varType.string, 'electron ID to be considered for MVA MET analysis ');
options.register ('tauTypeID',     "Loose",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'tau ID to be considered for MVA MET analysis ');
## selections
options.register ('applyZSelections',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply selection for Zll events');
options.register ('jetPtCut',1,VarParsing.multiplicity.singleton, VarParsing.varType.float, 'apply a jet pt cut for mva met input');
## JEC
options.register ('useJECFromDB',         False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'read JEC from the database for special JEC not in GT');
options.parseArguments()

## import the function to create the process
from PhysicsTools.PatAlgos.tools.tauTools import *


## create the process with all the information
process = createProcess(options.isMC, ## MC or data
                        "MVAMET",
                        options.globalTag, ## GT
                        options.muonTypeID, 0.4,## muons
                        options.electronTypeID, ## electrons
                        options.tauTypeID,## taus
                        options.applyZSelections,
                        options.jetPtCut,
                        options.useJECFromDB ## JEC
                        );
from JMEAnalysis.JMEValidator.FrameworkConfiguration import runMVAMET

runMVAMET( process, 
                  "MVAMET",
                  options.isMC,
                  srcMuons = "slimmedMuons", 
                  muonTypeID = "Tight", 
                  iso_map_muons = [], 
                  typeIsoMuons = "dBeta",
                  srcElectrons = "slimmedElectrons", 
                  electronTypeID = "Tight", 
                  electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight',
                  iso_map_electrons = [], 
                  typeIsoElectrons = "rhoCorr",
                  srcTaus = "slimmedTaus", 
                  tauTypeID = "Loose",
                  doTauCleaning = True,
                  jetCollectionPF    = "slimmedJets",
                  dRCleaning = 0.3
                  )
####### files
inputFiles = []
if options.isMC == True:
     inputFiles.append('/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/78750D2F-726D-E511-A7F7-0025905C2CEA.root')
elif options.isMC == False:
     inputFiles.append('/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/E42FEF61-6E27-E511-B93A-02163E0143C0.root')


## set input files
process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(inputFiles);
## output name
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName = cms.string('output.root')
process.TFileService.closeFileFast = cms.untracked.bool(True)


process.MAPAnalyzer =cms.EDAnalyzer('MAPAnalyzer',
                           srcVariables = cms.InputTag("MVAMET"),
                           srcVariableNames = cms.InputTag("MVAMET"),
                           variableNamesToSave = cms.vstring(
                                                            "Jet0_Eta",
                                                            "Jet0_M",
                                                            "Jet0_Phi",
                                                            "Jet0_Pt",
                                                            "Jet1_Eta",
                                                            "Jet1_M",
                                                            "Jet1_Phi",
                                                            "Jet1_Pt",
                                                            "Jet2_Eta",
                                                            "Jet2_Pt",
                                                            #"LongZCorrectedRecoil_Phi",
                                                            #"LongZCorrectedRecoil_Pt",
                                                            #"LongZCorrectedRecoil_sumEt",
                                                            #"LongZCorrectedRecoil_sumEtFraction",
                                                            "NCleanedJets",
                                                            "NVertex",
                                                            "PhiCorrectedRecoil_Phi",
                                                            "PhiCorrectedRecoil_Pt",
                                                            "PhiCorrectedRecoil_sumEt",
                                                            "PhiCorrectedRecoil_sumEtFraction",
                                                            "recoilpatpfMET_Cov00",
                                                            "recoilpatpfMET_Cov11",
                                                            "recoilpatpfMET_Phi",
                                                            "recoilpatpfMET_Pt",
                                                            "recoilpatpfMET_sumEt",
                                                            "recoilpatpfMET_sumEtFraction",
                                                            "recoilpatpfNoPUMET_Cov00",
                                                            "recoilpatpfNoPUMET_Cov11",
                                                            "recoilpatpfNoPUMET_Phi",
                                                            "recoilpatpfNoPUMET_Pt",
                                                            "recoilpatpfNoPUMET_sumEt",
                                                            "recoilpatpfNoPUMET_sumEtFraction",
                                                            "recoilpatpfPUCorrectedMET_Cov00",
                                                            "recoilpatpfPUCorrectedMET_Cov11",
                                                            "recoilpatpfPUCorrectedMET_Phi",
                                                            "recoilpatpfPUCorrectedMET_Pt",
                                                            "recoilpatpfPUCorrectedMET_sumEt",
                                                            "recoilpatpfPUCorrectedMET_sumEtFraction",
                                                            "recoilpatpfPUMET_Cov00",
                                                            "recoilpatpfPUMET_Cov11",
                                                            "recoilpatpfPUMET_Phi",
                                                            "recoilpatpfPUMET_Pt",
                                                            "recoilpatpfPUMET_sumEt",
                                                            "recoilpatpfPUMET_sumEtFraction",
                                                            "recoilpatpfTrackMET_Cov00",
                                                            "recoilpatpfTrackMET_Cov11",
                                                            "recoilpatpfTrackMET_Phi",
                                                            "recoilpatpfTrackMET_Pt",
                                                            "recoilpatpfTrackMET_sumEt",
                                                            "recoilpatpfTrackMET_sumEtFraction",
                                                            "recoilslimmedMETsPuppi_Cov00",
                                                            "recoilslimmedMETsPuppi_Cov11",
                                                            "recoilslimmedMETsPuppi_Phi",
                                                            "recoilslimmedMETsPuppi_Pt",
                                                            "recoilslimmedMETsPuppi_sumEt",
                                                            "recoilslimmedMETsPuppi_sumEtFraction",
                                                            "recoilslimmedMETs_Cov00",
                                                            "recoilslimmedMETs_Cov11",
                                                            "recoilslimmedMETs_Phi",
                                                            "recoilslimmedMETs_Pt",
                                                            "recoilslimmedMETs_sumEt",
                                                            "recoilslimmedMETs_sumEtFraction",
                                                            "Boson_Pt", "Boson_Phi", "Boson_M", "Boson_Eta", "Boson_sumET", "Boson_daughter",
                                                            "nCombinations"
                                                               ) )
process.p = cms.Path()
process.skimmvamet = cms.Sequence( process.MVAMET * process.MAPAnalyzer)
process.p *= (process.skimmvamet)
## logger
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

#! Output and Log                                                                                                                                                            
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
) 
process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('output_particles.root'),
                                      outputCommands = cms.untracked.vstring(
                                                                             'keep *_*_*_*'
                                                                             ),        
                                      SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p'))
                                      )
    
process.out = cms.EndPath(process.output)

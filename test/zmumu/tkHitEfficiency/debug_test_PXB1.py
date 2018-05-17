# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: debug_test_PXB1 --conditions auto:phase1_2018_realistic -n 10 --era Run2_2018 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended -n 3000 --filein file:/store/relval/CMSSW_10_1_0_pre3/RelValZMM_13/GEN-SIM-DIGI-RAW/PU25ns_101X_upgrade2018_realistic_v3-v1/10000/1C85BB66-7429-E811-9F68-0CC47A4C8F06.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/RAW/v1/000/315/506/00000/327B528D-2C4D-E811-B605-FA163E7FCEFA.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('debug_test_PXB1 nevts:3000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

import EventFilter.SiStripRawToDigi.SiStripDigis_cfi
process.siStripDigisRedone = EventFilter.SiStripRawToDigi.SiStripDigis_cfi.siStripDigis.clone(
    UnpackCommonModeValues = True,
)

# Output definition

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)

process.triggerResultsFilter = cms.EDFilter("TriggerResultsFilter",
    daqPartitions = cms.uint32(1),
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tIgnoreMask = cms.bool(False),
    l1tResults = cms.InputTag(""),
    l1techIgnorePrescales = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu20_v*','HLT_IsoMu24_v*','HLT_IsoMu27_v*')
)

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 15 && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)

process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)

process.pseudoProbeSta = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muonsSta"),
    cut = cms.string("outerTrack.isNonnull"),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('tagMuons@+ pseudoProbeSta@-')
)

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi
process.HitCollectorForDebug = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('HitCollectorForDebug'),
    MaxChi2 = cms.double(30.0), ## was 30 ## TO BE TUNED
    nSigma = cms.double(3.), ## was 3  ## TO BE TUNED 
)

process.preselSeq = cms.Sequence(process.tagMuons + process.muonsSta + process.pseudoProbeSta + process.tpPairs + process.onePairSta)

process.clusterInfo = cms.EDAnalyzer("PXB1Hits",
        pairs = cms.InputTag("tpPairs"),
        pixelClusters = cms.InputTag("siPixelClusters"),
        tracker = cms.InputTag("MeasurementTrackerEvent"),
        vertices = cms.InputTag("offlinePrimaryVertices"),
        primaryVertex = cms.InputTag("offlinePrimaryVertices"),
        lumiScalers = cms.InputTag("scalersRawToDigi"),
        # configuraton for refitter
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        RefitDirection = cms.string('oppositeToMomentum'),
        RefitRPCHits = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
        #Propagators
        PropagatorOpposite = cms.string("RungeKuttaTrackerPropagatorOpposite"),
        Chi2MeasurementEstimator = cms.string("HitCollectorForDebug"),
        #Error rescaling
        rescaleError = cms.double(1),
        #SiPixelQuality = cms.string(''),
        badComponentsFile = cms.string('/afs/cern.ch/user/d/ddicroce/work/TrackService/badComponents.txt'),
        ## https://github.com/cms-sw/cmssw/blob/9b7f92a91b55fe1bf3e38435a6afd5b97dea4c9f/RecoLocalTracker/SubCollectionProducers/src/JetCoreClusterSplitter.cc#L139-L153
        debug = cms.untracked.int32(0)
)

process.tagAndProbe = cms.Path(
    process.triggerResultsFilter +
    process.preselSeq +
    process.siStripDigisRedone +
    process.clusterInfo
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("debug_Zmm_lostHits.root"),
    outputCommands = cms.untracked.vstring("keep *",
                                           "drop *_*_*_TagProbe",
                                           "drop *_*_*_RECO",
                                           "keep recoTracks_generalTracks_*_*",
                                           "keep recoMuons_muons_*_*"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("tagAndProbe")),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_debug_Zmm_lostHits_PXB1.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.source.fileNames = [
 'root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/RAW/v1/000/315/506/00000/4E53F186-2C4D-E811-BB02-FA163E0C4B4F.root',
]

# Schedule definition
cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.tagAndProbe)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions 101X_dataRun2_Express_v7 -s RAW2DIGI,L1Reco,RECO,DQM:@trackingOnlyDQM --runUnscheduled --process RECO --data --era Run2_2018 --eventcontent DQM --scenario pp --datatier DQMIO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2018 -n 10 --no_exec --filein root://eoscms.cern.ch//eos/cms/store/data/Commissioning2018/Commissioning/RAW/v1/000/315/104/00000/041F059A-0B48-E811-9F9C-FA163E69FCF6.root --fileout file:step3_inDQM.root --nThreads 8
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.suppressError = cms.untracked.vstring("patTriggerFull")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/data/Commissioning2018/Commissioning/RAW/v1/000/315/104/00000/041F059A-0B48-E811-9F9C-FA163E69FCF6.root'),
    secondaryFileNames = cms.untracked.vstring()
)

JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/DCSOnly/json_DCSONLY_SpecialRun.txt'
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

import EventFilter.SiStripRawToDigi.SiStripDigis_cfi
process.siStripDigisRedone = EventFilter.SiStripRawToDigi.SiStripDigis_cfi.siStripDigis.clone(
    UnpackCommonModeValues = True,
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Express_v7', '')

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
    triggerConditions = cms.vstring('HLT_L1SingleMu*')
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
    fileName = cms.string("output_debug_Zmm_lostHits_lowPU_PXB1_MuEGTOTEM.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.source.fileNames = [
# 'root://cms-xrd-global.cern.ch//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/E6B2831E-7885-E811-A92B-FA163E896A2B.root',
# 'root://cms-xrd-global.cern.ch//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/B8455092-7885-E811-87D7-FA163E11ED71.root',
# 'root://cms-xrd-global.cern.ch//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/8ABBB599-7885-E811-8A86-FA163E9AA870.root',
# 'root://cms-xrd-global.cern.ch//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/ECD5E597-7885-E811-A81B-FA163EE6AFEB.root',
# 'root://cms-xrd-global.cern.ch//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/DA851DAE-7885-E811-95AD-02163E019FF7.root',
 'root://cmsxrootd.fnal.gov//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/E6B2831E-7885-E811-A92B-FA163E896A2B.root',
 'root://cmsxrootd.fnal.gov//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/B8455092-7885-E811-87D7-FA163E11ED71.root',
 'root://cmsxrootd.fnal.gov//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/8ABBB599-7885-E811-8A86-FA163E9AA870.root',
 'root://cmsxrootd.fnal.gov//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/ECD5E597-7885-E811-A81B-FA163EE6AFEB.root',
 'root://cmsxrootd.fnal.gov//store/data/Run2018C/MuonEGammaTOTEM/RAW/v1/000/319/468/00000/DA851DAE-7885-E811-95AD-02163E019FF7.root',
]

# Schedule definition
cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.tagAndProbe)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2018 

#call to customisation function customisePostEra_Run2_2018 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2018(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

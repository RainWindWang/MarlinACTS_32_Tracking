from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper
from k4MarlinWrapper.parseConstants import *
algList = []
evtsvc = EventDataSvc()


CONSTANTS = {
}

parseConstants(CONSTANTS)

read = LcioEvent()
read.OutputLevel = INFO
read.Files = ["input.slcio"]
algList.append(read)

Config = MarlinProcessorWrapper("Config")
Config.OutputLevel = INFO
Config.ProcessorType = "CLICRecoConfig"
Config.Parameters = {
                     "Overlay": ["False"],
                     "OverlayChoices": ["False", "BIB"],
                     "Tracking": ["Truth"],
                     "TrackingChoices": ["Truth", "Conformal"],
                     "VertexUnconstrained": ["OFF"],
                     "VertexUnconstrainedChoices": ["ON", "OFF"]
                     }

EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = INFO
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {
                          "HowOften": ["1"]
                          }

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = INFO
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
                              "FileName": ["lctuple_actsseededckf"],
                              "FileType": ["root"]
                              }

InitDD4hep = MarlinProcessorWrapper("InitDD4hep")
InitDD4hep.OutputLevel = INFO
InitDD4hep.ProcessorType = "InitializeDD4hep"
InitDD4hep.Parameters = {
                         "DD4hepXMLFile": ["../../luxegeo/compact/LUXETrackerAsEndcap.xml"],
                         "EncodingStringParameterName": ["GlobalTrackerReadoutID"]
                         }

MyDDPlanarDigiProcessor = MarlinProcessorWrapper("MyDDPlanarDigiProcessor")
MyDDPlanarDigiProcessor.OutputLevel = INFO
MyDDPlanarDigiProcessor.ProcessorType = "DDPlanarDigiProcessor"
MyDDPlanarDigiProcessor.Parameters = {
                                      "ResolutionU": ["0.005"],
                                      "ResolutionV": ["0.005"],
                                      "SimTrackHitCollectionName": ["SiHits"],
                                      "SimTrkHitRelCollection": ["SiTrackerHitRelations"],
                                      "SubDetectorName": ["LUXETrackerEndcap"],
                                      "TrackerHitCollectionName": ["SiTrackerHits"]
                                      }

MyCKFTracking = MarlinProcessorWrapper("MyCKFTracking")
MyCKFTracking.OutputLevel = INFO
MyCKFTracking.ProcessorType = "ACTSSeededCKFTrackingProc"
MyCKFTracking.Parameters = {
                            "CKF_Chi2CutOff": ["500"],
                            "CKF_NumMeasurementsCutOff": ["1"],
                            "RunCKF": ["True"],
                            "SeedFinding_RMax": ["600"],
                            "SeedFinding_ZMax": ["4300"],
                            "SeedingLayers": ["1", "2", "1", "4", "1", "6"],
                            "TGeoFile": ["../data/LUXETrackerAsEndcap.root"],
                            "TrackCollectionName": ["AllTracks"],
                            "TrackerHitCollectionNames": ["SiTrackerHits"]
                            }

MyTrackDeduper = MarlinProcessorWrapper("MyTrackDeduper")
MyTrackDeduper.OutputLevel = INFO
MyTrackDeduper.ProcessorType = "LUXEDuplicateRemoval"
MyTrackDeduper.Parameters = {
                             "InputTrackCollectionName": ["AllTracks"],
                             "OutputTrackCollectionName": ["Tracks"]
                             }

MyTrackTruth = MarlinProcessorWrapper("MyTrackTruth")
MyTrackTruth.OutputLevel = INFO
MyTrackTruth.ProcessorType = "TrackTruthProc"
MyTrackTruth.Parameters = {
                           "MCParticleCollection": ["MCParticle"],
                           "Particle2TrackRelationName": ["MCParticle_Tracks"],
                           "TrackCollection": ["Tracks"],
                           "TrackerHit2SimTrackerHitRelationName": ["SiTrackerHitRelations"]
                           }

MyACTSTuple = MarlinProcessorWrapper("MyACTSTuple")
MyACTSTuple.OutputLevel = INFO
MyACTSTuple.ProcessorType = "ACTSTuple"
MyACTSTuple.Parameters = {
                          "MCParticleCollection": ["MCParticle"],
                          "SimTrackerHitCollection": ["SiHits"],
                          "TrackerHitCollection": ["SiTrackerHits"]
                          }

MyLCTuple = MarlinProcessorWrapper("MyLCTuple")
MyLCTuple.OutputLevel = INFO
MyLCTuple.ProcessorType = "LCTuple"
MyLCTuple.Parameters = {
                        "LCRelationCollections": ["MCParticle_Tracks"],
                        "LCRelationPrefixes": ["mc2tr"],
                        "MCParticleCollection": ["MCParticle"],
                        "SimTrackerHitCollection": ["SiHits"],
                        "TrackCollection": ["Tracks"],
                        "TrackerHitCollection": ["SiTrackerHits"]
                        }

Output_REC = MarlinProcessorWrapper("Output_REC")
Output_REC.OutputLevel = INFO
Output_REC.ProcessorType = "LCIOOutputProcessor"
Output_REC.Parameters = {
                         "FullSubsetCollections": ["EfficientMCParticles", "InefficientMCParticles"],
                         "KeepCollectionNames": [],
                         "LCIOOutputFile": ["Output_REC.slcio"],
                         "LCIOWriteMode": ["WRITE_NEW"]
                         }

algList.append(MyAIDAProcessor)
algList.append(EventNumber)
algList.append(Config)
algList.append(InitDD4hep)
algList.append(MyDDPlanarDigiProcessor)
algList.append(MyCKFTracking)
algList.append(MyTrackDeduper)
algList.append(MyTrackTruth)
algList.append(MyLCTuple)
algList.append(Output_REC)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax   = 10,
                ExtSvc = [evtsvc],
                OutputLevel=INFO
              )

// -*- C++ -*-
//
// Package:    MuonAnalysis/TagAndProbe
// Class:      TkHits
// 
// Original Author:  Davide Di Croce
//         Created:  Wed, 10 Aug 2016 10:12:34 GMT
//

// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

namespace {
  static std::vector<std::string> sDETS{ "", "PXB", "PXF", "TIB", "TID", "TOB", "TEC" };
  static std::vector<unsigned>    iDETS{ 0,   1,    2,     3,      4,     5,     6   };
  static std::vector<unsigned>    iLAYS{ 0,   3,    2,     4,      3,     6,     9   };
  static std::vector<std::string> sLAYS{ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" };
}

class TkHits : public edm::EDAnalyzer {
  public:
    explicit TkHits(const edm::ParameterSet&);
    ~TkHits();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    const edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;    
    const edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> stripClusterLabel_;
    const edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> stripDigiLabel_;
    const edm::EDGetTokenT<edm::DetSetVector<SiStripRawDigi>> stripCommonModeLabel_;
    edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
    /// Layers to debug
    std::vector<std::string> layersToDebug_;
    /// Track Transformer
    TrackTransformer refitter_;
    std::string propagator_, propagatorOpposite_;
    std::string stripQualityLabel_;
    edm::ESHandle<TrackerTopology> theTrkTopo;
    edm::ESHandle<Propagator> thePropagator, thePropagatorOpposite;
    edm::ESHandle<SiStripQuality> theStripQuality;
    void bar(int i, int scale, int sat=254) {
      std::cout << std::setw(4) << i << " |" ;
      for (int k = 0; k < (i+scale/2)/scale; ++k) std::cout << "=";
      if (i >= sat) std::cout << "X";
      std::cout << std::endl; 
    }
    edm::EDGetTokenT<LumiScalersCollection> m_lumiScalerTag_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryInfoToken_;
    const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
    
    TTree *tree_;
    
    uint32_t run_, lumi_, bx_;
    uint64_t event_;
    double PU_, instLumi_,good_vertices_;
    //double TagMu_px_, TagMu_py_, TagMu_pz_, TagMu_pt_, TagMu_e_, TagMu_eta_, TagMu_phi_;
    //double Mu_px_, Mu_py_, Mu_pz_, Mu_pt_, Mu_e_, Mu_eta_, Mu_phi_;
    //int TagMu_charge_, Mu_charge_;
    double medch1_, bestch1_, medch2_, bestch2_;
    //double cmode1_, cmode2_;

    int missedhit_cent_1_, missedhit_cent_2_, missedhit_edg_1_, missedhit_edg_2_, validhits_pair1_, validhits_pair2_;
    int cluster0modules_pair1_, cluster0modules_pair2_, cluster1modules_pair1_, cluster1modules_pair2_, cluster2modules_pair1_, cluster2modules_pair2_;
    int validhit_edg_1_, validhit_edg_2_, validhit_cent_1_, validhit_cent_2_, ngoodvalid1_, ngoodvalid2_, nbadvalid1_, nbadvalid2_, diffstrip1_,diffstrip2_;
    double nextvhcharge1_, nextvhcharge2_, chargevalid_1_, chargevalid_2_;

    int diffmisstrip1_, diffmisstrip2_;
    double nextmischarge1_, nextmischarge2_;
};


// constructors and destructor
TkHits::TkHits(const edm::ParameterSet& iConfig):
  pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
  stripClusterLabel_(consumes<edmNew::DetSetVector<SiStripCluster>>(iConfig.getParameter<edm::InputTag>("stripClusters"))),
  stripDigiLabel_(consumes<edm::DetSetVector<SiStripDigi>>(iConfig.getParameter<edm::InputTag>("stripDigis"))),
  stripCommonModeLabel_(consumes<edm::DetSetVector<SiStripRawDigi>>(iConfig.getParameter<edm::InputTag>("stripCommonMode"))),
  tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
  layersToDebug_(iConfig.getUntrackedParameter<std::vector<std::string>>("layersToDebug", std::vector<std::string>())),
  refitter_(iConfig),
  propagator_(iConfig.getParameter<std::string>("PropagatorAlong")),
  propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),
  stripQualityLabel_(iConfig.getParameter<std::string>("SiStripQuality")),
  m_lumiScalerTag_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalerTag"))),
  puSummaryInfoToken_ ( consumes<std::vector<PileupSummaryInfo>>( iConfig.getParameter<edm::InputTag>("pileupInfoSummaryInputTag") ) ),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex")))
{ 
  //usesResource("TFileService");
}


TkHits::~TkHits()
{
}

// ------------ method called for each event  ------------
void
TkHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Candidate> > pairs;
  iEvent.getByToken(pairs_, pairs);
  if (pairs->empty()) return;
  iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
  refitter_.setServices(iSetup);
  iSetup.get<TrackingComponentsRecord>().get(propagator_, thePropagator);
  iSetup.get<TrackingComponentsRecord>().get(propagatorOpposite_, thePropagatorOpposite);
  iSetup.get<SiStripQualityRcd>().get(stripQualityLabel_, theStripQuality); 

  edm::Handle<edmNew::DetSetVector<SiStripCluster> > stripC; 
  iEvent.getByToken(stripClusterLabel_, stripC); 
  edm::Handle<edm::DetSetVector<SiStripDigi> > stripD; 
  iEvent.getByToken(stripDigiLabel_, stripD); 
  edm::Handle<edm::DetSetVector<SiStripRawDigi> > stripCM; 
  iEvent.getByToken(stripCommonModeLabel_, stripCM); 

  Handle<MeasurementTrackerEvent> tracker;
  iEvent.getByToken(tracker_, tracker);
  const StMeasurementConditionSet & stripConds = tracker->measurementTracker().stripDetConditions();
  instLumi_ = -999.;
  PU_ = -999.;
  
  validhits_pair1_ = 0;
  validhits_pair2_ = 0;
  validhit_edg_1_ = 0;
  validhit_edg_2_ = 0;
  validhit_cent_1_ = 0;
  validhit_cent_2_ = 0;
  ngoodvalid1_ = 0;
  ngoodvalid2_ = 0;
  nbadvalid1_ = 0;
  nbadvalid2_ = 0;
  nextvhcharge1_ = -999.;
  nextvhcharge2_ = -999.;
  diffstrip1_ = 0;
  diffstrip2_ = 0;
  chargevalid_1_ = -999.;
  chargevalid_2_ = -999.;
  
  missedhit_cent_1_ = 0;
  missedhit_cent_2_ = 0;
  missedhit_edg_1_ = 0;
  missedhit_edg_2_ = 0;
  cluster0modules_pair1_ = 0;
  cluster0modules_pair2_ = 0;
  cluster1modules_pair1_ = 0;
  cluster1modules_pair2_ = 0;
  cluster2modules_pair1_ = 0;
  cluster2modules_pair2_ = 0;
  
  medch1_ = -999.;
  medch2_ = -999.;
  bestch1_ = -999.;
  bestch2_ = -999.;
  //cmode1_ = -999.;
  //cmode2_ = -999.;

  diffmisstrip1_ = 0;
  diffmisstrip2_ = 0;
  nextmischarge1_ = -999.;
  nextmischarge2_ = -999.;


  //TagMu_px_ = TagMu_py_ = TagMu_pz_ = TagMu_pt_ = TagMu_e_ = TagMu_eta_ = TagMu_phi_ = -999.;
  //Mu_px_ = Mu_py_ = Mu_pz_ = Mu_pt_ = Mu_e_ = Mu_eta_ = Mu_phi_ = -999.;
  //TagMu_charge_ = Mu_charge_ = -999;
  
  //std::cout << std::endl;
  
  run_  = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();
  bx_ = iEvent.bunchCrossing();
  good_vertices_ = 0;
  edm::Handle<reco::VertexCollection> recoPrimaryVerticesHandle;
  iEvent.getByToken(pvToken_, recoPrimaryVerticesHandle);
  if (recoPrimaryVerticesHandle.isValid())
    if (recoPrimaryVerticesHandle->size() > 0)
      for (auto v : *recoPrimaryVerticesHandle)
        if (v.ndof() >= 4 && !v.isFake())
          ++good_vertices_;
  //std::cout << "# good_vertices" << good_vertices_ << std::endl; 
  if (iEvent.isRealData()){
    //std::cout << "Is Real Data" << std::endl;
    if (!m_lumiScalerTag_.isUninitialized()){
      edm::Handle<LumiScalersCollection> lumiScaler;
      iEvent.getByToken(m_lumiScalerTag_, lumiScaler);
      if (lumiScaler->begin() != lumiScaler->end())
      instLumi_ = lumiScaler->begin()->instantLumi();
    } else {
      throw cms::Exception("CorruptData") << "[AdditionalEventInfo::produce] AdditionalEventInfo requires a valid LumiScalerCollecion InpuTag" << std::endl;
    }
  } else {
    //std::cout << "Is not Real Data" << std::endl; 
    edm::Handle<std::vector<PileupSummaryInfo>> puSummaryInfoHandle;
    if (iEvent.getByToken(puSummaryInfoToken_, puSummaryInfoHandle)) {
      for (PileupSummaryInfo const & pileup: *puSummaryInfoHandle) {
        // only use the in-time pileup
        if (pileup.getBunchCrossing() == 0) {
          // use the per-event in-time pileup
          PU_ = pileup.getPU_NumInteractions();
        }
      }
    }
    /*int current_event = -1;
    for(size_t iSim=0; iSim != simVertices.size(); ++iSim) {
      const TrackingVertex& simVertex = simVertices[iSim];
      // Associate only to primary vertices of the in-time pileup
      // events (BX=0, first vertex in each of the events)
      if(simVertex.eventId().bunchCrossing() != 0) continue;
      if(simVertex.eventId().event() != current_event) {
        current_event = simVertex.eventId().event();
        PU_ = iSim;
      }
    }*/
  }

  //std::cout << "instLumi = " << instLumi_ << std::endl;
  //std::cout << "PU = " << PU_ << std::endl;
  //std::cout << "NVertixes = " << good_vertices_ << std::endl;

  unsigned int npairs = 0;
  const reco::Candidate *bestpair1 = nullptr, *bestpair2 = nullptr;
  const reco::Muon *mutag1 = nullptr, *mutag2 = nullptr, *mutag3 = nullptr, *mu1 = nullptr, *mu2 = nullptr, *mu3 = nullptr;
  for (const reco::Candidate & pair : *pairs) {
      ++npairs;
      double pair_mass = pair.mass();
      if ( npairs == 1) {
          bestpair1 = &pair;
          mutag1 = (const reco::Muon *) &(*pair.daughter(0)->masterClone());
          mu1 = (const reco::Muon *) &(*pair.daughter(1)->masterClone());
      }
      if ( npairs == 2 || (npairs > 2 && bestpair2 == nullptr) ) {
          mutag3 = (const reco::Muon *) &(*pair.daughter(0)->masterClone());
          mu3 = (const reco::Muon *) &(*pair.daughter(1)->masterClone());   
          if (abs(pair_mass - 91.1876) < abs(bestpair1->mass() - 91.1876)) {
              if ( mutag3->outerTrack().isNull() == 1 ) {
                  bestpair2 = nullptr;
                  mutag2 = nullptr;
                  mu2 = nullptr;
                  mutag1 = mutag3;
                  mu1 = mu3;
                  bestpair1 = &pair;
              } else if ( (mutag1->outerTrack().isNull()==0) && (mutag1->charge() == mu3->charge()) && 
                   (mutag1->outerTrack()->eta()==mu3->eta()) && (mutag3->outerTrack()->eta()==mu1->eta()) ) {
                  bestpair2 = bestpair1;
                  mutag2 = mutag1;
                  mu2 = mu1; 
                  mutag1 = mutag3;
                  mu1 = mu3;
                  bestpair1 = &pair;
              } else {   
                  bestpair2 = nullptr;
                  mutag2 = nullptr;
                  mu2 = nullptr;
                  mutag1 = mutag3;
                  mu1 = mu3;
                  bestpair1 = &pair;
              }
          } else if ( mutag3->outerTrack().isNull() == 1) {
              mutag2 = nullptr;
              mu2 = nullptr;
          } else if ( (mutag1->outerTrack().isNull()==0) && (mutag1->charge() == mu3->charge()) && 
                     (mutag1->outerTrack()->eta()==mu3->eta()) && (mutag3->outerTrack()->eta()==mu1->eta()) ) {
              bestpair2 = &pair;
              mutag2 = mutag3;
              mu2 = mu3;
          }
          mutag3 = nullptr;
          mu3 = nullptr;
      }
      else if ( npairs > 2 && bestpair2 != nullptr ) {
          mutag3 = (const reco::Muon *) &(*pair.daughter(0)->masterClone());
          mu3 = (const reco::Muon *) &(*pair.daughter(1)->masterClone());
          if ( abs(pair_mass - 91.1876) < abs(bestpair1->mass() - 91.1876) ) {
              if ( mutag3->outerTrack().isNull() == 1) {
                  mutag2 = nullptr;
                  mu2 = nullptr;
                  bestpair2 = nullptr;
                  mutag1 = mutag3;
                  mu1 = mu3;
                  bestpair1 = &pair;
              } else if ( (mutag1->outerTrack().isNull()==0) && (mutag1->charge() == mu3->charge()) && 
                          (mutag1->outerTrack()->eta()==mu3->eta()) && (mutag3->outerTrack()->eta()-mu1->eta()) ) {
                  bestpair2 = bestpair1;
                  mutag2 = mutag1;
                  mu2 = mu1;
                  bestpair1 = &pair;
                  mutag1 = mutag3;
                  mu1 = mu3;
              } else if ( (mutag2->outerTrack().isNull()==0) && (mutag2->charge() == mu3->charge()) && 
                          (mutag2->outerTrack()->eta()==mu3->eta()) && (mutag3->outerTrack()->eta()==mu2->eta()) ) {
                  bestpair1 = &pair;
                  mutag1 = mutag3;
                  mu1 = mu3;
              } else {
                  bestpair2 = nullptr;
                  mutag2 = nullptr;
                  mu2 = nullptr;
                  bestpair1 = &pair;
                  mutag1 = mutag3;
                  mu1 = mu3;
              }
          } else if ( mutag3->outerTrack().isNull() == 1) continue;
          else if (abs(pair_mass - 91.1876) < abs(bestpair2->mass() - 91.1876) ) {
              if ( (mutag1->outerTrack().isNull()==0) && (mutag1->charge() == mu3->charge()) && 
                   (mutag1->outerTrack()->eta()==mu3->eta()) && (mutag3->outerTrack()->eta()==mu1->eta()) ) {
                  bestpair2 = &pair;
                  mutag2 = mutag3;
                  mu2 = mu3;
              }
          }
      mutag3 = nullptr;
      mu3 = nullptr;
      }
  }
  
  unsigned int whichpair = 0;
  for (const reco::Candidate & pair : *pairs) {     
      if ( pair.mass()!=bestpair1->mass() && (bestpair2==nullptr || pair.mass()!=bestpair2->mass() ) ) continue;
      const reco::Muon &mu = dynamic_cast<const reco::Muon &>(*pair.daughter(1)->masterClone());
      //std::cout << "Probe Muon with pt " << mu.pt() << " | eta " << mu.eta() << std::endl;
      if (mu.innerTrack().isNull()) continue;
      whichpair++;
      const reco::Track & mutk = *mu.innerTrack();
      //std::cout << "Found hits: " << mutk.found() <<
      //             ", lost hits: " << mutk.hitPattern().numberOfLostStripHits(reco::HitPattern::TRACK_HITS) << std::endl;
      int nhits = mutk.recHitsSize();
      std::vector<Trajectory> traj  = refitter_.transform(mutk);
      if (traj.size() != 1) continue; 
      
      
      // --------  VALID HITS ANALYSIS:  -------------------------------

      const TrackingRecHit *valid = nullptr; DetId herevalid;
      for (int i = 1; i < nhits; ++i) {
        const TrackingRecHit *hitvalid = &* mutk.recHit(i);
        if (hitvalid->isValid()) {
          if (hitvalid->getType() != TrackingRecHit::missing) { 
            valid = hitvalid;
            if (valid == nullptr) continue;
            DetId wherevalid =  valid->geographicalId();
            int subdetvalid = wherevalid.subdetId(), layervalid = theTrkTopo->layer(wherevalid);
            if (!layersToDebug_.empty()) {
              if (std::find(layersToDebug_.begin(), layersToDebug_.end(), sDETS[subdetvalid]+sLAYS[layervalid]) == layersToDebug_.end()) {
                continue;
              }
            }
            //std::cout << "Found hit on" << sDETS[subdetvalid] << layervalid << ", detid " << wherevalid() <<  ", position " << valid->localPosition() << std::endl;
            if (whichpair == 1) validhits_pair1_++;
            if (whichpair == 2) validhits_pair2_++;
            MeasurementDetWithData mdetvalid = tracker->idToDet(wherevalid);
            herevalid = hitvalid->geographicalId();
            TrajectoryStateOnSurface tsosvalid;
            for (const auto &tm : traj.front().measurements()) {
              if (tm.recHit().get() && tm.recHitR().isValid()) {
                if (tm.recHitR().geographicalId() == herevalid) {
                  tsosvalid = tm.updatedState();
                }
              }
            }
            const GeomDet &gdetvalid = mdetvalid.fastGeomDet();
            std::vector<const GeomDet *> gdetsvalid;
            if (typeid(gdetvalid) == typeid(GluedGeomDet)) {
              gdetsvalid.push_back((static_cast<const GluedGeomDet &>(gdetvalid)).monoDet());
              gdetsvalid.push_back((static_cast<const GluedGeomDet &>(gdetvalid)).stereoDet());
            } else {
              gdetsvalid.push_back(&gdetvalid);
            }            
            for (const GeomDet *detvalid : gdetsvalid) {
              mdetvalid = tracker->idToDet(wherevalid);
              int denseIndex = stripConds.find(wherevalid());
              if (denseIndex == stripConds.nDet()) { std::cout << "Module missing in strip conditions set" << std::endl; continue; }
              int nStrips = stripConds.totalStrips(denseIndex);
              //std::cout << "Analyzing module at " << wherevalid() << ", isActive? " << mdetvalid.isActive() << ", strips " << nStrips << std::endl;
              float utrajfirstvalid = 0, utrajlastvalid = 0;
              if (!tsosvalid.isValid()) {
                //std::cout << "  Failed to find Valid Strip ??" << std::endl;
                continue;
              }
              if (!mdetvalid.isActive()) {
                //std::cout << "  Detector is inactive" << std::endl;
                continue;
              }
              //std::cout << "  Bad components on the detector" << std::endl;
              //std::cout << "  APVs (or fibers): ";
              for (unsigned int iapv = 0; iapv < (unsigned (nStrips)/128); ++iapv) {
                if (theStripQuality->IsApvBad(wherevalid(), iapv) || theStripQuality->IsFiberBad(wherevalid(), iapv/2)) std::cout << iapv << " " ;
              } 
              //std::cout << std::endl;
              //std::cout << "  Strips: ";
              SiStripQuality::Range rangevalid = theStripQuality->getRange(wherevalid());
              for (auto stripvalid = rangevalid.first; stripvalid < rangevalid.second; ++stripvalid) std::cout << (*stripvalid) << " " ;
              auto cl_itervalid = stripC->find(wherevalid);
              if (cl_itervalid == stripC->end()) {
                //std::cout << "  ... no strip clusters on this detid" << std::endl;
              } else {
                edmNew::DetSet<SiStripCluster> clustersvalid = *cl_itervalid; 
                auto stripHit = dynamic_cast<const TrackerSingleRecHit *>(valid);
                const SiStripCluster &clusterhit = stripHit->stripCluster();               
                utrajfirstvalid = clusterhit.firstStrip();
                utrajlastvalid = utrajfirstvalid + clusterhit.amplitudes().size() - 1;
                //std::cout << " Valid Strips: first at " << utrajfirstvalid << " , last at " << utrajlastvalid << std::endl;
                if (whichpair == 1){
                  int diffstrip = 6;
                  double nextvhcharge1 = -999.;
                  if ( utrajfirstvalid < 64 || utrajlastvalid > (nStrips-64) ){
                    validhit_edg_1_++;
                    //std::cout << "valid hit in lateral strip, first muon" << std::endl;
                  }
                  if ( utrajfirstvalid >= 64 && utrajlastvalid <= (nStrips-64) ){
                    validhit_cent_1_++;
                    //std::cout << "valid hit in central strip, first muon" << std::endl;
                  }
                  for (const SiStripCluster &clustervalid : clustersvalid) {
                    //std::cout << "  Cluster of " << clustervalid.amplitudes().size() << " strips: " << std::endl;
                    const std::vector<uint8_t> & ampsvalid = clustervalid.amplitudes();
                    for (unsigned int s = clustervalid.firstStrip(), i = 0, e  = ampsvalid.size(); i < e; ++s, ++i) {
                      float dQdx_fromOrigin = siStripClusterTools::chargePerCM(wherevalid, clustervalid, tsosvalid.localParameters());
                      //float invThick = siStripClusterTools::sensorThicknessInverse(wherevalid);
                      //float dQdx_fromOrigin = (ampsvalid[i])*(invThick)*(tsosvalid.localParameters().absdz()); 
                      if ( s == utrajfirstvalid ) {
                        chargevalid_1_ = dQdx_fromOrigin;
                        if ( dQdx_fromOrigin < 2000) ++nbadvalid1_;
                        else if ( dQdx_fromOrigin >= 2000) ++ngoodvalid1_;
                      }
                      //std::cout << "   " << std::setw(4) << s << " | " << (s/128) << " | " << dQdx_fromOrigin << " | "; bar(ampsvalid[i], 2);
                      if ( (s < utrajfirstvalid) && (std::abs(s-utrajfirstvalid) == diffstrip) && (dQdx_fromOrigin > nextvhcharge1) ){
                        nextvhcharge1 = dQdx_fromOrigin;
                      } else if ( (s < utrajfirstvalid) && (std::abs(s-utrajfirstvalid) < diffstrip) ){
                        diffstrip = std::abs(s-utrajfirstvalid);
                        nextvhcharge1 = dQdx_fromOrigin;
                      } else if ( (s > utrajlastvalid) && (std::abs(s-utrajlastvalid) == diffstrip) && (dQdx_fromOrigin > nextvhcharge1) ){
                        nextvhcharge1 = dQdx_fromOrigin;
                      } else if ( (s > utrajlastvalid) && (std::abs(s-utrajlastvalid) < diffstrip) ){
                        diffstrip = std::abs(s-utrajlastvalid);
                        nextvhcharge1 = dQdx_fromOrigin;
                      }
                    }
                  }
                  if ( diffstrip < 6) {
                    diffstrip1_ = diffstrip;
                    nextvhcharge1_ = nextvhcharge1;
                  }
                  //std::cout << " Strip Cluster 1: found hit charge = " << chargevalid_1_ << " | diff. num. strips EXTRA = " << diffstrip1_ << " | EXTRA charge = " << nextvhcharge1_ << std::endl;
                } else if (whichpair == 2){
                  int diffstrip = 6;
                  double nextvhcharge2 = -999.;
                  if ( utrajfirstvalid < 64 || utrajlastvalid > (nStrips-64) ){
                    validhit_edg_2_++;
                    //std::cout << "valid hit in lateral strip, second muon" << std::endl;
                  }
                  if ( utrajfirstvalid >= 64 && utrajlastvalid <= (nStrips-64) ){
                    validhit_cent_2_++;
                    //std::cout << "valid hit in central strip, second muon" << std::endl;
                  }
                  for (const SiStripCluster &clustervalid : clustersvalid) {
                    //std::cout << "  Cluster of " << clustervalid.amplitudes().size() << " strips: " << std::endl;
                    const std::vector<uint8_t> & ampsvalid = clustervalid.amplitudes();
                    for (unsigned int s = clustervalid.firstStrip(), i = 0, e  = ampsvalid.size(); i < e; ++s, ++i) {
                      float dQdx_fromOrigin = siStripClusterTools::chargePerCM(wherevalid, clustervalid, tsosvalid.localParameters());
                      //float invThick = siStripClusterTools::sensorThicknessInverse(wherevalid);
                      //float dQdx_fromOrigin = (ampsvalid[i])*(invThick)*(tsosvalid.localParameters().absdz());
                      if ( s == utrajfirstvalid ) {
                        chargevalid_2_ = dQdx_fromOrigin;
                        if ( dQdx_fromOrigin < 2000) ++nbadvalid2_;
                        else if ( dQdx_fromOrigin >= 2000) ++ngoodvalid2_;
                      }
                      //std::cout << "   " << std::setw(4) << s << " | " << (s/128) << " | " << dQdx_fromOrigin << " | "; bar(ampsvalid[i], 2);
                      if ( (s < utrajfirstvalid) && (std::abs(s-utrajfirstvalid) == diffstrip) && (dQdx_fromOrigin > nextvhcharge2) ){
                        nextvhcharge2 = dQdx_fromOrigin;
                      } else if ( (s < utrajfirstvalid) && (std::abs(s-utrajfirstvalid) < diffstrip) ){
                        diffstrip = std::abs(s-utrajfirstvalid);
                        nextvhcharge2 = dQdx_fromOrigin;
                      } else if ( (s > utrajlastvalid) && (std::abs(s-utrajlastvalid) == diffstrip) && (dQdx_fromOrigin > nextvhcharge2) ){
                        nextvhcharge2 = dQdx_fromOrigin;
                      } else if ( (s > utrajlastvalid) && (std::abs(s-utrajlastvalid) < diffstrip) ){
                        diffstrip = std::abs(s-utrajlastvalid);
                        nextvhcharge2 = dQdx_fromOrigin;
                      }
                    }
                  }
                  if ( diffstrip < 6) {
                    diffstrip2_ = diffstrip; 
                    nextvhcharge2_ = nextvhcharge2;
                  }
                  //std::cout << " Strip Cluster 2: found hit charge = " << chargevalid_2_ << " | diff. num. strips EXTRA = " << diffstrip2_ << " | EXTRA charge = " << nextvhcharge2_ << std::endl;
                }
              }
            }
          }
        } 
      }
          
      
      // --------  INVALID HITS ANALYSIS:  -----------------------------
      
      const TrackingRecHit *invalid = nullptr; DetId previous, next;
      for (int i = 1; i < nhits; ++i) {
        const TrackingRecHit *hit = &* mutk.recHit(i);
        if (hit->getType() == TrackingRecHit::missing) { 
            invalid = hit;
            for (int j = i-1; j >= 0; --j) { 
                hit = &* mutk.recHit(j);
                if (hit->isValid()) {
                    previous = hit->geographicalId();
                    break;
                }
            }
            for (int j = i+1; j < nhits; ++j) {
                hit = &* mutk.recHit(j);
                if (hit->isValid()) {
                    next = hit->geographicalId();
                    break;
                }
            }
            break;
        }
      }
      if (invalid == nullptr || previous.rawId() == 0 || next.rawId() == 0) continue;
      DetId where =  invalid->geographicalId();
      int subdet = where.subdetId(), layer = theTrkTopo->layer(where);
      if (!layersToDebug_.empty()) {
          if (std::find(layersToDebug_.begin(), layersToDebug_.end(), sDETS[subdet]+sLAYS[layer]) == layersToDebug_.end()) {
              continue;
          }
      }
      //std::cout << "Lost hit on " << sDETS[subdet] << layer << ", detid " << where() << ", previous hit on " << previous() << ", next on " << next()<< std::endl;
      MeasurementDetWithData mdet = tracker->idToDet(where);
      const GeomDet &gdet = mdet.fastGeomDet();
      //std::cout << "Lost hit det is a " << typeid(mdet.mdet()).name() << ", geom " << typeid(gdet).name() << std::endl;
      std::vector<const GeomDet *> gdets;
      if (typeid(gdet) == typeid(GluedGeomDet)) {
        gdets.push_back((static_cast<const GluedGeomDet &>(gdet)).monoDet());
        gdets.push_back((static_cast<const GluedGeomDet &>(gdet)).stereoDet());
      } else {
        gdets.push_back(&gdet);
      }
      TrajectoryStateOnSurface tsosBefore, tsosAfter;
      for (const auto &tm : traj.front().measurements()) {
          if (tm.recHit().get() && tm.recHitR().isValid()) {
              if (tm.recHitR().geographicalId() == previous) {
                  tsosBefore = tm.updatedState().isValid() ? tm.updatedState() : tm.forwardPredictedState();
              } else if (tm.recHitR().geographicalId() == next) {
                  tsosAfter = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
              }
          }
      }
      unsigned int foundcluster = 0;
      for (const GeomDet *det : gdets) {
          where = det->geographicalId(); mdet = tracker->idToDet(where);
          int denseIndex = stripConds.find(where());
          if (denseIndex == stripConds.nDet()) { std::cout << "Module missing in strip conditions set" << std::endl; continue; }
          int nStrips = stripConds.totalStrips(denseIndex);
          //std::cout << "Analyzing module at " << where() << ", isActive? " << mdet.isActive() << ", strips " << nStrips << std::endl;
          float utraj = 0; bool pred = false, hascluster = false; 
          //int cmode = -1;
          //float uerr = 0; bool anycluster = false, hasdigi = false, anydigi = false; unsigned int uapv = 0;
          if (!tsosBefore.isValid()) continue;
          TrajectoryStateOnSurface tsos = thePropagator->propagate(tsosBefore, det->surface());
             if (tsos.isValid()) {  
                pred = true;
                utraj = det->topology().measurementPosition( tsos.localPosition() ).x();
                //uerr  = std::sqrt( det->topology().measurementError( tsos.localPosition(), tsos.localError().positionError() ).uu() ); 
                //uapv = std::min<unsigned int>(nStrips-1,std::max<float>(0,utraj))/128;
                //std::cout << "  Searching around strip " << utraj << " +/- " << uerr << "    APV: " << uapv << std::endl;
             } else {
                //std::cout << "  Failed to propagate??" << std::endl;
                continue;
             }
          if (!mdet.isActive()) {
            //std::cout << "  Detector is inactive" << std::endl;
            continue;
          }
          //std::cout << "  Bad components on the detector" << std::endl;
          //std::cout << "  APVs (or fibers): ";
          for (unsigned int iapv = 0; iapv < (unsigned (nStrips)/128); ++iapv) {
            if (theStripQuality->IsApvBad(where(), iapv) || theStripQuality->IsFiberBad(where(), iapv/2)) std::cout << iapv << " " ;
          } 
          //std::cout << std::endl;
          //std::cout << "  Strips: ";
          SiStripQuality::Range range = theStripQuality->getRange(where());
          for (auto strip = range.first; strip < range.second; ++strip) std::cout << (*strip) << " " ;
          //std::cout << std::endl;
          if (whichpair == 1) {
            double beststrip1 = utraj;
            if ( beststrip1 < 64 || beststrip1 > (nStrips-64) ){
              //std::cout << "STRIP 1: " << beststrip1 << std::endl;
              missedhit_edg_1_++;
              //std::cout << "missed hit in lateral strip, first muon" << std::endl;
            }
            if ( (beststrip1 >= 64 && beststrip1 <= (nStrips-64)) ){
              //std::cout << "STRIP 1: " << beststrip1 << std::endl;
              missedhit_cent_1_++;
              //std::cout << "missed hit in central strip, first muon" << std::endl;
            }
          }
          if (whichpair == 2) {
            double beststrip2 = utraj;
            if ( beststrip2 < 64 || beststrip2 > (nStrips-64) ){
              //std::cout << "STRIP 2: " << beststrip2 << std::endl;
              missedhit_edg_2_++;
              //std::cout << "missed hit in lateral strip, second muon" << std::endl;
            }
            if ( beststrip2 >= 64 && beststrip2 <= (nStrips-64) ){
              //std::cout << "STRIP 2: " << beststrip2 << std::endl;
              missedhit_cent_2_++;
              //std::cout << "missed hit in central strip, second muon" << std::endl;
            }
          }
          auto cl_iter = stripC->find(where);
          if (cl_iter == stripC->end()) {
            //std::cout << "  ... no strip clusters on this detid" << std::endl;
          } else {
            edmNew::DetSet<SiStripCluster> clusters = *cl_iter;
            if (whichpair == 1){
              int i1 = 0;
              int diffstrip = 5;
              int diffmisstrip = 6;
              int firstmisstrip = -999;
              int lastmisstrip = -999;
              double nextmischarge1 = -999.;
              double sumch = 0;
              for (const SiStripCluster &cluster : clusters) {
                //std::cout << "  Cluster of " << cluster.amplitudes().size() << " strips: " << std::endl;
                const std::vector<uint8_t> & amps = cluster.amplitudes();
                for (unsigned int s = cluster.firstStrip(), i = 0, e  = amps.size(); i < e; ++s, ++i) {
                  float dQdx_fromOrigin = siStripClusterTools::chargePerCM(where, cluster, tsos.localParameters()); 
                  //std::cout << "   " << std::setw(4) << s << " | " << (s/128) << " | "; bar(amps[i], 2);
                  if (pred && std::abs(s-utraj) < diffstrip) { 
                    hascluster = true; 
                    //anycluster = true;
                    i1++;
                    sumch = sumch + dQdx_fromOrigin;
                    medch1_ = sumch/i1;
                    if ( (std::abs(s-utraj) == diffstrip) && (dQdx_fromOrigin > bestch1_) ){
                      bestch1_ = dQdx_fromOrigin;
                      firstmisstrip = s-i;
                      lastmisstrip = s-i+e-1;
                    } else if ( std::abs(s-utraj) < diffstrip ){
                      diffstrip = std::abs(s-utraj);
                      firstmisstrip = s-i;
                      lastmisstrip = s-i+e-1;
                      bestch1_ = dQdx_fromOrigin;
                    }
                  }
                  //if (pred && (s/128) == uapv) anycluster = true;
                }
                if (hascluster) {
                  for (int s = cluster.firstStrip(), i = 0, e  = amps.size(); i < e; ++s, ++i) {
                    float dQdx_fromOrigin = siStripClusterTools::chargePerCM(where, cluster, tsos.localParameters());
                    if ( (s < firstmisstrip) && (std::abs(s-firstmisstrip) == diffmisstrip) && (dQdx_fromOrigin > nextmischarge1) ){
                      nextmischarge1 = dQdx_fromOrigin;
                    } else if ( (s < firstmisstrip) && (std::abs(s-firstmisstrip) < diffmisstrip) ){
                      diffmisstrip = std::abs(s-firstmisstrip);
                      nextmischarge1 = dQdx_fromOrigin;
                    } else if ( (s > lastmisstrip) && (std::abs(s-lastmisstrip) == diffmisstrip) && (dQdx_fromOrigin > nextmischarge1) ){
                      nextmischarge1 = dQdx_fromOrigin;
                    } else if ( (s > lastmisstrip) && (std::abs(s-lastmisstrip) < diffmisstrip) ){
                      diffmisstrip = std::abs(s-lastmisstrip);
                      nextmischarge1 = dQdx_fromOrigin;
                    }
                  }
                  if ( diffmisstrip < 6) {
                    diffmisstrip1_ = diffmisstrip;
                    nextmischarge1_ = nextmischarge1;
                  }
                }
              }
              //std::cout << " Strip Cluster 1: Med. Charge = " << medch1_ << " | Closer Strip Cluster Charge = " << bestch1_ << std::endl;
            } else if (whichpair == 2){
              int i2 = 0;
              int diffstrip = 5;
              int diffmisstrip = 6;
              int firstmisstrip = -999;
              int lastmisstrip = -999;
              double nextmischarge2 = -999.;
              double sumch = 0;
              for (const SiStripCluster &cluster : clusters) {
                //std::cout << "  Cluster of " << cluster.amplitudes().size() << " strips: " << std::endl;
                const std::vector<uint8_t> & amps = cluster.amplitudes();
                for (unsigned int s = cluster.firstStrip(), i = 0, e  = amps.size(); i < e; ++s, ++i) {
                  float dQdx_fromOrigin = siStripClusterTools::chargePerCM(where, cluster, tsos.localParameters()); 
                  //std::cout << "   " << std::setw(4) << s << " | " << (s/128) << " | "; bar(amps[i], 2);
                  if (pred && std::abs(s-utraj) < diffstrip) {
                    hascluster = true; 
                    //anycluster = true;
                    i2++;
                    sumch = sumch + dQdx_fromOrigin;
                    medch2_ = sumch/i2;
                    if ( (std::abs(s-utraj) == diffstrip) && (dQdx_fromOrigin > bestch2_) ){
                      bestch2_ = dQdx_fromOrigin;
                    } else if ( std::abs(s-utraj) < diffstrip ){
                      diffstrip = std::abs(s-utraj);
                      bestch2_ = dQdx_fromOrigin;
                    }
                  }
                  //if (pred && (s/128) == uapv) anycluster = true;
                }
                if (hascluster) {
                  for (int s = cluster.firstStrip(), i = 0, e  = amps.size(); i < e; ++s, ++i) {
                    float dQdx_fromOrigin = siStripClusterTools::chargePerCM(where, cluster, tsos.localParameters());
                    if ( (s < firstmisstrip) && (std::abs(s-firstmisstrip) == diffmisstrip) && (dQdx_fromOrigin > nextmischarge2) ){
                      nextmischarge2 = dQdx_fromOrigin;
                    } else if ( (s < firstmisstrip) && (std::abs(s-firstmisstrip) < diffmisstrip) ){
                      diffmisstrip = std::abs(s-firstmisstrip);
                      nextmischarge2 = dQdx_fromOrigin;
                    } else if ( (s > lastmisstrip) && (std::abs(s-lastmisstrip) == diffmisstrip) && (dQdx_fromOrigin > nextmischarge2) ){
                      nextmischarge2 = dQdx_fromOrigin;
                    } else if ( (s > lastmisstrip) && (std::abs(s-lastmisstrip) < diffmisstrip) ){
                      diffmisstrip = std::abs(s-lastmisstrip);
                      nextmischarge2 = dQdx_fromOrigin;
                    }
                  }
                  if ( diffmisstrip < 6) {
                    diffmisstrip2_ = diffmisstrip;
                    nextmischarge2_ = nextmischarge2;
                  }
                }
              }
              //std::cout << " Strip Cluster 2: Med. Charge = " << medch2_ << " | Closer Strip Cluster Charge = " << bestch2_ << std::endl;
            }  
          }
          //auto di_iter = stripD->find(where);
          //if (di_iter == stripD->end()) {
              //std::cout << "  ... no strip digis on this detid" << std::endl;
          //} else {
              //std::cout << "  Digis on this detid" << std::endl;
              //const edm::DetSet<SiStripDigi> & digis = *di_iter;
              //for (unsigned int idigi = 0, ndigi = digis.size(); idigi < ndigi; ++idigi) {
                  //if (idigi > 0 && (digis[idigi].strip() > digis[idigi-1].strip()+1)) std::cout << "      ---------------------" << std::endl;
                  //std::cout << "   " << std::setw(4) << digis[idigi].strip() << " | " << (digis[idigi].strip()/128) << " | "; bar(digis[idigi].adc(), 2, 1024);
                  //if (pred && std::abs(digis[idigi].strip()-utraj) < 5) { hasdigi = true; anydigi = true; }
                  //if (pred && (digis[idigi].strip()/128) == uapv) anydigi = true;
              //}
          //}
          //auto cm_iter = stripCM->find(where);
          //if (cm_iter == stripCM->end()) {
          //    //std::cout << "  ... no strip common mode on this detid" << std::endl;
          //} else {
              //std::cout << "  Common mode on this detid" << std::endl;
          //    const edm::DetSet<SiStripRawDigi> & apvs = *cm_iter;
          //    for (unsigned int iapv = 0, napv = apvs.size(); iapv < napv; ++iapv) {
                  //std::cout << "   " << std::setw(4) << "..." << " | " << (iapv) << " | "; bar(apvs[iapv].adc(), 2, 1024);
          //        if (pred && (iapv == uapv)) cmode = apvs[iapv].adc();
          //    }
          //}
          if (pred) {
             //std::cout << "  Summary: " << ( hascluster  ? " cluster" : (anycluster ? " " : " no-clusters")) << 
             //                              ( hasdigi  ? " digi" : (anydigi ? " " : " no-digis"));
             //if (cmode >= 0) {
             //  if (whichpair == 1){
             //    cmode1_ = cmode;
                 //std::cout <<  " CM1 =" << cmode;
             //  } else if (whichpair == 2){
             //    cmode2_ = cmode;
                 //std::cout <<  " CM2 =" << cmode;
             //  }

            // }
          if (utraj < 0 || utraj > nStrips) std::cout << " maybe-outside" ;
             //std::cout << std::endl;
          }
          
          if (hascluster) ++foundcluster; 
      }
      if (whichpair == 1 && foundcluster == 0 ){
          //std::cout << "Missing clusters all modules" << std::endl;
          cluster0modules_pair1_++; 
      }
      if (whichpair == 1 && foundcluster == 1 ){
          //std::cout << "Found cluster in 1 module" << std::endl;
          cluster1modules_pair1_++;
      }
      if (whichpair == 1 && foundcluster == 2 ){
          //std::cout << "Found cluster in 2 modules" << std::endl;
          cluster2modules_pair1_++;
      }
      if (whichpair == 2 && foundcluster == 0 ){
          //std::cout << "Missing clusters all modules" << std::endl;
          cluster0modules_pair2_++;
      }
      if (whichpair == 2 && foundcluster == 1 ){
          //std::cout << "Found cluster in 1 module" << std::endl;
          cluster1modules_pair2_++;
      }
      if (whichpair == 2 && foundcluster == 2 ){
          //std::cout << "Found cluster in 2 modules" << std::endl;
          cluster2modules_pair2_++;
      }

      for (const auto &tm : traj.front().measurements()) {
          if (tm.recHit().get() && !tm.recHit()->isValid()) {
              DetId where =  tm.recHitR().geographicalId();
              int subdet = where.subdetId(), layer = theTrkTopo->layer(where);
              if (!layersToDebug_.empty()) {
                if (std::find(layersToDebug_.begin(), layersToDebug_.end(), sDETS[subdet]+sLAYS[layer]) == layersToDebug_.end()) {
                    continue;
                }
              }
              //std::cout << " missing hit on " << sDETS[subdet] << " layer " << layer << " type = " << tm.recHitR().getType() << ", detid " << where() << std::endl;
          }
      }
  }
  //if (validhits_pair1_ !=0 ||cluster0modules_pair1_ !=0 || cluster1modules_pair1_ != 0 || cluster2modules_pair1_ != 0){
    //std::cout << "LumiBlock: " << iEvent.luminosityBlock() << " | VH: " << validhits_pair1_ << " | 0CMH: " << cluster0modules_pair1_ << " | 1CMH: " << cluster1modules_pair1_ << " | 2CMH: " << cluster2modules_pair1_ << std::endl;
  //}

  //if (validhits_pair2_ !=0 ||cluster0modules_pair2_ !=0 || cluster1modules_pair2_ != 0 || cluster2modules_pair2_ != 0){
    //std::cout << "LumiBlock: " << iEvent.luminosityBlock() << " | VH: " << validhits_pair2_ << " | 0CMH: " << cluster0modules_pair2_ << " | 1CMH: " << cluster1modules_pair2_ << " | 2CMH: " << cluster2modules_pair2_ << std::endl;
  //}
  //std::cout << " Found HITS: VHC1 = " << validhit_cent_1_ << " | VHE1 = " << validhit_edg_1_ << " || VHC2 = " << validhit_cent_2_ << " | VHE2 = " << validhit_edg_2_ << std::endl;
  //std::cout << " Missed HITS: MHC1 = " << missedhit_cent_1_ << " | MHE1 = " << missedhit_edg_1_ << " || MHC2 = " << missedhit_cent_2_ << " | MHE2 = " << missedhit_edg_2_ << std::endl;
  
  if ( validhits_pair1_ != 0 || cluster0modules_pair1_ != 0 || cluster1modules_pair1_ != 0 || cluster2modules_pair1_ != 0 || validhits_pair2_ != 0 ||cluster0modules_pair2_ != 0 || cluster1modules_pair2_ != 0 || cluster2modules_pair2_ != 0 ) {
    tree_->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TkHits::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tkhits_analysis","tkhits_analysis");

  tree_->Branch("run",  &run_,  "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");
  tree_->Branch("bx", &bx_, "bx/i");
  tree_->Branch("instLumi", &instLumi_, "instLumi/d");
  tree_->Branch("PU", &PU_, "PU/d");
  tree_->Branch("good_vertices",  &good_vertices_, "good_vertices/d");

  //tree_->Branch("TagMu_px", &TagMu_px_, "TagMu_px/d");
  //tree_->Branch("TagMu_py", &TagMu_py_, "TagMu_py/d");
  //tree_->Branch("TagMu_pz", &TagMu_pz_, "TagMu_pz/d");
  //tree_->Branch("TagMu_pt", &TagMu_pt_, "TagMu_pt/d");
  //tree_->Branch("TagMu_e", &TagMu_e_, "TagMu_e/d");
  //tree_->Branch("TagMu_eta", &TagMu_eta_, "TagMu_eta/d");
  //tree_->Branch("TagMu_phi", &TagMu_phi_, "TagMu_phi/d");
  //tree_->Branch("TagMu_charge", &TagMu_charge_, "TagMu_charge/i");
  //tree_->Branch("Mu_px", &Mu_px_, "Mu_px/d");
  //tree_->Branch("Mu_py", &Mu_py_, "Mu_py/d");
  //tree_->Branch("Mu_pz", &Mu_pz_, "Mu_pz/d");
  //tree_->Branch("Mu_pt", &Mu_pt_, "Mu_pt/d");
  //tree_->Branch("Mu_e", &Mu_e_, "Mu_e/d");
  //tree_->Branch("Mu_eta", &Mu_eta_, "Mu_eta/d");
  //tree_->Branch("Mu_phi", &Mu_phi_, "Mu_phi/d");
  //tree_->Branch("Mu_charge", &Mu_charge_, "Mu_charge/i");
 
  tree_->Branch("validhits_pair1", &validhits_pair1_, "validhits_pair1/i");
  tree_->Branch("validhits_pair2", &validhits_pair2_, "validhits_pair2/i");
  
  tree_->Branch("ngoodvalid1", &ngoodvalid1_, "ngoodvalid1/i");
  tree_->Branch("ngoodvalid2", &ngoodvalid2_, "ngoodvalid2/i");
  tree_->Branch("nbadvalid1", &nbadvalid1_, "nbadvalid1/i");
  tree_->Branch("nbadvalid2", &nbadvalid2_, "nbadvalid2/i");
  tree_->Branch("chargevalid_1", &chargevalid_1_, "chargevalid_1/d");
  tree_->Branch("chargevalid_2", &chargevalid_2_, "chargevalid_2/d");
  tree_->Branch("diffstrip1", &diffstrip1_, "diffstrip1/i");
  tree_->Branch("diffstrip2", &diffstrip2_, "diffstrip2/i");
  tree_->Branch("nextvhcharge1", &nextvhcharge1_, "nextvhcharge1/d");
  tree_->Branch("nextvhcharge2", &nextvhcharge2_, "nextvhcharge2/d");
  
  tree_->Branch("validhit_cent_1", &validhit_cent_1_, "validhit_cent_1/i");
  tree_->Branch("validhit_cent_2", &validhit_cent_2_, "validhit_cent_2/i");
  tree_->Branch("validhit_edg_1", &validhit_edg_1_, "validhit_edg_1/i");
  tree_->Branch("validhit_edg_2", &validhit_edg_2_, "validhit_edg_2/i");
  tree_->Branch("missedhit_cent_1", &missedhit_cent_1_, "missedhit_cent_1/i");
  tree_->Branch("missedhit_cent_2", &missedhit_cent_2_, "missedhit_cent_2/i");
  tree_->Branch("missedhit_edg_1", &missedhit_edg_1_, "missedhit_edg_1/i");
  tree_->Branch("missedhit_edg_2", &missedhit_edg_2_, "missedhit_edg_2/i");
  
  tree_->Branch("cluster0modules_pair1", &cluster0modules_pair1_, "cluster0modules_pair1/i");
  tree_->Branch("cluster0modules_pair2", &cluster0modules_pair2_, "cluster0modules_pair2/i");
  tree_->Branch("cluster1modules_pair1", &cluster1modules_pair1_, "cluster1modules_pair1/i");
  tree_->Branch("cluster1modules_pair2", &cluster1modules_pair2_, "cluster1modules_pair2/i");
  tree_->Branch("cluster2modules_pair1", &cluster2modules_pair1_, "cluster2modules_pair1/i");
  tree_->Branch("cluster2modules_pair2", &cluster2modules_pair2_, "cluster2modules_pair2/i");
  
  tree_->Branch("medch1", &medch1_, "medch1/d");
  tree_->Branch("medch2", &medch2_, "medch2/d");
  tree_->Branch("bestch1", &bestch1_, "bestch1/d");
  tree_->Branch("bestch2", &bestch2_, "bestch2/d");
  //tree_->Branch("cmode1", &cmode1_, "cmode1/d");
  //tree_->Branch("cmode2", &cmode2_, "cmode2/d");

  tree_->Branch("diffmisstrip1", &diffmisstrip1_, "diffmisstrip1/i");
  tree_->Branch("diffmisstrip2", &diffmisstrip2_, "diffmisstrip2/i");
  tree_->Branch("nextmischarge1", &nextmischarge1_, "nextmischarge1/d");
  tree_->Branch("nextmischarge2", &nextmischarge2_, "nextmischarge2/d");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TkHits::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TkHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TkHits);

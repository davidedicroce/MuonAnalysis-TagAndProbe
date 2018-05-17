// -*- C++ -*-
//
// Package:    MuonAnalysis/TagAndProbe
// Class:      PixTkHits
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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
#include "CondFormats/SiPixelObjects/interface/SiPixelQuality.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "CondFormats/DataRecord/interface/SiPixelQualityFromDbRcd.h"
#include "TrackingTools/TransientTrackingRecHit/interface/InvalidTransientRecHit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDetException.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"

#include "TTree.h"

namespace {
static std::vector<std::string> sDETS{ "", "PXB", "PXF", "TIB", "TID", "TOB", "TEC" };
static std::vector<unsigned>    iDETS{ 0,   1,    2,     3,      4,     5,     6   };
static std::vector<unsigned>    iLAYS{ 0,   3,    2,     4,      3,     6,     9   };
static std::vector<std::string> sLAYS{ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" };
const float theRocWidth  = 8.1;
const float theRocHeight = 8.1;
}

class PixTkHits : public edm::EDAnalyzer {
public:
explicit PixTkHits(const edm::ParameterSet&);
~PixTkHits();

static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;

// ----------member data ---------------------------
const edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;    
const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;
const edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> stripDigiLabel_;
const edm::EDGetTokenT<edm::DetSetVector<SiStripRawDigi>> stripCommonModeLabel_;
edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
edm::EDGetTokenT<TrackerGeometry> trackergeo_;
/// Layers to debug
std::vector<std::string> layersToDebug_;
/// Track Transformer
TrackTransformer refitter_;
std::string propagator_, propagatorOpposite_;
std::string pixelQualityLabel_;
edm::ESHandle<TrackerTopology> theTrkTopo;
edm::ESHandle<Propagator> thePropagator, thePropagatorOpposite;
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
double bestch1_, bestch2_, cmode1_, cmode2_;

int missedhit_cent_1_, missedhit_cent_2_, missedhit_edg_1_, missedhit_edg_2_, validhits_pair1_, validhits_pair2_;
int cluster0modules_pair1_, cluster0modules_pair2_, cluster1modules_pair1_, cluster1modules_pair2_, cluster2modules_pair1_, cluster2modules_pair2_;
int validhit_edg_1_, validhit_edg_2_, validhit_cent_1_, validhit_cent_2_, ngoodvalid1_, ngoodvalid2_, nbadvalid1_, nbadvalid2_, diffpixel1x_,diffpixel2x_, diffpixel1y_,diffpixel2y_;
double nextvhcharge1_, nextvhcharge2_, chargevalid_1_, chargevalid_2_;

int diffmisspixel1x_, diffmisspixel2x_;
int diffmisspixel1y_, diffmisspixel2y_;
double nextmischarge1_, nextmischarge2_;
};


// constructors and destructor
PixTkHits::PixTkHits(const edm::ParameterSet& iConfig):
pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
stripDigiLabel_(consumes<edm::DetSetVector<SiStripDigi>>(iConfig.getParameter<edm::InputTag>("stripDigis"))),
stripCommonModeLabel_(consumes<edm::DetSetVector<SiStripRawDigi>>(iConfig.getParameter<edm::InputTag>("stripCommonMode"))),
tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
layersToDebug_(iConfig.getUntrackedParameter<std::vector<std::string>>("layersToDebug", std::vector<std::string>())),
refitter_(iConfig),
propagator_(iConfig.getParameter<std::string>("PropagatorAlong")),
propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),
m_lumiScalerTag_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalerTag"))),
puSummaryInfoToken_ ( consumes<std::vector<PileupSummaryInfo>>( iConfig.getParameter<edm::InputTag>("pileupInfoSummaryInputTag") ) ),
pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertex")))
{ 
//usesResource("TFileService");
}


PixTkHits::~PixTkHits()
{
}

// ------------ method called for each event  ------------
void
PixTkHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
//iSetup.get<SiPixelQualityRcd>().get(pixelQualityLabel_, thePixelQuality); 

edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelC; 
iEvent.getByToken(pixelClusterLabel_, pixelC); 
edm::Handle<edm::DetSetVector<SiStripDigi> > stripD; 
iEvent.getByToken(stripDigiLabel_, stripD); 
edm::Handle<edm::DetSetVector<SiStripRawDigi> > stripCM; 
iEvent.getByToken(stripCommonModeLabel_, stripCM); 

edm::ESHandle<TrackerGeometry> trackergeo;
iSetup.get<TrackerDigiGeometryRecord>().get(trackergeo);
const TrackerGeometry &trackergeom = *trackergeo;

Handle<MeasurementTrackerEvent> tracker;
iEvent.getByToken(tracker_, tracker);
const PxMeasurementConditionSet & pixelConds = tracker->measurementTracker().pixelDetConditions();
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
diffpixel1x_ = -1;
diffpixel1y_ = -1;
diffpixel2x_ = -1;
diffpixel2y_ = -1;
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

bestch1_ = -999.;
bestch2_ = -999.;
cmode1_ = -999.;
cmode2_ = -999.;

diffmisspixel1x_ = -1;
diffmisspixel1y_ = -1;
diffmisspixel2x_ = -1;
diffmisspixel2y_ = -1;
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
}

//std::cout << "instLumi = " << instLumi_ << std::endl;
//std::cout << "PU = " << PU_ << std::endl;

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
//	   ", lost hits: " << mutk.hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) << std::endl;
int nhits = mutk.recHitsSize();
std::vector<Trajectory> traj  = refitter_.transform(mutk);
if (traj.size() != 1) { std::cout << "Trajectory size != 1" << std::endl; continue; } 


// --------  VALID HITS ANALYSIS:  -------------------------------

const TrackingRecHit *valid = nullptr; DetId herevalid;
for (int i = 1; i < nhits; ++i) {
const TrackingRecHit *hitvalid = &* mutk.recHit(i);
if (hitvalid->isValid()) {
  if (hitvalid->getType() != TrackingRecHit::missing) { 
    valid = hitvalid;
    if (valid == nullptr) { std::cout << "valid = null pointer" << std::endl; continue; }
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
    auto temp = trackergeom.idToDet(wherevalid);
    const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*> (temp);
    const PixelTopology& topol = pixelDet->specificTopology();
    int nRowSiPixel = topol.nrows();
    int nColSiPixel = topol.ncolumns();
    //std::cout << " NRows (x) = " << nRowSiPixel << " | NCols (y) = " << nColSiPixel << std::endl;
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
    LocalVector localDir = tsosvalid.localMomentum()/tsosvalid.localMomentum().mag();
    float chargeCorr = std::abs(localDir.z());
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
      int denseIndex = pixelConds.find(wherevalid()); 
      if (denseIndex == pixelConds.nDet()) { std::cout << "Module missing in Pixel conditions set" << std::endl; continue; }
      
      //std::cout << "Analyzing module at " << wherevalid() << ", isActive? " << mdetvalid.isActive() << std::endl;
      float utrajfirstvalidX = 0, utrajlastvalidX = 0, utrajfirstvalidY = 0, utrajlastvalidY = 0;
      if (!tsosvalid.isValid()) {
	//std::cout << "  Failed to find Valid Pixel ??" << std::endl;
	continue;
      }
      if (!mdetvalid.isActive()) {
	//std::cout << "  Detector is inactive" << std::endl;
	continue;
      }
      auto cl_itervalid = pixelC->find(wherevalid);
      if (cl_itervalid == pixelC->end()) {
	//std::cout << "  ... no pixel clusters on this detid" << std::endl;
      } else {
	edmNew::DetSet<SiPixelCluster> clustersvalid = *cl_itervalid; 
	auto pixelHit = dynamic_cast<const TrackerSingleRecHit *>(valid);
	const SiPixelCluster &clusterhit = pixelHit->pixelCluster();               
	utrajfirstvalidX = clusterhit.minPixelCol();
	utrajfirstvalidY = clusterhit.minPixelRow();
	utrajlastvalidX = clusterhit.maxPixelCol();
	utrajlastvalidY = clusterhit.maxPixelRow();
	//std::cout << " X: Valid Pixel: first at " << utrajfirstvalidX << " , last at " << utrajlastvalidX << std::endl;
	//std::cout << " Y: Valid Pixel: first at " << utrajfirstvalidY << " , last at " << utrajlastvalidY << std::endl;
	int ValidPixels = clusterhit.size();
	//std::cout << "  Number of Valid Pixels : " << ValidPixels << std::endl; 
	//std::cout << "  x = " << clusterhit.x() << " y = " << clusterhit.y() << std::endl;
	//std::cout << "  size x = " << clusterhit.sizeX() << " size y = " << clusterhit.sizeY() << std::endl;
	//std::cout << "  min x = " << clusterhit.minPixelCol() << " max x = " << clusterhit.maxPixelCol() << std::endl;
	//std::cout << "  min y = " << clusterhit.minPixelRow() << " max y = " << clusterhit.maxPixelRow() << std::endl;
	if (whichpair == 1){
	  int diffpixelx = 10;
	  int diffpixely = 10;
	  double nextvhcharge1 = -999.;
	  chargevalid_1_ = clusterhit.charge()*chargeCorr;                   
	  if ( ( utrajfirstvalidX < (0.15*nRowSiPixel) || utrajlastvalidX > (nRowSiPixel - (0.15*nRowSiPixel)) ) || (utrajfirstvalidY < (0.15*nColSiPixel) || utrajlastvalidY > (nColSiPixel - 0.15*nColSiPixel) ) ){
	    validhit_edg_1_++;
	    //std::cout << "valid hit in lateral pixel, first muon" << std::endl;
	  }
	  if ( ( utrajfirstvalidX >= (0.15*nRowSiPixel) && utrajlastvalidX <= (nRowSiPixel - (0.15*nRowSiPixel)) ) && (utrajfirstvalidY >= (0.15*nColSiPixel) && utrajlastvalidY <= (nColSiPixel - (0.15*nColSiPixel)) ) ){
	    validhit_cent_1_++;
	    //std::cout << "valid hit in central pixel, first muon" << std::endl;
	  }
	  for (const SiPixelCluster &clustervalid : clustersvalid) {
	    //std::cout << "  Cluster of " << clustervalid.size() << " pixel: " << std::endl;
	    float charge = clustervalid.charge()*chargeCorr;                   
	    for (unsigned int sx = clustervalid.minPixelCol(), mx = clustervalid.maxPixelCol()+1; sx < mx; ++sx) {
	    for (unsigned int sy = clustervalid.minPixelRow(), my = clustervalid.maxPixelRow()+1; sy < my; ++sy) {
	      if ( charge == chargevalid_1_ ) {
		if ( charge < 20000) ++nbadvalid1_;
		else if ( charge >= 20000) ++ngoodvalid1_;
	      }
	      if ( ( charge != chargevalid_1_ ) && (sx < utrajfirstvalidX) && (std::abs(sx-utrajfirstvalidX) <= diffpixelx) && ( ( (sy < utrajfirstvalidY) && (std::abs(sy-utrajfirstvalidY) <= diffpixely) ) || ( (sy > utrajlastvalidY) && (std::abs(sy-utrajlastvalidY) <= diffpixely) ) ) ){
		diffpixelx = std::abs(sx-utrajfirstvalidX);
		if (sy < utrajfirstvalidY) diffpixely = std::abs(sy-utrajfirstvalidY);
		else diffpixely = std::abs(sy-utrajlastvalidY);
		nextvhcharge1 = charge;
	      } else if ( ( charge != chargevalid_1_ ) && (sx > utrajlastvalidX) && (std::abs(sx-utrajlastvalidX) <= diffpixelx) && ( ( (sy < utrajfirstvalidY) && (std::abs(sy-utrajfirstvalidY) <= diffpixely) ) || ( (sy > utrajlastvalidY) && (std::abs(sy-utrajlastvalidY) <= diffpixely) ) ) ){
		diffpixelx = std::abs(sx-utrajlastvalidX);
		if (sy < utrajfirstvalidY) diffpixely = std::abs(sy-utrajfirstvalidY);
		else diffpixely = std::abs(sy-utrajlastvalidY);
		nextvhcharge1 = charge;
	      }
	      //std::cout << "1st muon extra valid hit, min x = " << utrajfirstvalidX-sx << " | max x = " << sx-utrajlastvalidX << " | min y = " << utrajfirstvalidY-sy <<  " | max y = " << sy-utrajlastvalidY << " | charge = "<< charge << std::endl;
	    }
	    }
	    //std::cout << std::endl;
	  }
	  if ( diffpixelx < 10 && diffpixely < 10 && nextvhcharge1 != chargevalid_1_ ) {
	    diffpixel1x_ = diffpixelx;
	    diffpixel1y_ = diffpixely;
	    nextvhcharge1_ = nextvhcharge1;
	    //std::cout << " !!!!!!  EXTRA HIT !!!!!!! "<< std::endl;
	  }
	  //std::cout << " Pixel Cluster 1: found hit charge = " << chargevalid_1_ << " | diff. num. pixels EXTRA X = " << diffpixel1x_ << " | diff. num. pixels EXTRA Y = " << diffpixel1y_ << " | EXTRA charge = " << nextvhcharge1_ << std::endl;
	}
	else if (whichpair == 2){
	  int diffpixelx = 10;
	  int diffpixely = 10;
	  double nextvhcharge2 = -999.;
	  chargevalid_2_ = clusterhit.charge()*chargeCorr;
	  if ( ( utrajfirstvalidX < (0.15*nRowSiPixel) || utrajlastvalidX > (nRowSiPixel - (0.15*nRowSiPixel)) ) || (utrajfirstvalidY < (0.15*nColSiPixel) || utrajlastvalidY > (nColSiPixel - (0.15*nColSiPixel)) ) ){
	    validhit_edg_2_++;
	    //std::cout << "valid hit in lateral pixel, second muon" << std::endl;
	  }
	  if ( ( utrajfirstvalidX >= (0.15*nRowSiPixel) && utrajlastvalidX <= (nRowSiPixel - (0.15*nRowSiPixel)) ) && (utrajfirstvalidY >= (0.15*nColSiPixel) && utrajlastvalidY <= (nColSiPixel - (0.15*nColSiPixel)) ) ){
	    validhit_cent_2_++;
	    //std::cout << "valid hit in central pixel, second muon" << std::endl;
	  }
	  for (const SiPixelCluster &clustervalid : clustersvalid) {
	    //std::cout << "  Cluster of " << clustervalid.size() << " pixel: " << std::endl;
	    float charge = clustervalid.charge()*chargeCorr;                   
	    for (unsigned int sx = clustervalid.minPixelCol(), mx = clustervalid.maxPixelCol()+1; sx < mx; ++sx) {
	    for (unsigned int sy = clustervalid.minPixelRow(), my = clustervalid.maxPixelRow()+1; sy < my; ++sy) {
	      if ( charge == chargevalid_2_ ) {
		if ( charge < 20000) ++nbadvalid2_;
		else if ( charge >= 20000) ++ngoodvalid2_;
	      }
	      if ( ( charge != chargevalid_2_ ) && (sx < utrajfirstvalidX) && (std::abs(sx-utrajfirstvalidX) <= diffpixelx) && ( ( (sy < utrajfirstvalidY) && (std::abs(sy-utrajfirstvalidY) <= diffpixely) ) || ( (sy > utrajlastvalidY) && (std::abs(sy-utrajlastvalidY) <= diffpixely) ) ) ){
		nextvhcharge2 = charge;
		diffpixelx = std::abs(sx-utrajfirstvalidX);
		if (sy < utrajfirstvalidY) diffpixely = std::abs(sy-utrajfirstvalidY);
		else diffpixely = std::abs(sy-utrajlastvalidY);
	      } else if ( ( charge != chargevalid_2_ ) && (sx > utrajlastvalidX) && (std::abs(sx-utrajlastvalidX) <= diffpixelx) && ( ( (sy < utrajfirstvalidY) && (std::abs(sy-utrajfirstvalidY) <= diffpixely) ) || ( (sy > utrajlastvalidY) && (std::abs(sy-utrajlastvalidY) <= diffpixely) ) ) ){
		nextvhcharge2 = charge;
		diffpixelx = std::abs(sx-utrajlastvalidX);
		if (sy < utrajfirstvalidY) diffpixely = std::abs(sy-utrajfirstvalidY);
		else diffpixely = std::abs(sy-utrajlastvalidY);
	      }
	      //std::cout << "2st muon extra valid hit, min x = " << utrajfirstvalidX-sx << " | max x = " << sx-utrajlastvalidX << " | min y = " << utrajfirstvalidY-sy <<  " | max y = " << sy-utrajlastvalidY << " | charge = "<< nextvhcharge2 << std::endl;
	    }
	    }
	    //std::cout << std::endl;
	  }
	  if ( diffpixelx < 10 && diffpixely < 10 && nextvhcharge2 != chargevalid_2_ ) {
	    diffpixel2x_ = diffpixelx; 
	    diffpixel2y_ = diffpixely; 
	    nextvhcharge2_ = nextvhcharge2;
	    //std::cout << " !!!!!!  EXTRA HIT !!!!!!! "<< std::endl;
	  }
	  //std::cout << " Strip Cluster 2: found hit charge = " << chargevalid_2_ << " | diff. num. pixels EXTRA = " << diffpixel2x_ << " | diff. num. pixels EXTRA Y = " << diffpixel2y_ << " | EXTRA charge = " << nextvhcharge2_ << std::endl;
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
  int denseIndex = pixelConds.find(where());
  if (denseIndex == pixelConds.nDet()) { std::cout << "Module missing in pixel conditions set" << std::endl; continue; }

  //std::cout << "Analyzing module at " << where() << ", isActive? " << mdet.isActive() << std::endl;
  float utrajx = 0, utrajy = 0; bool pred = false, hascluster = false;
  //float uerr = 0;
  if (!tsosBefore.isValid()) continue;
  TrajectoryStateOnSurface tsos = thePropagator->propagate(tsosBefore, det->surface());
     if (tsos.isValid()) {  
	pred = true;
	utrajx = det->topology().measurementPosition( tsos.localPosition() ).x();
	utrajy = det->topology().measurementPosition( tsos.localPosition() ).y();
	//uerr  = std::sqrt( det->topology().measurementError( tsos.localPosition(), tsos.localError().positionError() ).uu() ); 
	//std::cout << "  Searching around pixel col" << utrajx << " +/- " << uerr << std::endl;
	//std::cout << "  Searching around pixel row" << utrajy << " +/- " << uerr << std::endl;
     } else {
	//std::cout << "  Failed to propagate??" << std::endl;
	continue;
     }
  if (!mdet.isActive()) {
    //std::cout << "  Detector is inactive" << std::endl;
    continue;
  }
  LocalVector localDir = tsos.localMomentum()/tsos.localMomentum().mag();
  float chargeCorr = std::abs(localDir.z());

  ///////////////////////////////////////////////
  auto temp = trackergeom.idToDet(where);
  const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*> (temp);
  const PixelTopology& topol = pixelDet->specificTopology();
  int nRowSiPixel = topol.nrows();
  int nColSiPixel = topol.ncolumns();
  //std::cout << " NRows (x) = " << nRowSiPixel << " | NCols (y) = " << nColSiPixel << std::endl;
  //////////////////////////////////////////////
          if (whichpair == 1) {
            double bestpixel1x = utrajx;
            double bestpixel1y = utrajy;
            if ( ( bestpixel1x < (0.15*nRowSiPixel) || bestpixel1x > (nRowSiPixel-(0.15*nRowSiPixel)) ) || (bestpixel1y < (0.15*nColSiPixel) || bestpixel1y > (nColSiPixel-(0.15*nColSiPixel)) ) ){
              //std::cout << "STRIP 1: x = " << bestpixel1x << " , y = " << bestpixel1y << std::endl;
              missedhit_edg_1_++;
              //std::cout << "missed hit in lateral strip, first muon" << std::endl;
            }
            if ( ( bestpixel1x >= (0.15*nRowSiPixel) && bestpixel1x <= (nRowSiPixel-(0.15*nRowSiPixel)) ) && (bestpixel1y >= (0.15*nColSiPixel) && bestpixel1y <= (nColSiPixel-(0.15*nColSiPixel)) ) ){
              //std::cout << "STRIP 1: x = " << bestpixel1x << " , y = " << bestpixel1y << std::endl;
              missedhit_cent_1_++;
              //std::cout << "missed hit in central strip, first muon" << std::endl;
            }
          }
          if (whichpair == 2) {
            double bestpixel2x = utrajx;
            double bestpixel2y = utrajy;
            if ( ( bestpixel2x < (0.15*nRowSiPixel) || bestpixel2x > (nRowSiPixel-(0.15*nRowSiPixel)) ) || (bestpixel2y < (0.15*nColSiPixel) || bestpixel2y > (nColSiPixel-(0.15*nColSiPixel)) ) ){
              //std::cout << "STRIP 2: x = " << bestpixel2x << " , y = " << bestpixel2y << std::endl;
              missedhit_edg_2_++;
              //std::cout << "missed hit in lateral strip, second muon" << std::endl;
            }
            if ( ( bestpixel2x >= (0.15*nRowSiPixel) && bestpixel2x <= (nRowSiPixel-(0.15*nRowSiPixel)) ) && (bestpixel2y >= (0.15*nColSiPixel) && bestpixel2y <= (nColSiPixel-(0.15*nColSiPixel)) ) ){
              //std::cout << "STRIP 2: x = " << bestpixel2x << " , y = " << bestpixel2y << std::endl;
              missedhit_cent_2_++;
              //std::cout << "missed hit in central strip, second muon" << std::endl;
            }
          }
          auto cl_iter = pixelC->find(where);
          if (cl_iter == pixelC->end()) {
            //std::cout << "  ... no pixels clusters on this detid" << std::endl;
          } else {
            edmNew::DetSet<SiPixelCluster> clusters = *cl_iter;
            if (whichpair == 1){
              int diffpixelx = 10;
              int diffpixely = 10;
              int diffmisspixelx = 10;
              int diffmisspixely = 10;
              int firstmisspixelx = -999;
              int firstmisspixely = -999;
              int lastmisspixelx = -999;
              int lastmisspixely = -999;
              double nextmischarge1 = -999.;
              for (const SiPixelCluster &cluster : clusters) {
                //std::cout << "  Cluster of " << cluster.size() << " pixels: " << std::endl;
                //std::cout << "  x = " << cluster.x() << " y = " << cluster.x() << std::endl;
                //std::cout << "  size x = " << cluster.sizeX() << " size y = " << cluster.sizeY() << std::endl;
                //std::cout << "  min x = " << cluster.minPixelCol() << " max x = " << cluster.maxPixelCol() << std::endl;
                //std::cout << "  min y = " << cluster.minPixelRow() << " max y = " << cluster.maxPixelRow() << std::endl;
                for (int sx = cluster.minPixelCol(), mx = cluster.maxPixelCol()+1; sx < mx; ++sx) {
                for (int sy = cluster.minPixelRow(), my = cluster.maxPixelRow()+1; sy < my; ++sy) {
                  float charge = cluster.charge()*chargeCorr;
		  if (pred && std::abs(sx-utrajx) <= diffpixelx && std::abs(sy-utrajy) <= diffpixely && (sqrt(pow(std::abs(sx-utrajx),2)+pow(std::abs(sy-utrajy),2))<10)) { 
                    hascluster = true;
                    //std::cout << "  Missed Cluster FOUND: charge = " << charge << " , previous charge = " << bestch1_ << std::endl;
                    //std::cout << "  x distance = " << std::abs(sx-utrajx) << " , previous x dist = " << diffpixelx << std::endl;
                    //std::cout << "  y distance = " << std::abs(sy-utrajy) << " , previous y dist = " << diffpixelx << std::endl;
                    if (std::abs(sx-utrajx) == diffpixelx && std::abs(sy-utrajy) == diffpixely && charge < bestch1_ ) {
                      //std::cout << "same distance and lower charge" << std::endl;
                      continue;
                    }
                    //std::cout << "NEW BEST MISSED CHARGE" << std::endl;
                    diffpixelx = std::abs(sx-utrajx);
                    diffpixely = std::abs(sy-utrajy);
                    firstmisspixelx = cluster.minPixelCol();
                    firstmisspixely = cluster.minPixelRow();
                    lastmisspixelx = cluster.maxPixelCol();
                    lastmisspixely = cluster.maxPixelRow();
                    bestch1_ = charge;
                  }
                }
                }
                if (hascluster) {
                  for (int sx = cluster.minPixelCol(), mx = cluster.maxPixelCol()+1; sx < mx; ++sx) {
                  for (int sy = cluster.minPixelRow(), my = cluster.maxPixelRow()+1; sy < my; ++sy) {
                    float charge = cluster.charge()*chargeCorr;
                    if ( (charge != bestch1_) && (sx < firstmisspixelx) && (std::abs(sx-firstmisspixelx) <= diffmisspixelx) && ( ( (sy < firstmisspixely) && (std::abs(sy-firstmisspixely) <= diffmisspixely) ) || ( (sy > lastmisspixely) && (std::abs(sy-lastmisspixely) <= diffmisspixely) ) ) ){
                      if ( charge < nextmischarge1 && ( (sy < firstmisspixely && std::abs(sy-firstmisspixely) == diffmisspixely) || (sy > lastmisspixely && std::abs(sy-lastmisspixely) == diffmisspixely) ) ) continue;
                      nextmischarge1 = charge;
                      diffmisspixelx = std::abs(sx-firstmisspixelx);
                      if (sy < firstmisspixely) diffmisspixely = std::abs(sy-firstmisspixely);
                      else diffmisspixely = std::abs(sy-lastmisspixely);
                    } else if ( (charge != bestch1_) && (sx > lastmisspixelx) && (std::abs(sx-lastmisspixelx) <= diffmisspixelx) && ( ( (sy < firstmisspixely) && (std::abs(sy-firstmisspixely) <= diffmisspixely) ) || ( (sy > lastmisspixely) && (std::abs(sy-lastmisspixely) <= diffmisspixely) ) ) ){
                      if ( charge < nextmischarge1 && ( (sy < firstmisspixely && std::abs(sy-firstmisspixely) == diffmisspixely) || (sy > lastmisspixely && std::abs(sy-lastmisspixely) == diffmisspixely) ) ) continue;
                      nextmischarge1 = charge;
                      diffmisspixelx = std::abs(sx-lastmisspixelx);
                      if (sy < firstmisspixely) diffmisspixely = std::abs(sy-firstmisspixely);
                      else diffmisspixely = std::abs(sy-lastmisspixely);
                    }
                  }
                  }
                  if ( diffmisspixelx < 10 && diffmisspixely < 10 && nextmischarge1 != bestch1_ ) {
                    diffmisspixel1x_ = diffmisspixelx;
                    diffmisspixel1y_ = diffmisspixely;
                    nextmischarge1_ = nextmischarge1;
                  }
                }
              }
              //std::cout << " Missed HIT: Strip Cluster Charge 1 = " << bestch1_ << std::endl;
            } 
            else if (whichpair == 2){
              int diffpixelx = 10;
              int diffpixely = 10;
              int diffmisspixelx = 10;
              int diffmisspixely = 10;
              int firstmisspixelx = -999;
              int firstmisspixely = -999;
              int lastmisspixelx = -999;
              int lastmisspixely = -999;
              double nextmischarge2 = -999.;
              for (const SiPixelCluster &cluster : clusters) {
                //std::cout << std::endl;
                //std::cout << "  Cluster of " << cluster.size() << " pixels: " << std::endl;
                for (int sx = cluster.minPixelCol(), mx = cluster.maxPixelCol()+1; sx < mx; ++sx) {
                for (int sy = cluster.minPixelRow(), my = cluster.maxPixelRow()+1; sy < my; ++sy) {
                  float charge = cluster.charge()*chargeCorr;
                  if (pred && std::abs(sx-utrajx) <= diffpixelx && std::abs(sy-utrajy) <= diffpixely && (sqrt(pow(std::abs(sx-utrajx),2)+pow(std::abs(sy-utrajy),2))<10)) {
                    hascluster = true;
                    if (std::abs(sx-utrajx) == diffpixelx && std::abs(sy-utrajy) == diffpixely && charge < bestch2_ ) continue;
                    diffpixelx = std::abs(sx-utrajx);
                    diffpixely = std::abs(sy-utrajy);
                    firstmisspixelx = cluster.minPixelCol();
                    firstmisspixely = cluster.minPixelRow();
                    lastmisspixelx = cluster.maxPixelCol();
                    lastmisspixely = cluster.maxPixelRow();
                    bestch2_ = charge;
                  }
                }
                }
                if (hascluster) {
                  for (int sx = cluster.minPixelCol(), mx = cluster.maxPixelCol()+1; sx < mx; ++sx) {
                  for (int sy = cluster.minPixelRow(), my = cluster.maxPixelRow()+1; sy < my; ++sy) {
                    float charge = cluster.charge()*chargeCorr;
                    if ( (charge != bestch2_) && (sx < firstmisspixelx) && (std::abs(sx-firstmisspixelx) <= diffmisspixelx) && ( ( (sy < firstmisspixely) && (std::abs(sy-firstmisspixely) <= diffmisspixely) ) || ( (sy > lastmisspixely) && (std::abs(sy-lastmisspixely) <= diffmisspixely) ) ) ){
                      if ( charge < nextmischarge2 && ( (sy < firstmisspixely && std::abs(sy-firstmisspixely) == diffmisspixely) || (sy > lastmisspixely && std::abs(sy-lastmisspixely) == diffmisspixely) ) ) continue;
                      nextmischarge2 = charge;
                      diffmisspixelx = std::abs(sx-firstmisspixelx);
                      if (sy < firstmisspixely) diffmisspixely = std::abs(sy-firstmisspixely);
                      else diffmisspixely = std::abs(sy-lastmisspixely);
                    } else if ( (charge != bestch2_) && (sx > lastmisspixelx) && (std::abs(sx-lastmisspixelx) <= diffmisspixelx) && ( ( (sy < firstmisspixely) && (std::abs(sy-firstmisspixely) <= diffmisspixely) ) || ( (sy > lastmisspixely) && (std::abs(sy-lastmisspixely) <= diffmisspixely) ) ) ){
                      if ( charge < nextmischarge2 && ( (sy < firstmisspixely && std::abs(sy-firstmisspixely) == diffmisspixely) || (sy > lastmisspixely && std::abs(sy-lastmisspixely) == diffmisspixely) ) ) continue;
                      nextmischarge2 = charge;
                      diffmisspixelx = std::abs(sx-lastmisspixelx);
                      if (sy < firstmisspixely) diffmisspixely = std::abs(sy-firstmisspixely);
                      else diffmisspixely = std::abs(sy-lastmisspixely);
                    }
                  }
                  }
                  if ( diffmisspixelx < 10 && diffmisspixely < 10 && nextmischarge2 != bestch2_ ) {
                    diffmisspixel2x_ = diffmisspixelx;
                    diffmisspixel2y_ = diffmisspixely;
                    nextmischarge2_ = nextmischarge2;
                  }
                }
              }
              //std::cout << " Missed HIT: Strip Cluster Charge 2 = " << bestch2_ << std::endl;
            }  
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
  if (validhits_pair1_ !=0 ||cluster0modules_pair1_ !=0 || cluster1modules_pair1_ != 0 || cluster2modules_pair1_ != 0){
    //std::cout << "1st muon| LumiBlock: " << iEvent.luminosityBlock() << " | VH: " << validhits_pair1_ << " | 0CMH: " << cluster0modules_pair1_ << " | 1CMH: " << cluster1modules_pair1_ << " | 2CMH: " << cluster2modules_pair1_ << std::endl;
  }

  if (validhits_pair2_ !=0 ||cluster0modules_pair2_ !=0 || cluster1modules_pair2_ != 0 || cluster2modules_pair2_ != 0){
    //std::cout << "2nd muon | LumiBlock: " << iEvent.luminosityBlock() << " | VH: " << validhits_pair2_ << " | 0CMH: " << cluster0modules_pair2_ << " | 1CMH: " << cluster1modules_pair2_ << " | 2CMH: " << cluster2modules_pair2_ << std::endl;
  }
  //std::cout << " Found HITS: VHC1 = " << validhit_cent_1_ << " | VHE1 = " << validhit_edg_1_ << " || VHC2 = " << validhit_cent_2_ << " | VHE2 = " << validhit_edg_2_ << std::endl;
  //std::cout << " Missed HITS: MHC1 = " << missedhit_cent_1_ << " | MHE1 = " << missedhit_edg_1_ << " || MHC2 = " << missedhit_cent_2_ << " | MHE2 = " << missedhit_edg_2_ << std::endl;
  
  if ( validhits_pair1_ != 0 || cluster0modules_pair1_ != 0 || cluster1modules_pair1_ != 0 || cluster2modules_pair1_ != 0 || validhits_pair2_ != 0 ||cluster0modules_pair2_ != 0 || cluster1modules_pair2_ != 0 || cluster2modules_pair2_ != 0 ) {
    tree_->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
PixTkHits::beginJob()
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
  tree_->Branch("diffpixel1x", &diffpixel1x_, "diffpixel1x/i");
  tree_->Branch("diffpixel1y", &diffpixel1y_, "diffpixel1y/i");
  tree_->Branch("diffpixel2x", &diffpixel2x_, "diffpixel2x/i");
  tree_->Branch("diffpixel2y", &diffpixel2y_, "diffpixel2y/i");
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
  
  tree_->Branch("bestch1", &bestch1_, "bestch1/d");
  tree_->Branch("bestch2", &bestch2_, "bestch2/d");
  tree_->Branch("cmode1", &cmode1_, "cmode1/d");
  tree_->Branch("cmode2", &cmode2_, "cmode2/d");

  tree_->Branch("diffmisspixel1x", &diffmisspixel1x_, "diffmisspixel1x/i");
  tree_->Branch("diffmisspixel1y", &diffmisspixel1y_, "diffmisspixel1y/i");
  tree_->Branch("diffmisspixel2x", &diffmisspixel2x_, "diffmisspixel2x/i");
  tree_->Branch("diffmisspixel2y", &diffmisspixel2y_, "diffmisspixel2y/i");
  tree_->Branch("nextmischarge1", &nextmischarge1_, "nextmischarge1/d");
  tree_->Branch("nextmischarge2", &nextmischarge2_, "nextmischarge2/d");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixTkHits::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PixTkHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixTkHits);

// Standard C++ includes
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>

// Boost includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "TTree.h"

// CMSSW framework includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

class WeightsTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit WeightsTreeMaker(const edm::ParameterSet&);
  ~WeightsTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;
  
  virtual void beginRun(const edm::Run &, const edm::EventSetup & ) override;
  virtual void endRun(const edm::Run & , const edm::EventSetup & ) override;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock & , const edm::EventSetup & ) override;
  virtual void endLuminosityBlock(const edm::LuminosityBlock & , const edm::EventSetup & ) override;
  
  // Tokens -- get the event weight from GenEventInfoProduct, and cross-section from LHEEventProduct 
  const edm::EDGetTokenT<LHEEventProduct>     lheInfoToken;
  const edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken;

  // Event weights
  float  wgt;
  unsigned int putrue, event, lumi;  
  TTree* tree;
  TH1F* putrue_distribution;
  TH1F* weight_distribution;

};

WeightsTreeMaker::WeightsTreeMaker(const edm::ParameterSet& iConfig): 
  lheInfoToken      (mayConsume<LHEEventProduct>    (iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken      (mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  pileupInfoToken   (mayConsume<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileupInfo"))){
  usesResource("TFileService");    
}


WeightsTreeMaker::~WeightsTreeMaker() {}

void WeightsTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace boost::algorithm;

  // Get handles to all the requisite collections
  edm::Handle<LHEEventProduct> lheInfoH;
  iEvent.getByToken(lheInfoToken, lheInfoH);  
  edm::Handle<GenEventInfoProduct> genInfoH;
  iEvent.getByToken(genInfoToken, genInfoH);  
  edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoH;
  iEvent.getByToken(pileupInfoToken,pileupInfoH);

  event = iEvent.id().event();
  lumi  = iEvent.luminosityBlock();

  // Pileup information
  putrue = 0;
  if (pileupInfoH.isValid()) {
    for (auto pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
      if (pileupInfo_iter->getBunchCrossing() == 0) // take current bunch crossing
	putrue = (unsigned int) pileupInfo_iter->getTrueNumInteractions();
    }
  }
  putrue_distribution->Fill(putrue);
  
  // Read the event weights
  wgt = 1.0;
  if(lheInfoH.isValid() and lheInfoH->weights().size() > 0)
    wgt = lheInfoH->weights()[0].wgt;
  else if(genInfoH.isValid())
    wgt = genInfoH->weight();
  weight_distribution->Fill(1,wgt);

  // Fill the tree
  tree->Fill();
    
}

void WeightsTreeMaker::beginJob() {

  edm::Service<TFileService> fs;  
  putrue_distribution = fs->make<TH1F>("putrue_distribution","putrue_distribution",99,0,99);
  weight_distribution = fs->make<TH1F>("weight_distribution","weight_distribution",1,0,2);

  tree = fs->make<TTree>("tree", "tree");
  tree->Branch("event", &event,"event/I");
  tree->Branch("lumi", &lumi,"lumi/I");
  tree->Branch("putrue", &putrue, "putrue/I");
  tree->Branch("wgt", &wgt, "wgt/F");
}

void WeightsTreeMaker::endJob() {}

void WeightsTreeMaker::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup){
}

void WeightsTreeMaker::endRun(const edm::Run & , const edm::EventSetup & ) {
}

void WeightsTreeMaker::beginLuminosityBlock(const edm::LuminosityBlock &  iLumi, const edm::EventSetup & iSetup) {  
}

void WeightsTreeMaker::endLuminosityBlock(const edm::LuminosityBlock & , const edm::EventSetup &y) {}


void WeightsTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(WeightsTreeMaker);


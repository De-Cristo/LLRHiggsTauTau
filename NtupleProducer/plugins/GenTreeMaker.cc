// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
// FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// Gen Info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
// DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

class GenTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns,edm::one::WatchLuminosityBlocks> {

public:
  explicit GenTreeMaker(const edm::ParameterSet&);
  ~GenTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool isMiniAOD;
  const edm::EDGetTokenT<LHEEventProduct> lheInfoToken;
  const edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;
  const edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken;
  const edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  genMetToken_1;
  edm::EDGetTokenT<edm::View<reco::GenMET> > genMetToken_2;
  
  TTree*   tree;
  bool     addEvtInfo;

  // event information
  uint32_t event, run, lumi;  
  float    xsec, wgtsign, wgtxsec;

  // gen jets
  std::vector<TLorentzVector> genjets;

  // gen particles
  std::vector<TLorentzVector> gens;
  std::vector<int>            gid;
  std::vector<int>            gstatus;

  // lhe particles
  std::vector<TLorentzVector> gens_lhe;
  std::vector<TLorentzVector> quark_lhe;
  std::vector<int>            gid_lhe;
  std::vector<int>            gstatus_lhe;

  // gen tau decay products
  std::vector<TLorentzVector> gentau;
  std::vector<int>            gtauid;
  std::vector<int>            gtaustatus;

  // gen met
  float genmet, genmetphi;
  
  template<typename T> 
  class PtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }    
  };  
  PtSorter<reco::GenJetRef> jetSorter;
};


GenTreeMaker::GenTreeMaker(const edm::ParameterSet& iConfig):
  isMiniAOD         (iConfig.getParameter<bool>("isMiniAOD")),
  lheInfoToken      (mayConsume<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("lheEvent"))),
  genEventInfoToken (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genEvent"))),
  genParticlesToken (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken      (consumes<std::vector<reco::GenJet> >   (iConfig.getParameter<edm::InputTag>("genJets"))),
  addEvtInfo        (iConfig.existsAs<bool>("addEventInfo")  ? iConfig.getParameter<bool>("addEventInfo")  : true),
  xsec              (iConfig.existsAs<double>("xsec")        ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.){  

  if(isMiniAOD)
    genMetToken_1   = consumes<edm::View<pat::MET> >         (iConfig.getParameter<edm::InputTag>("genMet"));
  else
    genMetToken_2   = consumes<edm::View<reco::GenMET> >     (iConfig.getParameter<edm::InputTag>("genMet"));  
  usesResource();
  usesResource("TFileService");    

}


GenTreeMaker::~GenTreeMaker() {}

void GenTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace boost::algorithm;
    
  Handle<LHEEventProduct>      lheInfoH;
  Handle<GenEventInfoProduct>  genEventInfoH;
  Handle<View<GenParticle> >   genParticlesH;
  Handle<vector<GenJet> >      genJetsH;
  Handle<View<pat::MET> >      genMetH_1;
  Handle<View<reco::GenMET> >  genMetH_2;

  if(lheInfoToken.index())
    iEvent.getByToken(lheInfoToken,lheInfoH);

  iEvent.getByToken(genEventInfoToken, genEventInfoH);
  iEvent.getByToken(genParticlesToken , genParticlesH);
  iEvent.getByToken(genJetsToken      , genJetsH);
  if(isMiniAOD)
    iEvent.getByToken(genMetToken_1, genMetH_1);
  else
    iEvent.getByToken(genMetToken_2, genMetH_2);
  
  /// event info
  event = iEvent.id().event();
  run   = iEvent.id().run();
  lumi  = iEvent.luminosityBlock();
  
  //// event weight
  wgtsign = 1.0;
  if(lheInfoH.isValid() and lheInfoH->weights().size() > 1) // there are also other variations
    wgtsign = lheInfoH->weights()[0].wgt;
  else if(genEventInfoH.isValid())
    wgtsign = genEventInfoH->weight();
  if(lheInfoH.isValid()) 
    wgtxsec = lheInfoH->originalXWGTUP();
  else 
    wgtxsec = 1.0;
 

  /// Missing energy
  if(isMiniAOD){
    if(genMetH_1.isValid()) {      
      genmet    = genMetH_1->front().genMET()->pt();
      genmetphi = genMetH_1->front().genMET()->phi();
    }
    else {
      genmet    = -99.;
      genmetphi = -99.;
    }
  }
  else{
    if(genMetH_2.isValid()){
      genmet    = genMetH_2->front().pt();
      genmetphi = genMetH_2->front().phi();
    }
    else {
      genmet    = -99.;
      genmetphi = -99.;
    }
  }

  gens.clear();
  gid.clear();
  gstatus.clear();
  gentau.clear();
  gtauid.clear();
  gtaustatus.clear();
  gens_lhe.clear();
  quark_lhe.clear();
  gid_lhe.clear();
  gstatus_lhe.clear();

  // GEN information
  if (genParticlesH.isValid()) {
    for (auto gens_iter = genParticlesH->begin(); gens_iter != genParticlesH->end(); ++gens_iter) {       
      TLorentzVector g4;

      if((gens_iter->pdgId() == 25 or gens_iter->pdgId() ==  23  or abs(gens_iter->pdgId()) == 24 or abs(gens_iter->pdgId()) == 35 or abs(gens_iter->pdgId()) == 45) 
	 and gens_iter->numberOfDaughters() > 1){ 
	g4.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	gens.push_back(g4);
	gid.push_back(gens_iter->pdgId());
	gstatus.push_back(gens_iter->status());
      }      

      if (abs(gens_iter->pdgId()) > 10 and abs(gens_iter->pdgId()) < 17 and gens_iter->fromHardProcessFinalState()) { 
	g4.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	gens.push_back(g4);
	gid.push_back(gens_iter->pdgId());
	gstatus.push_back(gens_iter->status());
      }      

      if (((abs(gens_iter->pdgId()) >= 1 and abs(gens_iter->pdgId()) <= 5) or abs(gens_iter->pdgId()) == 21) and 
	  (gens_iter->fromHardProcessFinalState() or gens_iter->status() == 23)) { 
	g4.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	gens.push_back(g4);
	gid.push_back(gens_iter->pdgId());
	gstatus.push_back(gens_iter->status());
      }      

      // look for decayed taus
      if (abs(gens_iter->pdgId()) == 15 and gens_iter->isPromptDecayed() and gens_iter->numberOfDaughters() > 1){
	g4.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	gentau.push_back(g4);
	gtauid.push_back(gens_iter->pdgId());
	gtaustatus.push_back(gens_iter->status());
      }

      if(gens_iter->isDirectHardProcessTauDecayProductFinalState()){
	g4.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	gentau.push_back(g4);
	gtauid.push_back(gens_iter->pdgId());
	gtaustatus.push_back(gens_iter->status());
      }
    }
  }

  if(lheInfoH.isValid()){

    const auto& hepeup = lheInfoH->hepeup();
    const auto& pup = hepeup.PUP;

    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int status = hepeup.ISTUP[i];
      int id     = hepeup.IDUP[i];
      TLorentzVector p4(pup[i][0], pup[i][1], pup[i][2], pup[i][3]);
      
      if(abs(id) == 25 or abs(id) == 23 or abs(id) == 24){ // Higgs, Z or W
	gens_lhe.push_back(p4);
	gstatus_lhe.push_back(status);
	gid_lhe.push_back(id);
      }
      else if(abs(id) >= 10 and abs(id) <= 17 and status == 1){ // final state leptons
	gens_lhe.push_back(p4);
	gstatus_lhe.push_back(status);
	gid_lhe.push_back(id);
      }
      else if(abs(id) >= 1 and abs(id) <= 5 and status == 1){ // final state quarks
	gens_lhe.push_back(p4);
	gstatus_lhe.push_back(status);
	gid_lhe.push_back(id);
	quark_lhe.push_back(p4);
      }
      else if(abs(id) == 21 and status == 1){ // final state gluons
	gens_lhe.push_back(p4);
	gstatus_lhe.push_back(status);
	gid_lhe.push_back(id);
	quark_lhe.push_back(p4);
      }
    }
  }

  // gen jets
  vector<GenJetRef> jets;
  if(genJetsH.isValid()) {
    for (auto jets_iter = genJetsH->begin(); jets_iter != genJetsH->end(); ++jets_iter) {
      GenJetRef jetref(genJetsH, jets_iter - genJetsH->begin());
      if (jetref.isAvailable() and jetref.isNonnull()) {
	bool isCleanJet = true;
	for(size_t igen = 0; igen < gens.size(); igen++){ // clean gen gets from gen-leptons with status 1
	  if(abs(int(gid.at(igen))) >= 11 and abs(int(gid.at(igen))) <= 16 and int(gstatus.at(igen)) == 1 and 
	     deltaR(jets_iter->eta(), jets_iter->phi(), gens.at(igen).Eta(), gens.at(igen).Phi()) < 0.4)
	    isCleanJet = false;
	}
	if(isCleanJet){
	  if(isMiniAOD) jets.push_back(jetref);
	  else if(not isMiniAOD and jetref->pt() > 8) jets.push_back(jetref);
	}
      }
    }
  }
  
  // sort in pt
  sort(jets.begin(), jets.end(), jetSorter);
  tree->Fill();    
} 

void GenTreeMaker::beginJob() {
  
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree","tree");
  
  if (addEvtInfo) { 
    tree->Branch("event", &event  , "event/i");
    tree->Branch("run"  , &run    , "run/i");
    tree->Branch("lumi" , &lumi   , "lumi/i");
  }
  
  tree->Branch("xsec"     , &xsec             , "xsec/F");
  tree->Branch("wgtsign"  , &wgtsign          , "wgtsign/F");
  tree->Branch("wgtxsec"  , &wgtxsec          , "wgtxsec/F");
   
  tree->Branch("genmet"    , &genmet    , "genmet/F");
  tree->Branch("genmetphi" , &genmetphi , "genmetphi/F");

  tree->Branch("genjets", "std::vector<TLorentzVector>", &genjets, 32000, 0);

  tree->Branch("gens"   , "std::vector<TLorentzVector>", &gens, 32000, 0);
  tree->Branch("gid"    , "std::vector<int>", &gid);
  tree->Branch("gstatus", "std::vector<int>", &gstatus);

  tree->Branch("gentau", "std::vector<TLorentzVector>", &gentau, 32000, 0);
  tree->Branch("gtauid", "std::vector<int>", &gtauid);
  tree->Branch("gtaustatus", "std::vector<int>", &gtaustatus);

  tree->Branch("gens_lhe"   , "std::vector<TLorentzVector>", &gens_lhe, 32000, 0);
  tree->Branch("quark_lhe"   , "std::vector<TLorentzVector>", &quark_lhe, 32000, 0);
  tree->Branch("gid_lhe"    , "std::vector<int>", &gid_lhe);
  tree->Branch("gstatus_lhe", "std::vector<int>", &gstatus_lhe);
 
}

void GenTreeMaker::endJob() {
}

void GenTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void GenTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void GenTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void GenTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void GenTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(GenTreeMaker);

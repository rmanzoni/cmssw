// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

//
// class declaration
//

namespace l1t {
  
  class L1TGlobalNtuplizer : public edm::EDAnalyzer {
public:
  explicit L1TGlobalNtuplizer(const edm::ParameterSet&);
  ~L1TGlobalNtuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::EDGetToken m_gtAlgToken;
  edm::EDGetToken m_emulGtAlgToken;

  int event_, run_, lumi_;
  std::vector<int> algobit_;
 
  bool m_doGtAlg;
  
  TTree* tree_; 
 
};
  
  //
  // constants, enums and typedefs
  //
  
  //
  // static data member definitions
  //
  
  //
  // constructors and destructor
  //
  L1TGlobalNtuplizer::L1TGlobalNtuplizer(const edm::ParameterSet& iConfig)
  {
  //now do what ever initialization is needed
  edm::InputTag nullTag("None");

  edm::InputTag gtAlgTag     = iConfig.getParameter<edm::InputTag>("gtAlgToken");
  m_gtAlgToken               = consumes<GlobalAlgBlkBxCollection>(gtAlgTag);
  m_doGtAlg                  = !(gtAlgTag==nullTag); 

}
  
  
  L1TGlobalNtuplizer::~L1TGlobalNtuplizer()
  {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
  
  
  //
  // member functions
  //
  
  // ------------ method called for each event  ------------
  void L1TGlobalNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
  {

  using namespace edm;
 
  run_   = iEvent.id().run();
  lumi_  = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();
  algobit_.clear();
   
  //Algorithm Bits
  if (m_doGtAlg) {
    Handle< BXVector<GlobalAlgBlk> > algs;
    iEvent.getByToken(m_gtAlgToken,algs);


    for ( auto itr = algs->begin(0); itr != algs->end(0); ++itr ) {
      for(int algBit=0; algBit<211; algBit++) {  //Fix Me: Should access size of algo vector...need method in GlobalAlgBlk class
        if(itr->getAlgoDecisionFinal(algBit)) {
          algobit_.push_back(algBit);
        }  
      }
    }

    
//     for ( int ibx=algs->getFirstBX(); ibx<=algs->getLastBX(); ++ibx) {
// 
//       for ( auto itr = algs->begin(ibx); itr != algs->end(ibx); ++itr ) {
//            
//         for(int algBit=0; algBit<211; algBit++) {  //Fix Me: Should access size of algo vector...need method in GlobalAlgBlk class
//           if(itr->getAlgoDecisionFinal(algBit)) {
//             algobit_.push_back(algBit);
//           }  
//         }
//       }
//     }
  }
  
  tree_ -> Fill();
  
}
  
  
  
  // ------------ method called once each job just before starting event loop  ------------
  void L1TGlobalNtuplizer::beginJob()
  {
  edm::Service<TFileService> outfile_;

  TH1::SetDefaultSumw2();
  
  tree_ = outfile_ -> make<TTree>("tree", "tree");
  tree_ -> Branch("event"  , &event_  );
  tree_ -> Branch("run"    , &run_    );
  tree_ -> Branch("lumi"   , &lumi_   );
  tree_ -> Branch("algobit", &algobit_);
}
  
  // ------------ method called once each job just after ending the event loop  ------------
  void L1TGlobalNtuplizer::endJob() 
  {
  }
  
  // ------------ method called when starting to processes a run  ------------
  /*
  void 
  L1TGlobalNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
  */
  
  // ------------ method called when ending the processing of a run  ------------
  /*
  void 
  L1TGlobalNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
  */
  
  // ------------ method called when starting to processes a luminosity block  ------------
  /*
  void 
  L1TGlobalNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
  */
  
  // ------------ method called when ending the processing of a luminosity block  ------------
  /*
  void 
  L1TGlobalNtuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
  */
  
  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void L1TGlobalNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

}

using namespace l1t;

//define this as a plug-in
DEFINE_FWK_MODULE(L1TGlobalNtuplizer);

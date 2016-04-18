
// system include files
#include <boost/shared_ptr.hpp>

// user include files

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"

#include "DataFormats/L1Trigger/interface/Tau.h"

//
// class declaration
//

using namespace l1t;

class FilterL1TTaus : public edm::EDProducer {
  public:
    explicit FilterL1TTaus(const edm::ParameterSet& ps);
    ~FilterL1TTaus();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  private:
    virtual void beginJob() override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
  
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
    // ----------member data ---------------------------
  
    // input token
    edm::EDGetToken tauToken_;
    float minPt_, maxEta_;
    int iso_, maxNumber_;
    bool verbose_;
};



FilterL1TTaus::FilterL1TTaus(const edm::ParameterSet& ps) {
  // register what you produce
  produces<l1t::TauBxCollection> ();

  // register what you consume and keep token for later access:
  tauToken_  = consumes<l1t::TauBxCollection>(ps.getParameter<edm::InputTag>("L1TauCollection"));
  minPt_     =                               (ps.getParameter<double>       ("minPt"          ));
  maxEta_    =                               (ps.getParameter<double>       ("maxEta"         ));
  iso_       =                               (ps.getParameter<int>          ("iso"            ));
  maxNumber_ =                               (ps.getParameter<int>          ("maxNumber"      ));
  verbose_   =                               (ps.getParameter<bool>         ("verbose"        ));
}

FilterL1TTaus::~FilterL1TTaus() {

}

// ------------ method called to produce the data  ------------
void
FilterL1TTaus::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace l1t;
  
  //inputs
  Handle< l1t::TauBxCollection > itaus;
  iEvent.getByToken(tauToken_,itaus);

  int bxFirst = itaus->getFirstBX();
  int bxLast  = itaus->getLastBX ();

  //outputs
  std::auto_ptr<l1t::TauBxCollection> ftaus (new l1t::TauBxCollection(0, bxFirst, bxLast));

  // loop over BX
  for(int ibx = bxFirst; ibx < bxLast+1; ++ibx) {
    int i = 1;
    std::cout << "bx " << ibx;  
    for(auto tau = itaus->begin(ibx); tau != itaus->end(ibx); ++tau, ++i){
      
      /*
      std::cout << "\ttau pt  " << tau->pt()  
                << "\ttau iso " << tau->hwIso() 
                << "\ttau eta " << tau->eta() 
                << "\ttau #   " << i            << std::endl;
      //*/
      if (tau->pt()  < minPt_)        continue;
      if (tau->hwIso() < iso_  )      continue;
      if (fabs(tau->eta()) > maxEta_) continue;
      if (i > maxNumber_)             continue;

      ftaus->push_back(ibx, *tau);
    }
  }
  iEvent.put(ftaus);
}

// ------------ method called once each job just before starting event loop  ------------
void
FilterL1TTaus::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
FilterL1TTaus::endJob() {
}

// ------------ method called when starting to processes a run  ------------
// void
// FilterL1TTaus::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
// {
// }


// ------------ method called when ending the processing of a run  ------------
// void
// FilterL1TTaus::endRun(edm::Run const&, edm::EventSetup const&)
// {
// }

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
FilterL1TTaus::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup cons
t&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
FilterL1TTaus::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&
)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FilterL1TTaus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterL1TTaus);
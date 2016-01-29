/////////////////////////////////////////////////////////////////////////
//
// METAnalyzer
// ------------------
//
//                        05/2015    Sébastien Brochet <sebastien.brochet@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#include "JMEAnalysis/JMEValidator/interface/METAnalyzer.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include <vector>
#include <iostream>
#include <regex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
METAnalyzer::METAnalyzer(const edm::ParameterSet& iConfig):
  JME::Analyzer(iConfig){
  if (iConfig.existsAs<edm::InputTag>("srcJet"))
    srcJet_ = iConfig.getParameter<edm::InputTag>("srcJet");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input jet collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcJetPF"))
    srcJetPF_ = iConfig.getParameter<edm::InputTag>("srcJetPF");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input PF jet collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcVertex"))
    srcVertex_ = iConfig.getParameter<edm::InputTag>("srcVertex");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input vertex collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcZboson"))
    srcZboson_ = iConfig.getParameter<edm::InputTag>("srcZboson");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input Zboson collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcLeptons"))
    srcLeptons_ = iConfig.getParameter<edm::InputTag>("srcLeptons");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input Leptons collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenJets"))
    srcGenJets_ = iConfig.getParameter<edm::InputTag>("srcGenJets");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input gen jet collection not given \n";

//  if (iConfig.existsAs<edm::InputTag>("srcGenJetsCleaned"))//    srcGenJetsCleaned_ = iConfig.getParameter<edm::InputTag>("srcGenJetsCleaned");
//    srcGenJetsCleaned_ = iConfig.getParameter<edm::InputTag>("srcGenJetsCleaned");
//  else throw cms::Exception("Configuration")<<"[METAnalyzer] input gen jet cleaned collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenMet"))
    srcGenMet_ = iConfig.getParameter<edm::InputTag>("srcGenMet");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input gen MET location not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenParticles"))
    srcGenParticles_ = iConfig.getParameter<edm::InputTag>("srcGenParticles");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input gen particle collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenEventInfo"))
    srcGenEventInfo_ = iConfig.getParameter<edm::InputTag>("srcGenEventInfo");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input no gen event info \n";

  if (iConfig.existsAs<double>("dRgenMatching"))
    dRgenMatching_ = iConfig.getParameter<double>("dRgenMatching");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input dR for gen jet matching not given \n";

  if (iConfig.existsAs<bool>("isMC"))
    isMC_ = iConfig.getParameter<bool>("isMC");
  else throw cms::Exception("Configuration")<<"[METAnalyzer] input isMC \n";

  
  if(!(srcJet_ == edm::InputTag("")))
    srcJetToken_ = consumes<pat::JetCollection>(srcJet_);

  if(!(srcJetPF_ == edm::InputTag("")))
    srcJetPFToken_ = consumes<pat::JetCollection>(srcJetPF_);

  if(!(srcZboson_ == edm::InputTag("")))
    srcZbosonToken_ = consumes<std::vector<reco::Particle>>(srcZboson_);

  if(!(srcLeptons_ == edm::InputTag("")))
    srcLeptonsToken_ = consumes<reco::CandidateView>(srcLeptons_);

  if(!(srcVertex_ == edm::InputTag("")))
    srcVertexToken_ = consumes<reco::VertexCollection>(srcVertex_);

  if(!(srcGenJets_ == edm::InputTag("")) and isMC_)
    srcGenJetsToken_ = consumes<reco::GenJetCollection>(srcGenJets_);

//  if(!(srcGenJetsCleaned_ == edm::InputTag("")) and isMC_)
//    srcGenJetsCleanedToken_ = consumes<pat::JetCollection>(srcGenJetsCleaned_);

  if(!(srcGenMet_ == edm::InputTag("")) and isMC_)
    srcGenMetToken_ = consumes<pat::METCollection>(srcGenMet_);

  if(!(srcGenParticles_ == edm::InputTag("")) and isMC_)
    srcGenParticlesToken_ = consumes<reco::GenParticleCollection>(srcGenParticles_);

  if(!(srcGenEventInfo_ == edm::InputTag("")) and isMC_)
    srcGenEventInfoToken_ = consumes<GenEventInfoProduct>(srcGenEventInfo_);


  // save all recoils as specified in config 
  srcRecoilTags_ = iConfig.getParameter<vInputTag>("srcRecoils");
  for(vInputTag::const_iterator it=srcRecoilTags_.begin(); it!=srcRecoilTags_.end(); it++)
  {
    srcRecoils_.push_back( consumes<pat::METCollection >( *it ) );
    // set references to write out n-tuple
    for(auto recoilAttribute : recoilAttributes_)
    {
      recoilReferences_.push_back(tree[it->instance() + "_" + recoilAttribute].write<float>());
    }
  }
  // save all METs as specified in config 
  srcMETTags_ = iConfig.getParameter<vInputTag>("srcMETs");
  for(vInputTag::const_iterator it=srcMETTags_.begin(); it!=srcMETTags_.end(); it++)
  {
    srcMETs_.push_back( consumes<pat::METCollection >( *it ) );
    // set references to write out n-tuple
    for(auto METAttribute : METAttributes_)
    {
      METReferences_.push_back(tree[it->label() + "_" + METAttribute].write<float>());
    }
  }


}


//______________________________________________________________________________
METAnalyzer::~METAnalyzer(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________

//______________________________________________________________________________
void METAnalyzer::analyze(const edm::Event& iEvent,
			     const edm::EventSetup& iSetup){





  // take gen level info weight  
  int ijet = 0;

  if(isMC_){
    edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
    iEvent.getByToken(srcGenEventInfoToken_,GenEventInfoHandle);
    eventMCWeight = GenEventInfoHandle->weight()/fabs(GenEventInfoHandle->weight());

    // find generator level Z kinematics
    edm::Handle<reco::GenParticleCollection> GenParticlesHandle;
    iEvent.getByToken(srcGenParticlesToken_, GenParticlesHandle);
    reco::GenParticle GenBoson;

    for(auto aGenParticle : *GenParticlesHandle){
      if(aGenParticle.pdgId() == 23 or abs(aGenParticle.pdgId()) == 24 or aGenParticle.pdgId() == 22){
	if(aGenParticle.numberOfDaughters() !=2) continue;
	bool goodBoson = false;
	int numberOfLeptonDaughter   = 0;
	int leptonDaughterPdgId      = 0;
	int numberOfNeutrinoDaughter = 0;
	for(unsigned int i0 = 0; i0 < aGenParticle.numberOfDaughters(); i0++) {
	  const reco::GenParticle *daughter = aGenParticle.daughterRef(i0).get();
	  if(abs(daughter->pdgId()) == 11 or abs(daughter->pdgId()) == 13 or abs(daughter->pdgId()) == 15){
	    numberOfLeptonDaughter ++;
	    leptonDaughterPdgId = abs(daughter->pdgId());
	}
	  else if(abs(daughter->pdgId()) == 12 or abs(daughter->pdgId()) == 14 or abs(daughter->pdgId()) == 16)
	    numberOfNeutrinoDaughter++;
	}
	
	if(numberOfLeptonDaughter == 1 and numberOfNeutrinoDaughter == 1)
	  goodBoson = true;

	if(numberOfLeptonDaughter == 2 and numberOfNeutrinoDaughter == 0)
	  goodBoson = true;
	
	if(goodBoson == false) continue;

	GenBoson_Pt_  = aGenParticle.pt();
	GenBoson_Eta_ = aGenParticle.eta();
	GenBoson_Phi_ = aGenParticle.phi();
	GenBoson_M_   = aGenParticle.mass();
	GenBoson.setP4(aGenParticle.p4());
	GenBoson_daughter_ = leptonDaughterPdgId;
      }
    }
  

    // store gen jets and gen jets cleaned info
    edm::Handle<reco::GenJetCollection> GenJetsHandle;
    iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
    
//    edm::Handle<pat::JetCollection> GenJetsCleanedHandle;
//    iEvent.getByToken(srcGenJetsCleanedToken_, GenJetsCleanedHandle);
    
    NGenJets_        = GenJetsHandle->size();
    //NGenJetsCleaned_ = GenJetsHandle->size();

    for( auto GenJet : *GenJetsHandle){     
      if(ijet == 0){
	GenLeadingJet_Pt_  = GenJet.pt();
	GenLeadingJet_Eta_ = GenJet.eta();
	GenLeadingJet_Phi_ = GenJet.phi();
	GenLeadingJet_M_   = GenJet.mass();
	ijet++;
	continue;
      }
      else if(ijet == 1){
	GenTrailingJet_Pt_  = GenJet.pt();
	GenTrailingJet_Eta_ = GenJet.eta();
	GenTrailingJet_Phi_ = GenJet.phi();
	GenTrailingJet_M_   = GenJet.mass();
	break;
      }
    }
    
    ijet = 0;
    
//    for( auto GenJet : *GenJetsCleanedHandle){     
 //     if(ijet == 0){
//	GenLeadingJetCleaned_Pt_  = GenJet.pt();
//	GenLeadingJetCleaned_Eta_ = GenJet.eta();
//	GenLeadingJetCleaned_Phi_ = GenJet.phi();
//	GenLeadingJetCleaned_M_   = GenJet.mass();
//	ijet++;
//	continue;
 //     }
//      else if(ijet == 1){
//	GenTrailingJetCleaned_Pt_  = GenJet.pt();
//	GenTrailingJetCleaned_Eta_ = GenJet.eta();
//	GenTrailingJetCleaned_Phi_ = GenJet.phi();
//	GenTrailingJetCleaned_M_   = GenJet.mass();
//	break;
 //     }
  //  }


    edm::Handle<std::vector<pat::MET>> GenMetHandle;
    iEvent.getByToken(srcGenMetToken_, GenMetHandle);
    
    const reco::GenMET* genMET_ = GenMetHandle->at(0).genMET();
    GenRecoil_sumEt_ = genMET_->sumEt();
    
    pat::MET genRecoil;
    genRecoil.setP4(-GenBoson.p4()-genMET_->p4());
    GenRecoil_Pt_    = genRecoil.pt();
    GenRecoil_Phi_   = genRecoil.phi();
    RecoilVec.SetMagPhi(GenRecoil_Pt_,reco::deltaPhi(GenRecoil_Phi_,GenBoson.phi()));
    GenRecoil_PerpZ_ = RecoilVec.Py();
    GenRecoil_LongZ_ = RecoilVec.Px();
    
    GenBosonVec.SetMagPhi(GenBoson.pt(),reco::deltaPhi(GenBoson.phi(),GenRecoil_Phi_));
    GenBoson_PerpU_ = GenBosonVec.Py();
    GenBoson_LongU_ = GenBosonVec.Px() - GenBoson.pt();
  }


  // store jet info
  edm::Handle<std::vector<pat::Jet>> jetHandle;
  iEvent.getByToken(srcJetToken_, jetHandle);

  // store jet PF info
  edm::Handle<std::vector<pat::Jet>> jetPFHandle;
  iEvent.getByToken(srcJetPFToken_, jetPFHandle);

  NCleanedJets_  = jetHandle->size();
  NCleanedJetsPF_  = jetPFHandle->size();

  ijet = 0;
  for( auto jet : *jetHandle){     
    if(ijet == 0){
      LeadingJet_Pt_  = jet.pt();
      LeadingJet_Eta_ = jet.eta();
      LeadingJet_Phi_ = jet.phi();
      LeadingJet_M_   = jet.mass();
      ijet++;
      continue;
    }
    else if(ijet == 1){
      TrailingJet_Pt_  = jet.pt();
      TrailingJet_Eta_ = jet.eta();
      TrailingJet_Phi_ = jet.phi();
      TrailingJet_M_   = jet.mass();
      break;
    }
  }


  ijet = 0;
  for( auto jet : *jetPFHandle){     
    if(ijet == 0){
      LeadingJetPF_Pt_  = jet.pt();
      LeadingJetPF_Eta_ = jet.eta();
      LeadingJetPF_Phi_ = jet.phi();
      LeadingJetPF_M_   = jet.mass();
      ijet++;
      continue;
    }
    else if(ijet == 1){
      TrailingJetPF_Pt_  = jet.pt();
      TrailingJetPF_Eta_ = jet.eta();
      TrailingJetPF_Phi_ = jet.phi();
      TrailingJetPF_M_   = jet.mass();
      break;
    }
  }

  // vertexes
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(srcVertexToken_, vertexHandle);
  NVertex_  = vertexHandle->size();

  // Zboson
  edm::Handle<std::vector<reco::Particle>> BosonHandle;
  iEvent.getByToken(srcZbosonToken_, BosonHandle);

  const reco::Particle& Boson = BosonHandle->at(0);
  Boson_Pt_  = Boson.pt();
  Boson_Eta_ = Boson.eta();
  Boson_Phi_ = Boson.phi();
  Boson_M_   = Boson.p4().M();
  Boson_daughter_ = Boson.pdgId();
  BosonVec.SetMagPhi(Boson.pt(), TVector2::Phi_mpi_pi(Boson.phi()));

  rotationMatrix(0,0) = rotationMatrix(1,1) = std::cos( GenBosonVec.Phi());
  rotationMatrix(0,1) =   std::sin( GenBosonVec.Phi());
  rotationMatrix(1,0) = - std::sin( GenBosonVec.Phi());
  // Leptons
  edm::Handle<reco::CandidateView> LeptonHandle;
  iEvent.getByToken(srcLeptonsToken_, LeptonHandle);

  int iLep = 0;
  for(reco::CandidateView::const_iterator lepton = LeptonHandle->begin(); lepton != LeptonHandle->end(); ++lepton){

    if(iLep == 0){
      LeadingLepton_Pt_  = (*lepton).pt();
      LeadingLepton_Eta_ = (*lepton).eta();
      LeadingLepton_Phi_ = (*lepton).phi();
      LeadingLepton_M_   = (*lepton).p4().M();
    }
    else if(iLep == 1){

      if((*lepton).pt() > LeadingLepton_Pt_){

	TrailingLepton_Pt_  = LeadingLepton_Pt_;
	TrailingLepton_Eta_ = LeadingLepton_Eta_;
	TrailingLepton_Phi_ = LeadingLepton_Phi_;
	TrailingLepton_M_   = LeadingLepton_M_;
	
	LeadingLepton_Pt_  = (*lepton).pt();
	LeadingLepton_Eta_ = (*lepton).eta();
	LeadingLepton_Phi_ = (*lepton).phi();
	LeadingLepton_M_   = (*lepton).p4().M();
      }
      else {
	TrailingLepton_Pt_  = (*lepton).pt();
	TrailingLepton_Eta_ = (*lepton).eta();
	TrailingLepton_Phi_ = (*lepton).phi();
	TrailingLepton_M_   = (*lepton).p4().M();
      }
    }
    iLep++;
  }
    

  // save recoils to Map
  size_t recoilAttributeCounter = 0;
  float refRecoil = 1;
  bool first = true;
  for ( std::vector<edm::EDGetTokenT<pat::METCollection> >::const_iterator srcRecoil = srcRecoils_.begin(); srcRecoil != srcRecoils_.end(); ++srcRecoil )
  {
    //get inputs
    edm::Handle<pat::METCollection> vRecoil;
    iEvent.getByToken(*srcRecoil, vRecoil);
    assert((*vRecoil).size() == 1 );
    auto recoil = (*vRecoil)[0];
    recoilReferences_[recoilAttributeCounter++].get() = recoil.sumEt();
    //if(srcRecoil == srcRecoils_.begin())
    if(first)
    {
      refRecoil = recoil.sumEt();
      first = false;
    }
    recoilReferences_[recoilAttributeCounter++].get() = recoil.sumEt() / refRecoil;
    recoilReferences_[recoilAttributeCounter++].get() = recoil.p4().Pt();
    recoilReferences_[recoilAttributeCounter++].get() = recoil.p4().Phi();

    RecoilVec.SetMagPhi(recoil.p4().Pt(), reco::deltaPhi(recoil.p4().Phi(),Boson_Phi_));

    TVector2 recoilTV2( recoil.p4().Pt(), recoil.p4().Phi());
    TVector2 recoilOnRecoBoson = recoilTV2.Rotate( - BosonVec.Phi());
    recoilReferences_[recoilAttributeCounter++].get() = RecoilVec.Py();
    recoilReferences_[recoilAttributeCounter++].get() = RecoilVec.Px();

    TVector2 recoilOnGenBoson = recoilTV2.Rotate( - GenBosonVec.Phi());
    recoilReferences_[recoilAttributeCounter++].get() = recoilOnGenBoson.Py();
    recoilReferences_[recoilAttributeCounter++].get() = recoilOnGenBoson.Px() - recoil.p4().Pt();

    // todo: check if to rotate matrix to different reference frame
    ROOT::Math::SMatrix<double,2> rotatedMatrix = recoil.getSignificanceMatrix();//rotationMatrix * recoil.getSignificanceMatrix();
    recoilReferences_[recoilAttributeCounter++].get() = rotatedMatrix(1,1);
    recoilReferences_[recoilAttributeCounter++].get() = rotatedMatrix(0,0);
    recoilReferences_[recoilAttributeCounter++].get() = rotatedMatrix(0,1);
    recoilReferences_[recoilAttributeCounter++].get() = rotatedMatrix(1,0);
  }
  // save METs to Map
  size_t METAttributeCounter = 0;
  for ( std::vector<edm::EDGetTokenT<pat::METCollection> >::const_iterator srcMET = srcMETs_.begin(); srcMET != srcMETs_.end(); ++srcMET )
  {
    //get inputs
    edm::Handle<pat::METCollection> vMET;
    iEvent.getByToken(*srcMET, vMET);
    assert((*vMET).size() == 1 );
    auto Met = (*vMET)[0];
    METReferences_[METAttributeCounter++].get() = Met.sumEt();
    METReferences_[METAttributeCounter++].get() = Met.p4().Pt();
    METReferences_[METAttributeCounter++].get() = Met.p4().Phi();
    // todo: check if to rotate matrix to different reference frame
    ROOT::Math::SMatrix<double,2> rotatedMatrix = Met.getSignificanceMatrix(); //rotationMatrix * Met.getSignificanceMatrix();
    METReferences_[METAttributeCounter++].get() = rotatedMatrix(1,1);
    METReferences_[METAttributeCounter++].get() = rotatedMatrix(0,0);
    METReferences_[METAttributeCounter++].get() = rotatedMatrix(0,1);
    METReferences_[METAttributeCounter++].get() = rotatedMatrix(1,0);
  }

  // dump all jet info
  AllJets_Pt_.clear();
  AllJets_Eta_.clear();
  AllJets_Phi_.clear();
  AllJets_M_.clear();
  GenMatchedJets_Pt_.clear();
  GenMatchedJets_Eta_.clear();
  GenMatchedJets_Phi_.clear();
  GenMatchedJets_M_.clear();

  for(auto jet : *jetHandle){
    AllJets_Pt_.push_back(jet.pt());
    AllJets_Eta_.push_back(jet.eta());
    AllJets_Phi_.push_back(jet.phi());
    AllJets_M_.push_back(jet.mass());
    if(isMC_){
      edm::Handle<reco::GenJetCollection> GenJetsHandle;
      iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
      for(auto GenJet : *GenJetsHandle){
	if(reco::deltaR(jet.eta(),jet.phi(),GenJet.eta(),GenJet.phi()) < dRgenMatching_){
	  GenMatchedJets_Pt_.push_back(GenJet.pt());
	  GenMatchedJets_Eta_.push_back(GenJet.eta());
	  GenMatchedJets_Phi_.push_back(GenJet.phi());
	  GenMatchedJets_M_.push_back(GenJet.mass());
	  break;
	}
      }
    }
  }

  NGenMatchedJets_ = GenMatchedJets_Pt_.size();

  // dump all PF jet info
  AllJetsPF_Pt_.clear();
  AllJetsPF_Eta_.clear();
  AllJetsPF_Phi_.clear();
  AllJetsPF_M_.clear();
  GenMatchedJetsPF_Pt_.clear();
  GenMatchedJetsPF_Eta_.clear();
  GenMatchedJetsPF_Phi_.clear();
  GenMatchedJetsPF_M_.clear();

  for(auto jet : *jetPFHandle){
    AllJetsPF_Pt_.push_back(jet.pt());
    AllJetsPF_Eta_.push_back(jet.eta());
    AllJetsPF_Phi_.push_back(jet.phi());
    AllJetsPF_M_.push_back(jet.mass());
    if(isMC_){
      edm::Handle<reco::GenJetCollection> GenJetsHandle;
      iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
      for(auto GenJet : *GenJetsHandle){
	if(reco::deltaR(jet.eta(),jet.phi(),GenJet.eta(),GenJet.phi()) < dRgenMatching_){
	  GenMatchedJetsPF_Pt_.push_back(GenJet.pt());
	  GenMatchedJetsPF_Eta_.push_back(GenJet.eta());
	  GenMatchedJetsPF_Phi_.push_back(GenJet.phi());
	  GenMatchedJetsPF_M_.push_back(GenJet.mass());
	  break;
	}
      }
    }
  }

  NGenMatchedJetsPF_ = GenMatchedJetsPF_Pt_.size();

  tree.fill();
  
}

pat::MET METAnalyzer::getUncorrectedRecoil(const pat::MET& input)
{
  // turn recoil in direction of MET
  TVector2 helperVector2(input.px(), input.py());
  helperVector2 = helperVector2.Rotate(M_PI);
  reco::Candidate::LorentzVector helperVector4( helperVector2.Px(), helperVector2.Py(), 0, input.sumEt());
  pat::MET negRecoil(input);
  negRecoil.setP4(helperVector4);

  //retrieve uncorrected values and rotate them back by PI
  //helperVector2.SetMagPhi( negRecoil.uncorrectedPt(), negRecoil.uncorrectedPhi());
  helperVector2 = helperVector2.Rotate(- M_PI);
  helperVector4.SetXYZT(helperVector2.Px(), helperVector2.Py(), 0, negRecoil.sumEt());
  pat::MET uncorrectedRecoil(input);
  uncorrectedRecoil.setP4(helperVector4);
  //uncorrectedRecoil.setSumEt(negRecoil.uncorrectedSumEt());
  return uncorrectedRecoil;
}
////////////////////////////////////////////////////////////////////////////////
// define METAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METAnalyzer);

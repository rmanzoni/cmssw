#ifndef HLTrigger_Muon_HLTTriMuonIsolation_h
#define HLTrigger_Muon_HLTTriMuonIsolation_h

#include <iostream>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


class HLTTriMuonIsolation : public edm::global::EDProducer<> {
    public:
        explicit HLTTriMuonIsolation(const edm::ParameterSet& iConfig);
        ~HLTTriMuonIsolation();
        virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> L3MuonsToken_        ;
        const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> AllMuonsToken_       ;
        const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> L3DiMuonsFilterToken_;
        const edm::EDGetTokenT<reco::TrackCollection>                IsoTracksToken_      ;

        edm::Handle<reco::RecoChargedCandidateCollection> L3MuCands           ;
        edm::Handle<trigger::TriggerFilterObjectWithRefs> L3DiMuonsFilterCands;
        edm::Handle<reco::RecoChargedCandidateRef>        PassedL3Muons       ;
        edm::Handle<reco::RecoChargedCandidateCollection> AllMuCands          ;
        edm::Handle<reco::TrackCollection>                IsoTracks           ;

        static bool ptComparer(const reco::RecoChargedCandidate mu_1, const reco::RecoChargedCandidate mu_2) { return mu_1.pt() > mu_2.pt(); }
        
        const double Muon1PtCut_      ;
        const double Muon2PtCut_      ;
        const double Muon3PtCut_      ;
        const double TriMuonPtCut_    ;
        const double TriMuonEtaCut_   ;
        const double ChargedRelIsoCut_;
        const double ChargedAbsIsoCut_;
        const double IsoConeSize_     ;
        const double MinTriMuonMass_  ;
        const double MaxTriMuonMass_  ;
        const int    TriMuonAbsCharge_;
        const double MaxDZ_           ;
        const bool   EnableRelIso_    ;
        const bool   EnableAbsIso_    ;
};

HLTTriMuonIsolation::HLTTriMuonIsolation(const edm::ParameterSet& iConfig):
    L3MuonsToken_        (consumes<reco::RecoChargedCandidateCollection> (iConfig.getParameter<edm::InputTag>("L3MuonsSrc"        ))),
    AllMuonsToken_       (consumes<reco::RecoChargedCandidateCollection> (iConfig.getParameter<edm::InputTag>("AllMuonsSrc"       ))),
    L3DiMuonsFilterToken_(consumes<trigger::TriggerFilterObjectWithRefs> (iConfig.getParameter<edm::InputTag>("L3DiMuonsFilterSrc"))),
    IsoTracksToken_      (consumes<reco::TrackCollection>                (iConfig.getParameter<edm::InputTag>("IsoTracksSrc"      ))),
    Muon1PtCut_                                                          (iConfig.getParameter<double>       ("Muon1PtCut"        )) ,
    Muon2PtCut_                                                          (iConfig.getParameter<double>       ("Muon2PtCut"        )) ,
    Muon3PtCut_                                                          (iConfig.getParameter<double>       ("Muon3PtCut"        )) ,
    TriMuonPtCut_                                                        (iConfig.getParameter<double>       ("TriMuonPtCut"      )) ,
    TriMuonEtaCut_                                                       (iConfig.getParameter<double>       ("TriMuonEtaCut"     )) ,
    ChargedRelIsoCut_                                                    (iConfig.getParameter<double>       ("ChargedRelIsoCut"  )) ,
    ChargedAbsIsoCut_                                                    (iConfig.getParameter<double>       ("ChargedAbsIsoCut"  )) ,
    IsoConeSize_                                                         (iConfig.getParameter<double>       ("IsoConeSize"       )) ,
    MinTriMuonMass_                                                      (iConfig.getParameter<double>       ("MinTriMuonMass"    )) ,
    MaxTriMuonMass_                                                      (iConfig.getParameter<double>       ("MaxTriMuonMass"    )) , 
    TriMuonAbsCharge_                                                    (iConfig.getParameter<int>          ("TriMuonAbsCharge"  )) ,
    MaxDZ_                                                               (iConfig.getParameter<double>       ("MaxDZ"             )) , 
    EnableRelIso_                                                        (iConfig.getParameter<bool>         ("EnableRelIso"      )) ,
    EnableAbsIso_                                                        (iConfig.getParameter<bool>         ("EnableAbsIso"      )) 
{
    //register products
    produces<reco::CompositeCandidateCollection>("Taus");
    produces<reco::CompositeCandidateCollection>("SelectedTaus");
}

HLTTriMuonIsolation::~HLTTriMuonIsolation(){ }

void
HLTTriMuonIsolation::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & iSetup) const
{
    std::unique_ptr<reco::CompositeCandidateCollection> Taus        (new reco::CompositeCandidateCollection);
    std::unique_ptr<reco::CompositeCandidateCollection> SelectedTaus(new reco::CompositeCandidateCollection);

    // Get the L3 muon candidates
    edm::Handle<reco::RecoChargedCandidateCollection> L3MuCands;
    iEvent.getByToken(L3MuonsToken_, L3MuCands);

    // Get the L3 muon candidates that passed the filter
    edm::Handle<trigger::TriggerFilterObjectWithRefs> L3DiMuonsFilterCands;
    iEvent.getByToken(L3DiMuonsFilterToken_, L3DiMuonsFilterCands);
    
    std::vector<reco::RecoChargedCandidateRef> PassedL3Muons;
    L3DiMuonsFilterCands->getObjects(trigger::TriggerMuon, PassedL3Muons);
  
    // Get the Trk + L3 muon candidates (after merging)
    edm::Handle<reco::RecoChargedCandidateCollection> AllMuCands;
    iEvent.getByToken(AllMuonsToken_, AllMuCands);
    
    // Get iso tracks
    edm::Handle<reco::TrackCollection> IsoTracks;
    iEvent.getByToken(IsoTracksToken_, IsoTracks);
          
    std::cout << "\n====================================== We have  " << AllMuCands->size() << " L3 + Trk muons " << std::endl;
    std::cout << "                                       of which "   << L3MuCands->size()  << " L3       muons " << std::endl;
    // edm::LogDebug("HLTTriMuonIsolation") << "\n======================================" << std::endl;
    LogDebug("HLTTriMuonIsolation") << "======================================" << std::endl;
    
    std::cout << "L3 muons" << std::endl;
    for (reco::RecoChargedCandidateCollection::const_iterator imu = L3MuCands->begin(); imu != L3MuCands->end(); ++imu){
        std::cout << "\t pt " << imu->pt() << "\teta " << imu->eta() << "\tphi " << imu->phi() << std::endl;
    }

    std::cout << "L3 + Trk after merging muons" << std::endl;
    for (reco::RecoChargedCandidateCollection::const_iterator jmu = AllMuCands->begin(); jmu != AllMuCands->end(); ++jmu){
        std::cout << "\t pt " << jmu->pt() << "\teta " << jmu->eta() << "\tphi " << jmu->phi() << std::endl;
    }

    if (AllMuCands->size() >= 3){

        // Create the 3-muon candidates
        // loop over L3/Trk muons and create all combinations
        for (unsigned int i = 0; i != AllMuCands->size(); ++i){
            const reco::TrackRef &tk_i = (*AllMuCands)[i].track();
            for (unsigned int j = i+1; j != AllMuCands->size(); ++j){
                const reco::TrackRef &tk_j = (*AllMuCands)[j].track();
                for (unsigned int k = j+1; k != AllMuCands->size(); ++k){
                    const reco::TrackRef &tk_k = (*AllMuCands)[k].track();

                    std::cout << "\n" << __LINE__ << "]\t" << "using muon indices i " << i << " j " << j << " k " << k << std::endl;

                    // consider triples only when at least two muons pass the previous filter 
                    // Match by dR and not by TrackRef because muons mix L3 & Trk reconstructions and tracks are reshuffled...

/*
                    int passingPreviousFilter = 0;
                    for (std::vector<reco::RecoChargedCandidateRef>::const_iterator imu = PassedL3Muons.begin(); imu != PassedL3Muons.end(); ++imu){
                        reco::TrackRef candTrkRef = (*imu)->get<reco::TrackRef>();
                        // if (tk_i == candTrkRef) passingPreviousFilter++;
                        // if (tk_j == candTrkRef) passingPreviousFilter++; 
                        // if (tk_k == candTrkRef) passingPreviousFilter++;
                        if (reco::deltaR2(tk_i->momentum(), candTrkRef->momentum()) < (0.05*0.05)) passingPreviousFilter++;
                        if (reco::deltaR2(tk_j->momentum(), candTrkRef->momentum()) < (0.05*0.05)) passingPreviousFilter++; 
                        if (reco::deltaR2(tk_k->momentum(), candTrkRef->momentum()) < (0.05*0.05)) passingPreviousFilter++;
                    }                    
                    std::cout << __LINE__ << "]\t" << "how many muons pass the previous filter?\t" << passingPreviousFilter << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "how many muons pass the previous filter?\t" << passingPreviousFilter << std::endl;
                    if (passingPreviousFilter < 2) continue;
//*/

                    int passingPreviousFilter_1 = 0;
                    int passingPreviousFilter_2 = 0;
                    int passingPreviousFilter_3 = 0;
    
                    for (std::vector<reco::RecoChargedCandidateRef>::const_iterator imu = PassedL3Muons.begin(); imu != PassedL3Muons.end(); ++imu){
                        reco::TrackRef candTrkRef = (*imu)->get<reco::TrackRef>();
                        if (reco::deltaR2(tk_i->momentum(), candTrkRef->momentum()) < (0.03*0.03)) passingPreviousFilter_1++;
                        if (reco::deltaR2(tk_j->momentum(), candTrkRef->momentum()) < (0.03*0.03)) passingPreviousFilter_2++; 
                        if (reco::deltaR2(tk_k->momentum(), candTrkRef->momentum()) < (0.03*0.03)) passingPreviousFilter_3++;
                    }                    
                    std::cout << __LINE__ << "]\t" << "how many muon_1 pass the previous filter?\t" << passingPreviousFilter_1 << std::endl;
                    std::cout << __LINE__ << "]\t" << "how many muon_2 pass the previous filter?\t" << passingPreviousFilter_2 << std::endl;
                    std::cout << __LINE__ << "]\t" << "how many muon_3 pass the previous filter?\t" << passingPreviousFilter_3 << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "how many muons pass the previous filter?\t" << passingPreviousFilter << std::endl;
                    if (((passingPreviousFilter_1 < 1) & (passingPreviousFilter_2 < 1)) ||
                        ((passingPreviousFilter_1 < 1) & (passingPreviousFilter_3 < 1)) ||
                        ((passingPreviousFilter_2 < 1) & (passingPreviousFilter_3 < 1))) continue;

                    // Create a composite candidate to be a tau
                    reco::CompositeCandidate Tau;

                    // sort the muons by pt and add them to the tau
                    reco::RecoChargedCandidateCollection Daughters;
                    
                    Daughters.push_back((*AllMuCands)[i]);
                    Daughters.push_back((*AllMuCands)[j]);
                    Daughters.push_back((*AllMuCands)[k]);

                    // start building the tau
                    int                      charge   = Daughters[0].charge() + Daughters[1].charge() + Daughters[2].charge();
                    math::XYZTLorentzVectorD taup4    = Daughters[0].p4()     + Daughters[1].p4()     + Daughters[2].p4()    ;
                    int                      tauPdgId = charge > 0? 15 : -15;
                    
                    for (reco::RecoChargedCandidateCollection::const_iterator idau = Daughters.begin(); idau != Daughters.end(); ++idau){
                        std::cout <<  __LINE__ << "]\t" << "\tdaughter pt\t" << idau->pt() << std::endl; 
                    }
                                        
                    std::sort(Daughters.begin(), Daughters.end(), ptComparer);

                    std::cout << __LINE__ << "]\t" << "...sorting..." << std::endl;

                    for (reco::RecoChargedCandidateCollection::const_iterator jdau = Daughters.begin(); jdau != Daughters.end(); ++jdau){
                        std::cout <<  __LINE__ << "]\t" << "\tdaughter pt\t" << jdau->pt() << std::endl; 
                    }

                    Tau.addDaughter((Daughters)[0], "Muon_1");
                    Tau.addDaughter((Daughters)[1], "Muon_2");
                    Tau.addDaughter((Daughters)[2], "Muon_3");

                    Tau.setP4(taup4);
                    Tau.setCharge(charge);
                    Tau.setPdgId(tauPdgId);
                    Tau.setVertex((Daughters)[0].vertex()); // assign the leading muon vertex as tau vertex FIXME!
                    
                    std::cout << __LINE__ << "]\t" << "filled tau" << std::endl;
                    std::cout << __LINE__ << "]\t" << "tau pt\t" << Tau.pt() << "\ttau eta\t" << Tau.eta() << "\ttau phi\t" << Tau.phi() << "\ttau mass\t" << Tau.mass() << "\ttau charge\t" << Tau.charge() << std::endl;
                    std::cout << __LINE__ << "]\t" << "tau vtx x\t" << Tau.vx() << "\ttau vtx y\t" << Tau.vy() << "\ttau vtx z\t" << Tau.vz()  << std::endl;
                    std::cout << __LINE__ << "]\t" << "\tdaugther(0) pt\t" << Tau.daughter(0)->pt() << "\teta\t" << Tau.daughter(0)->eta() << "\tphi\t" << Tau.daughter(0)->phi() << "\tcharge\t" << Tau.daughter(0)->charge() << "\tvtx x\t" << Tau.daughter(0)->vx() << "\tvtx y\t" << Tau.daughter(0)->vy() << "\tvtx z\t" << Tau.daughter(0)->vz() << std::endl;
                    std::cout << __LINE__ << "]\t" << "\tdaugther(1) pt\t" << Tau.daughter(1)->pt() << "\teta\t" << Tau.daughter(1)->eta() << "\tphi\t" << Tau.daughter(1)->phi() << "\tcharge\t" << Tau.daughter(1)->charge() << "\tvtx x\t" << Tau.daughter(1)->vx() << "\tvtx y\t" << Tau.daughter(1)->vy() << "\tvtx z\t" << Tau.daughter(1)->vz() << std::endl;
                    std::cout << __LINE__ << "]\t" << "\tdaugther(2) pt\t" << Tau.daughter(2)->pt() << "\teta\t" << Tau.daughter(2)->eta() << "\tphi\t" << Tau.daughter(2)->phi() << "\tcharge\t" << Tau.daughter(2)->charge() << "\tvtx x\t" << Tau.daughter(2)->vx() << "\tvtx y\t" << Tau.daughter(2)->vy() << "\tvtx z\t" << Tau.daughter(2)->vz() << std::endl;
                    
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "filled tau" << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "tau pt\t" << Tau.pt() << "\ttau eta\t" << Tau.eta() << "\ttau phi\t" << Tau.phi() << "\ttau mass\t" << Tau.mass() << "\ttau charge\t" << Tau.charge() << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "\tdaugther(0) pt\t" << Tau.daughter(0)->pt() << "\teta\t" << Tau.daughter(0)->eta() << "\tphi\t" << Tau.daughter(0)->phi() << "\tcharge\t" << Tau.daughter(0)->charge() << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "\tdaugther(1) pt\t" << Tau.daughter(1)->pt() << "\teta\t" << Tau.daughter(1)->eta() << "\tphi\t" << Tau.daughter(1)->phi() << "\tcharge\t" << Tau.daughter(1)->charge() << std::endl;
                    // edm::LogDebug("HLTTriMuonIsolation") << "]\t" << "\tdaugther(2) pt\t" << Tau.daughter(2)->pt() << "\teta\t" << Tau.daughter(2)->eta() << "\tphi\t" << Tau.daughter(2)->phi() << "\tcharge\t" << Tau.daughter(2)->charge() << std::endl;                    
                                        
                    Taus->push_back(Tau);
                }
            }
        }

        // Loop over taus and further select
        for (reco::CompositeCandidateCollection::const_iterator itau = Taus->begin(); itau != Taus->end(); ++itau){
            if (    itau->pt()   < TriMuonPtCut_  ) continue;
            if (    itau->mass() < MinTriMuonMass_) continue;
            if (    itau->mass() > MaxTriMuonMass_) continue;
            if (abs(itau->eta()) > TriMuonEtaCut_ ) continue;
            if (itau->daughter(0)->pt() < Muon1PtCut_) continue;
            if (itau->daughter(1)->pt() < Muon2PtCut_) continue;
            if (itau->daughter(2)->pt() < Muon3PtCut_) continue;
            if ((abs(itau->charge()) != TriMuonAbsCharge_) & (TriMuonAbsCharge_ >= 0)) continue;
            if (abs(itau->daughter(0)->vz() - itau->vz()) > MaxDZ_) continue;
            if (abs(itau->daughter(1)->vz() - itau->vz()) > MaxDZ_) continue;
            if (abs(itau->daughter(2)->vz() - itau->vz()) > MaxDZ_) continue;
            
            double sumPt = -itau->pt(); // remove the candidate pt from the iso sum

            std::cout << __LINE__ << "]\t" << "\tinitial sum pt\t" << sumPt << std::endl;
            
            // FIXME add track quality requirements
            for (reco::TrackCollection::const_iterator itrk = IsoTracks->begin(); itrk != IsoTracks->end(); ++itrk){
                if (reco::deltaR2(itrk->momentum(), itau->p4()) > IsoConeSize_*IsoConeSize_) continue;
                if (abs(itrk->vz() - itau->vz()) > MaxDZ_) continue;
                sumPt += itrk->pt();
                std::cout << __LINE__ << "]\t\t" << "\tcounting sum pt\t" << sumPt << std::endl;
            }
            std::cout << __LINE__ << "]\t" << "\tfinal sum pt\t" << sumPt << std::endl;
            
            if (std::max(0., sumPt) > (EnableAbsIso_ * ChargedAbsIsoCut_) || 
                std::max(0., sumPt) > (EnableRelIso_ * ChargedRelIsoCut_ * itau->pt())) continue;
            
            SelectedTaus->push_back(*itau); 
        }
    }
            
    // Put the vector of 3-muon candidates in the event 
    iEvent.put(std::move(Taus)        , "Taus"        );
    iEvent.put(std::move(SelectedTaus), "SelectedTaus");
}

void
HLTTriMuonIsolation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("L3MuonsSrc"        , edm::InputTag("hltIterL3FromL2MuonCandidates"  ));
  desc.add<edm::InputTag>("AllMuonsSrc"       , edm::InputTag("hltGlbTrkMuonCands"             ));
  desc.add<edm::InputTag>("L3DiMuonsFilterSrc", edm::InputTag("hltDiMuonForTau3MuDzFiltered0p3"));
  desc.add<edm::InputTag>("IsoTracksSrc"      , edm::InputTag("hltIter2L3FromL2MuonMerged"     ));
  desc.add<double>("Muon1PtCut"      , 5.   );
  desc.add<double>("Muon2PtCut"      , 3.   );
  desc.add<double>("Muon3PtCut"      , 0.   );
  desc.add<double>("TriMuonPtCut"    , 8.   );
  desc.add<double>("TriMuonEtaCut"   , 2.5  );
  desc.add<double>("ChargedAbsIsoCut", 3.0  );
  desc.add<double>("ChargedRelIsoCut", 0.1  );
  desc.add<double>("IsoConeSize"     , 0.5  );
  desc.add<double>("MinTriMuonMass"  , 0.5  );
  desc.add<double>("MaxTriMuonMass"  , 2.8  );
  desc.add<int>   ("TriMuonAbsCharge", -1   );
  desc.add<double>("MaxDZ"           , 0.3  );
  desc.add<bool>  ("EnableRelIso"    , false);
  desc.add<bool>  ("EnableAbsIso"    , true );
  descriptions.add("hltTriMuonIsolationProducer",desc);
}



#endif
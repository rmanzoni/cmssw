#ifndef METAnalyzer_H
#define METAnalyzer_H

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class METAnalyzer : public JME::Analyzer {
    public:
        // construction/destruction
        explicit METAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~METAnalyzer();
        typedef std::vector<edm::InputTag> vInputTag;
    private:

        // member functions
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
        pat::MET getUncorrectedRecoil(const pat::MET& input);

    private:

        double dRgenMatching_;
        bool   isMC_;
        TVector2 RecoilVec;
        TVector2 BosonVec;
        TVector2 GenBosonVec;
        ROOT::Math::SMatrix<double,2> rotationMatrix;

        edm::InputTag srcJet_;
        edm::InputTag srcJetPF_;

        edm::InputTag srcZboson_;
        edm::InputTag srcLeptons_;
        edm::InputTag srcVertex_;

        edm::InputTag srcGenMet_;

        edm::InputTag srcGenJets_;
        edm::InputTag srcGenJetsCleaned_;

        edm::InputTag srcGenParticles_;
        edm::InputTag srcGenEventInfo_;


        edm::InputTag srcMetFiltersBits_;
        edm::InputTag srcTriggerBits_;
        edm::InputTag srcTriggerPrescales_;

        std::vector<edm::EDGetTokenT<pat::METCollection > >  srcRecoils_;
        vInputTag srcRecoilTags_;
        std::vector<edm::EDGetTokenT<pat::METCollection > >  srcMETs_;
        vInputTag srcMETTags_;

        edm::EDGetTokenT<std::vector<pat::Jet>> srcJetToken_;
        edm::EDGetTokenT<std::vector<pat::Jet>> srcJetPFToken_;
        edm::EDGetTokenT<std::vector<reco::Particle>>   srcZbosonToken_;
        edm::EDGetTokenT<reco::CandidateView>           srcLeptonsToken_;
        edm::EDGetTokenT<reco::VertexCollection>        srcVertexToken_;

        edm::EDGetTokenT<std::vector<pat::MET>>         srcGenMetToken_;

        edm::EDGetTokenT<std::vector<reco::GenJet>>     srcGenJetsToken_;
        edm::EDGetTokenT<pat::JetCollection>            srcGenJetsCleanedToken_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> srcGenParticlesToken_;
        edm::EDGetTokenT<GenEventInfoProduct> srcGenEventInfoToken_;

        edm::EDGetTokenT<edm::TriggerResults>  srcMetFiltersBitsToken_;
        edm::EDGetTokenT<edm::TriggerResults>  srcTriggerBitsToken_;
        edm::EDGetTokenT<pat::PackedTriggerPrescales>  srcTriggerPrescalesToken_;

        // MC weight factor
        float& eventMCWeight = tree["eventMCWeight"].write<float>(); 
        // Generator level boson 
        float& GenBoson_Pt_  = tree["GenBoson_Pt"].write<float>();
        float& GenBoson_Eta_ = tree["GenBoson_Eta"].write<float>();
        float& GenBoson_Phi_ = tree["GenBoson_Phi"].write<float>();
        float& GenBoson_M_   = tree["GenBoson_M"].write<float>();
        int&   GenBoson_daughter_ = tree["GenBoson_daughter"].write<int>();

        // Generator level jets
        int& NGenJets_        = tree["NGenJets"].write<int>();
        int& NGenJetsCleaned_ = tree["NGenJetsCleaned"].write<int>();
        int& NGenMatchedJets_ = tree["NGenMatchedJets"].write< int>();
        int& NGenMatchedJetsPF_ = tree["NGenMatchedJetsPF"].write< int>();

        float& GenLeadingJet_Pt_  = tree["GenJet1_Pt"].write<float>();
        float& GenLeadingJet_Eta_ = tree["GenJet1_Eta"].write<float>();
        float& GenLeadingJet_Phi_ = tree["GenJet1_Phi"].write<float>();
        float& GenLeadingJet_M_   = tree["GenJet1_M"].write<float>();

        float& GenLeadingJetCleaned_Pt_  = tree["GenJet1Cleaned_Pt"].write<float>();
        float& GenLeadingJetCleaned_Eta_ = tree["GenJet1Cleaned_Eta"].write<float>();
        float& GenLeadingJetCleaned_Phi_ = tree["GenJet1Cleaned_Phi"].write<float>();
        float& GenLeadingJetCleaned_M_   = tree["GenJet1Cleaned_M"].write<float>();

        float& GenTrailingJet_Pt_  = tree["GenJet2_Pt"].write<float>();
        float& GenTrailingJet_Eta_ = tree["GenJet2_Eta"].write<float>();
        float& GenTrailingJet_Phi_ = tree["GenJet2_Phi"].write<float>();
        float& GenTrailingJet_M_   = tree["GenJet2_M"].write<float>();

        float& GenTrailingJetCleaned_Pt_  = tree["GenJet2Cleaned_Pt"].write<float>();
        float& GenTrailingJetCleaned_Eta_ = tree["GenJet2Cleaned_Eta"].write<float>();
        float& GenTrailingJetCleaned_Phi_ = tree["GenJet2Cleaned_Phi"].write<float>();
        float& GenTrailingJetCleaned_M_   = tree["GenJet2Cleaned_M"].write<float>();

	// Generator level recoil

        float& GenRecoil_sumEt_ = tree["GenRecoil_sumEt"].write<float>();
        float& GenRecoil_Pt_    = tree["GenRecoil_Pt"].write<float>();
        float& GenRecoil_Phi_   = tree["GenRecoil_Phi"].write<float>();
        float& GenRecoil_PerpZ_ = tree["GenRecoil_PerpZ"].write<float>();
        float& GenRecoil_LongZ_ = tree["GenRecoil_LongZ"].write<float>();

        float& GenBoson_PerpU_  = tree["GenBoson_PerpZ"].write<float>();
        float& GenBoson_LongU_  = tree["GenBoson_LongZ"].write<float>();

        // recoils
        std::vector<std::reference_wrapper<float>> recoilReferences_; // assign via function
        std::vector<std::string> recoilAttributes_ = {"sumEt", "sumEtFraction", "Pt", "Phi", "PerpZ", "LongZ", "Boson_PerpU", "Boson_LongU", "Cov11", "Cov00", "Cov01", "Cov10"}; // to be extended?
    // METs
        std::vector<std::reference_wrapper<float>> METReferences_; // assign via function
        std::vector<std::string> METAttributes_ = {"sumEt", "Pt", "Phi", "Cov11", "Cov00", "Cov01", "Cov10"}; // to be extended?


        // jets and boson
        int& NCleanedJets_ =  tree["NCleanedJets"].write<int>();
        int& NCleanedJetsPF_ =  tree["NCleanedJetsPF"].write<int>();
        int& NVertex_      =  tree["NVertex"].write<int>();

        float& LeadingJet_Pt_  = tree["Jet1_Pt"].write<float>();
        float& LeadingJet_Eta_ = tree["Jet1_Eta"].write<float>();
        float& LeadingJet_Phi_ = tree["Jet1_Phi"].write<float>();
        float& LeadingJet_M_   = tree["Jet1_M"].write<float>();

        float& TrailingJet_Pt_  = tree["Jet2_Pt"].write<float>();
        float& TrailingJet_Eta_ = tree["Jet2_Eta"].write<float>();
        float& TrailingJet_Phi_ = tree["Jet2_Phi"].write<float>();
        float& TrailingJet_M_   = tree["Jet2_M"].write<float>();

        float& LeadingJetPF_Pt_  = tree["Jet1PF_Pt"].write<float>();
        float& LeadingJetPF_Eta_ = tree["Jet1PF_Eta"].write<float>();
        float& LeadingJetPF_Phi_ = tree["Jet1PF_Phi"].write<float>();
        float& LeadingJetPF_M_   = tree["Jet1PF_M"].write<float>();

        float& TrailingJetPF_Pt_  = tree["Jet2PF_Pt"].write<float>();
        float& TrailingJetPF_Eta_ = tree["Jet2PF_Eta"].write<float>();
        float& TrailingJetPF_Phi_ = tree["Jet2PF_Phi"].write<float>();
        float& TrailingJetPF_M_   = tree["Jet2PF_M"].write<float>();

        float& Boson_Pt_     =  tree["Boson_Pt"].write<float>();
        float& Boson_Phi_    =  tree["Boson_Phi"].write<float>();
        float& Boson_Eta_    =  tree["Boson_Eta"].write<float>();
        float& Boson_M_      =  tree["Boson_M"].write<float>();
        int& Boson_daughter_ =  tree["Boson_daughter"].write<int>();

        float& LeadingLepton_Pt_     =  tree["LeadingLepton_Pt"].write<float>();
        float& LeadingLepton_Phi_    =  tree["LeadingLepton_Phi"].write<float>();
        float& LeadingLepton_Eta_    =  tree["LeadingLepton_Eta"].write<float>();
        float& LeadingLepton_M_      =  tree["LeadingLepton_M"].write<float>();

        float& TrailingLepton_Pt_     =  tree["TrailingLepton_Pt"].write<float>();
        float& TrailingLepton_Phi_    =  tree["TrailingLepton_Phi"].write<float>();
        float& TrailingLepton_Eta_    =  tree["TrailingLepton_Eta"].write<float>();
        float& TrailingLepton_M_      =  tree["TrailingLepton_M"].write<float>();

        std::vector<float>& AllJets_Pt_  = tree["AllJets_Pt"].write<std::vector<float>>();
        std::vector<float>& AllJets_Eta_ = tree["AllJets_Eta"].write<std::vector<float>>();
        std::vector<float>& AllJets_Phi_ = tree["AllJets_Phi"].write<std::vector<float>>();
        std::vector<float>& AllJets_M_   = tree["AllJets_M"].write<std::vector<float>>();

        std::vector<float>& AllJetsPF_Pt_  = tree["AllJetsPF_Pt"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_Eta_ = tree["AllJetsPF_Eta"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_Phi_ = tree["AllJetsPF_Phi"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_M_   = tree["AllJetsPF_M"].write<std::vector<float>>();

        std::vector<float>& GenMatchedJets_Pt_  = tree["GenMatchedJets_Pt"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Eta_ = tree["GenMatchedJets_Eta"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Phi_ = tree["GenMatchedJets_Phi"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_M_   = tree["GenMatchedJets_M"].write<std::vector<float>>();

        std::vector<float>& GenMatchedJetsPF_Pt_  = tree["GenMatchedJets_Pt"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_Eta_ = tree["GenMatchedJets_Eta"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_Phi_ = tree["GenMatchedJets_Phi"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_M_   = tree["GenMatchedJets_M"].write<std::vector<float>>();

        /// met filter
        int& flag_HBHENoiseFilter_      =  tree["flag_HBHENoiseFilter"].write<int>();
        int& flag_CSCTightHaloFilter_   =  tree["flag_CSCTightHaloFilter"].write<int>();
        int& flag_hcalLaserEventFilter_ =  tree["flag_hcalLaserEventFilter"].write<int>();
        int& flag_EcalDeadCellTriggerPrimitiveFilter_ = tree["flag_EcalDeadCellTriggerPrimitiveFilter"].write<int>();
        int& flag_goodVertices_          = tree["flag_goodVertices"].write<int>();
        int& flag_trackingFailureFilter_ = tree["flag_trackingFailureFilter"].write<int>();
        int& flag_eeBadScFilter_         = tree["flag_eeBadScFilter"].write<int>();
        int& flag_ecalLaserCorrFilter_   = tree["flag_ecalLaserCorrFilter"].write<int>();
        int& flag_trkPOGFilters_         =  tree["flag_trkPOGFilters"].write<int>();  
        int& flag_trkPOG_manystripclus53X_    = tree["flag_trkPOG_manystripclus53X"].write<int>();
        int& flag_trkPOG_toomanystripclus53X_ = tree["flag_trkPOG_toomanystripclus53X"].write<int>();
        int& flag_trkPOG_logErrorTooManyClusters_ = tree["flag_trkPOG_logErrorTooManyClusters"].write<int>();
        int& flag_METFilters_ = tree["flag_METFilters"].write<int>();


        std::vector<std::string>& DoubleMuPaths_    = tree["DoubleMuPaths"].write<std::vector<std::string> >();
        std::vector<std::string>& DoubleElePaths_   = tree["DoubleElePaths"].write<std::vector<std::string> >();
        std::vector<std::string>& DoubleTauPaths_   = tree["DoubleTauPaths"].write<std::vector<std::string> >();

};

#endif

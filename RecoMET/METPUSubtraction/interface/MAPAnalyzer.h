#ifndef MAPAnalyzer_H
#define MAPAnalyzer_H


#include "RecoMET/METPUSubtraction/interface/Analyzer.h"
class MAPAnalyzer : public JME::Analyzer {
    public:
        // construction/destruction
        explicit MAPAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~MAPAnalyzer();
    private:
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
    private:
        edm::InputTag srcVariableNames_;
        edm::InputTag srcVariables_;
        std::vector<std::string> variableNamesToSave_;
        std::map<std::string, std::reference_wrapper<Float_t>> var_;
        edm::EDGetTokenT<std::vector<Float_t>> srcVariablesToken_;
        edm::EDGetTokenT<std::vector<std::string>> srcVariableNamesToken_;
};

#endif

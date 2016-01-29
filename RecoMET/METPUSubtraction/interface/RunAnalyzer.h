#pragma once

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

class RunAnalyzer: public JME::Analyzer {
    public:
        explicit RunAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~RunAnalyzer();

        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void analyze(edm::Event const&, edm::EventSetup const&) override {};

    private:
        edm::EDGetTokenT<GenRunInfoProduct> genRunInfoToken_;

    private:
        ULong64_t& run_ = tree["run"].write<ULong64_t>();
        float& xsec_ = tree["cross_section"].write<float>(false);
};

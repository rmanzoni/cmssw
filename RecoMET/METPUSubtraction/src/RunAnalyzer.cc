#include "FWCore/Utilities/interface/BranchType.h"

#include "JMEAnalysis/JMEValidator/interface/RunAnalyzer.h"

RunAnalyzer::RunAnalyzer(const edm::ParameterSet& iConfig): JME::Analyzer(iConfig),
    genRunInfoToken_ (consumes<GenRunInfoProduct, edm::BranchType::InRun>(edm::InputTag("generator")))
{
    // Empty
}

RunAnalyzer::~RunAnalyzer() {
    // Empty
}

void RunAnalyzer::endRun(const edm::Run& run, const edm::EventSetup& iSetup) {
    run_ = run.id().run();

    edm::Handle<GenRunInfoProduct> genInfo;
    if (run.getByToken(genRunInfoToken_, genInfo)) {
        xsec_ = genInfo->crossSection();
    }

    tree.fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RunAnalyzer);

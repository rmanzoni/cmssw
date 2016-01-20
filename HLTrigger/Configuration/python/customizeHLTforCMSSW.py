import FWCore.ParameterSet.Config as cms

#
# reusable functions

def producers_by_type(process, *types):
    return (module for module in process._Process__producers.values() if module._TypedParameterizable__type in types)

def esproducers_by_type(process, *types):
    return (module for module in process._Process__esproducers.values() if module._TypedParameterizable__type in types)

#
# one action function per PR - put the PR number into the name of the function

def customiseFor1nnnn(process):
    return process

#
# CMSSW version specific customizations

def customiseHLTforCMSSW(process, menuType="GRun", fastSim=False):
    import os
    cmsswVersion = os.environ['CMSSW_VERSION']

    if cmsswVersion >= "CMSSW_7_6":
        process = customiseFor1nnnn(process)

    return process

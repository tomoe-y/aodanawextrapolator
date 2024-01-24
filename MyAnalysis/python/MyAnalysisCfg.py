# Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def MyxAODAnalysisCfg(ConfigFlags, name = "MyxAODAnalysis", **kwargs):
    acc= ComponentAccumulator()

    from TrigDecisionTool.TrigDecisionToolConfig import TrigDecisionToolCfg
    trig_dec_tool = acc.getPrimaryAndMerge(TrigDecisionToolCfg(ConfigFlags))
    kwargs.setdefault("trigDecisionTool", trig_dec_tool)

    from TrkConfig.AtlasExtrapolatorConfig import MuonExtrapolatorCfg
    kwargs.setdefault("Extrapolator", acc.popToolsAndMerge(MuonExtrapolatorCfg(ConfigFlags)))

    from MuonSelectorTools.MuonSelectorToolsConfig import MuonSelectionToolCfg
    sel_tool = acc.getPrimaryAndMerge(MuonSelectionToolCfg(ConfigFlags, MaxEta=2.5, MuQuality=1, TurnOffMomCorr=True))
    kwargs.setdefault("MuonSelectionTool", sel_tool)

    #histSvc.Output += ["%s DATAFILE='%s' OPT='RECREATE'" % ("analysis", flags.Output.HISTFileName)]
    #histSvc = CompFactory.THistSvc(Output=["ANALYSIS DATAFILE='analysis.root', OPT='RECREATE'"])
    histSvc = CompFactory.THistSvc(Output=["%s DATAFILE='%s' OPT='RECREATE'" % ("ANALYSIS", ConfigFlags.Output.HISTFileName)])
    acc.addService(histSvc)

    the_alg = CompFactory.MyxAODAnalysis(name = name,**kwargs)
    acc.addEventAlgo(the_alg, primary = True)
    return acc

# Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration

###################################################
# job options file to run the MyAnalysis based on component accumulator     #
###################################################

def SetupArgParser():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--inputFile", "-i", default=["/gpfs/fs8001/youhei/L2MuonSA/dataset_aod_official/valid1.801164.P8B_A14_CTEQ6L1_bb_Jpsi1S_mu6mu4.recon.AOD.e8514_e8528_s4111_s4114_r14781_tid33979104_00/AOD.33979104._000001.pool.root.1"], 
                        help="Input file to run on ", nargs="+")
    parser.add_argument("--maxEvents", default=-1, type=int, help="How many events shall be run maximally")
    parser.add_argument("--outputFile", default=["MyAnalysis.root"], help="Output file to run on")
    return parser

def setupServicesCfg(flags):
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    result = MainServicesCfg(flags)
    ### Setup the file reading
    from AthenaConfiguration.Enums import Format
    if flags.Input.Format == Format.POOL:
        from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
        result.merge(PoolReadCfg(flags))
    elif flags.Input.Format == Format.BS:
        from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg
        result.merge(ByteStreamReadCfg(flags))

    from MuonConfig.MuonGeometryConfig import MuonGeoModelCfg
    result.merge(MuonGeoModelCfg(flags))
    # from MuonConfig.MuonGeometryConfig import MuonIdHelperSvcCfg
    # result.merge(MuonIdHelperSvcCfg(flags))

    return result

if __name__ == "__main__":
    from AthenaConfiguration.AllConfigFlags import initConfigFlags
    args = SetupArgParser().parse_args()

    flags = initConfigFlags()
    flags.Input.Files = args.inputFile
    flags.Exec.MaxEvents = args.maxEvents
    flags.Output.HISTFileName = args.outputFile
    flags.lock()

    cfg = setupServicesCfg(flags)
    msgService = cfg.getService('MessageSvc')
    msgService.Format = "S:%s E:%e % F%128W%S%7W%R%T  %0W%M"

    flags.dump()

    #from AthenaConfiguration.ComponentFactory import CompFactory
    #histSvc = CompFactory.THistSvc()
    #histSvc.MaxFileSize = -1   # speeds up jobs that output lots of histograms
    #histSvc.Output += ["%s DATAFILE='%s' OPT='RECREATE'" % ("analysis", flags.Output.HISTFileName)]
    #cfg.addService(histSvc)

    from MyAnalysis.MyAnalysisCfg import MyxAODAnalysisCfg
    cfg.merge(MyxAODAnalysisCfg(flags))
    cfg.printConfig(withDetails=True, summariseProps=True)
    flags.dump()


    sc = cfg.run(flags.Exec.MaxEvents)
    if not sc.isSuccess():
        import sys
        sys.exit("Execution failed")

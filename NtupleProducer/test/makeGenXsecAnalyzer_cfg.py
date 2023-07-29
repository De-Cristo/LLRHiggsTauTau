from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
options = VarParsing ('python')
options.parseArguments()

process = cms.Process('makeGenXsecAnalyzer')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 25000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",    
    fileNames  = cms.untracked.vstring(options.inputFiles))

process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")
process.p = cms.Path(process.xsec)
process.s = cms.Schedule(process.p);

### CMSSW command line parameter parser                                                                                                                                                               
from FWCore.ParameterSet.VarParsing import VarParsing
import os
options = VarParsing ('python')

## data or MC options                                                                                                                                                                                 
options.register (
    'outputName',"gentree.root",VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'name for the outputfile');

options.register (
    'inputFileList',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'list with files');

options.register(
    'xsec',-1.,VarParsing.multiplicity.singleton, VarParsing.varType.float,
    'external value for sample cross section, in case of data it is fixed to 0.001');

options.register(
    'isMiniAOD',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'flag to tell if one is running on miniAOD or GEN files');

options.parseArguments()

# Define the CMSSW process
import FWCore.ParameterSet.Config as cms
process = cms.Process("GENTREE")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

# Define the input source
if options.inputFileList:
    with open(options.inputFileList,"r") as f:
        for lines in f:
            options.inputFiles.append(lines)
process.source = cms.Source("PoolSource", 
        fileNames = cms.untracked.vstring(options.inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputName)
)

# Processing setup
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents)
)

# Make the tree 
process.gentree = cms.EDAnalyzer("GenTreeMaker",
        xsec       = cms.double(options.xsec),   
        isMiniAOD  = cms.bool(options.isMiniAOD),
        genMet     = cms.InputTag("slimmedMETs"),
        genJets    = cms.InputTag("slimmedGenJets"),
        genParticles  = cms.InputTag("prunedGenParticles"),
        lheEvent = cms.InputTag("externalLHEProducer"),
        genEvent = cms.InputTag("generator"),
)

if not options.isMiniAOD:
    process.gentree.genParticles   = cms.InputTag("genParticles");
    process.gentree.genJets    = cms.InputTag("ak4GenJetsNoNu");
    process.gentree.genMet     = cms.InputTag("genMetTrue");
    process.gentree.lheEvent   = cms.InputTag("source");
    process.gentree.lheRunInfo = cms.InputTag("source");
    
process.gentreePath = cms.Path(process.gentree)


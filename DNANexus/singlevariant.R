#!/usr/bin/env Rscript

args=(commandArgs(TRUE))
num=as.numeric(args[1])# 1
gdsfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])    

#num=22
#stats=as.character(args[3])# 1
stats="Score"
#####
##### single variant association tests
#library(SeqArray)
#library(SeqVarTools)
#library(data.table)
#library(GENESIS)
#library(GWASTools)

cat('\nReading in packages for analysis...\n')
.libPaths(c("rpackages4_1_3",.libPaths()))

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("TOPMed_AFib_pipeline/DNANexus/singlevariant_modified.R")
print(num)

#perform test
singlevariant(num=num,gdsfile=gdsfile,varfile=NULL,varidfile=NULL,phenfile=phenfile,nullfile=nullfile,stat=stats,outfile=outfile,mcount=10)

sessionInfo()
quit("no")

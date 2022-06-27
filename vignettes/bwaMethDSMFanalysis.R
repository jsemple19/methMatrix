## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(magrittr)
library(methMatrix)

## -----------------------------------------------------------------------------
#genomeFile="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"
#bamFile="/Users/semple/Documents/MeisterLab/sequencingData/20190218_dSMFv016v020_N2gw_2x150pe/properpair/aln/dS16N2gw_20190218.noOL.bam"
bamFile<-system.file("extdata", "aln/dS03-N2_20181119.noOL.bam",
                       package="methMatrix",mustWork=TRUE)
genomeFile<-system.file("extdata",
                          "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa",
                          package="methMatrix",mustWork=TRUE)
amplicons<-readRDS(system.file("extdata", "genome/ampliconGR.RDS",
                         package="methMatrix", mustWork=TRUE))
#amplicons=readRDS("/Users/semple/Documents/MeisterLab/dSMF/PromoterPrimerDesign/usefulObjects/ampliconGR.RDS")
names(GenomicRanges::mcols(amplicons))<-"ID"
GenomeInfoDb::seqlevelsStyle(amplicons)<-"ensembl"
regionGR<-amplicons[1]


## ----eval=FALSE---------------------------------------------------------------
#  makeCGorGCbed(genomeFile,"CG")
#  makeCGorGCbed(genomeFile,"GC")

## -----------------------------------------------------------------------------
#bedFileCG="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.CG.bed"
bedFileCG<-system.file("extdata",
                       "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CG.bed",
                       package="methMatrix",mustWork=TRUE)
#bedFileGC="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.GC.bed"
bedFileGC<-system.file("extdata",
                       "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.GC.bed",
                       package="methMatrix",mustWork=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  genomeMotifGR<-nanodsmf::findGenomeMotifs(genomeFile)

## -----------------------------------------------------------------------------
#genomeMotifGR<-readRDS("/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJ#NA13758.WS250.genomic.CGGC_motifs.RDS")
genomeMotifGR<-readRDS(system.file("extdata",
              "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CGGC_motifs.RDS",
                      package="methMatrix", mustWork=TRUE))

## -----------------------------------------------------------------------------
matCG<-getReadMatrix(bamFile=bamFile,genomeFile=genomeFile,bedFile=bedFileCG,region=regionGR, samtoolsPath="~/miniconda3/bin/")
matGC<-getReadMatrix(bamFile,genomeFile,bedFileGC,regionGR,samtoolsPath="~/miniconda3/bin/")

matCG[2:5,18:26]
matGC[2:5,8:16]

## -----------------------------------------------------------------------------
methMat<-combineCGandGCmatrices(matCG,matGC,regionGR,genomeMotifGR)
methMat[2:5,10:18]
colSums(methMat,na.rm=T)


## -----------------------------------------------------------------------------
# read in the table with the list of matrix files
matTableFile<-system.file("extdata", "csv/MatrixLog_ampTSS.csv",
                          package="methMatrix", mustWork=TRUE)
matTable<-read.csv(matTableFile,stringsAsFactors=F)

# modify the matTable path to data
matTable$filename <- system.file("extdata", 
                      gsub("^.*/rds/","rds/", 
                           matTable[,"filename"]), 
                      package="methMatrix", mustWork=TRUE)

# read in one matrix
i=2
dataMatrix<-readRDS(matTable[i,"filename"])

# read in the genomic ranges for TSSs
tssFile<-system.file("extdata","genome/ampliconMaxTSSgr.RDS",
                     package="methMatrix",mustWork=TRUE)
tssWin<-readRDS(tssFile)
names(GenomicRanges::mcols(tssWin))<-"ID"
idx<-match(unique(matTable$region),tssWin$ID)
tssWin<-tssWin[idx]
#saveRDS(tssWin,"./inst/extData/genome/TSSgr.RDS")

regionType<-"TSS"
winSize<-500

matTableFile<-system.file("extdata", 
                          "csv/MatrixLog_relCoord_ampTSS.csv",
                          package="methMatrix", mustWork=TRUE)
matTable<-read.csv(matTableFile,stringsAsFactors=F)
matTable$filename <-system.file("extdata", gsub("^\\./","",matTable[,"filename"]), package="methMatrix", mustWork=TRUE)

allSampleRelCoordMats<-getRelativeCoordMats(matList=matTable,
                                            regionGRs=tssWin,
                                            regionType=regionType, 
                                            anchorCoord=winSize/2)

TSSrelCoord<-convertGRtoRelCoord(tssWin,1,anchorPoint="middle")
tssWinRelCoord<-convertGRtoRelCoord(tssWin,winSize,anchorPoint="middle")

minConversionRate=0.8
maxNAfraction=0.2
seqDate=""
plotAllMatrices(allSampleRelCoordMats, samples, regionGRs=tssWinRelCoord,  
                featureGRs=TSSrelCoord, featureLabel="TSS", 
                regionType=regionType,
                maxNAfraction=maxNAfraction, withAvr=FALSE,
                includeInFileName=seqDate, drawArrow=FALSE)



context("Convert bam files to methylation matrices")


testthat::test_that("getReadMatrix from bam works", {
  bamFile<-"./inst/extData/aln/dS03-N2_20181119.noOL.bam"
  genomeFile<-"./inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa"
  bedFile<-"./inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CG.bed"
  amplicons<-readRDS(system.file("extdata", "genome/ampliconGR.RDS",
                         package="methMatrix", mustWork=TRUE))
  GenomeInfoDb::seqlevels(amplicons)<-gsub("chr","",
                                            GenomeInfoDb::seqlevels(amplicons))
  methMat<-getReadMatrix(bamFile, genomeFile, bedFile, amplicons[1],
                samtoolsPath="/Applications/anaconda3/bin/")
  testthat::expect_equal(dim(methMat),c(742,38))
})


#library(tictoc)
#tic()
testthat::test_that("getSingleMoleculeMatrices  works", {
  genomeFile<-"./inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa"
  amplicons<-readRDS(system.file("extdata", "genome/ampliconGR.RDS",
                                 package="methMatrix", mustWork=TRUE))
  names(GenomicRanges::mcols(amplicons))<-"ID"
  GenomeInfoDb::seqlevels(amplicons)<-gsub("chr","",
                                            GenomeInfoDb::seqlevels(amplicons))
  genomeMotifGR<-readRDS(system.file("extdata",
              "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CGGC_motifs.RDS",
                      package="methMatrix", mustWork=TRUE))
  # read in the table with the list of matrix files
  sampleTableFile<-system.file("extdata", "txt/bwameth_Aligned.txt",
                            package="methMatrix", mustWork=TRUE)
  sampleTable<-read.delim(sampleTableFile, stringsAsFactors=F)
  sampleTable$FileName<-gsub("^\\.","\\./inst/extData",sampleTable$FileName)
  #matTable$filename <-system.file("extdata", matTable[,"filename"],
  #package="methMatrix", mustWork=TRUE)

  matTable<-getSingleMoleculeMatrices(sampleTable, genomeFile,
                                      amplicons[c(1,10:12)],
                            "rawAmp", genomeMotifGR, path="./inst/extData",
                            samtoolsPath="/Applications/anaconda3/bin/",
                            convRatePlots=TRUE, nThreads=2)
  testthat::expect_equal(matTable$fewNAreads,c(502,1075,93,275))
})
#toc()


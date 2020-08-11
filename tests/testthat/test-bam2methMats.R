context("Convert bam files to methylation matrices")


testthat::test_that("getReadMatrix from bam works", {
  #packageDir<-"/Users/semple/Documents/MeisterLab/myPackages/methMatrix"
  #bamFile<-paste0(packageDir,"/inst/extData/aln/dS03-N2_20181119.noOL.bam")
  bamFile<-system.file("extdata", "aln/dS03-N2_20181119.noOL.bam",
                       package="methMatrix",mustWork=TRUE)
  genomeFile<-system.file("extdata",
                          "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa",
                          package="methMatrix",mustWork=TRUE)
  #genomeFile<-paste0(packageDir,
  #                   "/inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa")
  bedFile<-system.file("extdata",
                       "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CG.bed",
                       package="methMatrix",mustWork=TRUE)
  #bedFile<-paste0(packageDir,
  #                "/inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CG.bed")
  amplicons<-readRDS(system.file("extdata", "genome/ampliconGR.RDS",
                         package="methMatrix", mustWork=TRUE))
  GenomeInfoDb::seqlevels(amplicons)<-gsub("chr","",
                                            GenomeInfoDb::seqlevels(amplicons))
  methMat<-getReadMatrix(bamFile, genomeFile, bedFile, amplicons[1],
                samtoolsPath="/Applications/anaconda3/bin/")
  testthat::expect_setequal(dim(methMat),c(742,38))
})



testthat::test_that("getSingleMoleculeMatrices works", {
  #packageDir<-"/Users/semple/Documents/MeisterLab/myPackages/methMatrix"
  #genomeFile<-paste0(packageDir,
  #                   "/inst/extData/genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa")
  genomeFile<-system.file("extdata",
                          "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.fa",
                          package="methMatrix",mustWork=TRUE)
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
  matrixFilePath<-system.file("extdata","", package="methMatrix",mustWork=TRUE)
  sampleTable$FileName<-gsub("^\\.",matrixFilePath,sampleTable$FileName)
  #matTable$filename <-system.file("extdata", matTable[,"filename"],
  #package="methMatrix", mustWork=TRUE)

  matTable<-getSingleMoleculeMatrices(sampleTable, genomeFile,
                                      amplicons[c(1,10:12)],
                            "rawAmp", genomeMotifGR, path=matrixFilePath,
                            samtoolsPath="/Applications/anaconda3/bin/",
                            convRatePlots=TRUE, nThreads=2)
  testthat::expect_setequal(matTable$fewNAreads,c(502,1075,93,275))
})



context("test matrix manipulation")

testthat::test_that("getRelativeCoord orders columns correctly", {
  tss=33450
  dataMatrix<-matrix(c(1,1,0,0,0,1,1,0,0,0,0,0),nrow=3)
  colnames(dataMatrix)<-c(-200,-100,150,230)+tss
  regionGR<-GenomicRanges::GRanges(seqnames=c("chrI"),
                                  ranges=IRanges::IRanges(start=tss-250,
                                                          end=tss+250),
                                  strand="-")
  relMat<-getRelativeCoord(dataMatrix, regionGR,
                           invert=ifelse(GenomicRanges::strand(regionGR)=="+",
                                         F, T))
  testthat::expect_equal(as.numeric(colnames(relMat)),c(20,100,350,450))
})


testthat::test_that("getMatrices subsets the matrix table correctly", {
  #tss=33450
  # read in the table with the list of matrix files
  matTableFile<-system.file("extdata", "csv/MatrixLog_relCoord_ampTSS.csv",
                            package="methMatrix", mustWork=TRUE)
  matTable<-read.csv(matTableFile,stringsAsFactors=F)

  # modify the matTable path to data
  matTable$filename <-system.file("extdata", matTable[,"filename"],
                                  package="methMatrix",
                                  mustWork=TRUE)
  # add matrices from another sample
  matTable_dS03N2<-matTable
  matTable_dS03N2$sample<-"dS03-N2"
  matTable_dS03N2$filename<-gsub("dS02-182","dS03-N2",matTable$filename)

  matTable<-rbind(matTable,matTable_dS03N2)

  # read in the genomic ranges for TSSs
  tssFile<-"./inst/extData/genome/ampliconMaxTSSgr.RDS"
  tssWin<-readRDS(tssFile)
  names(GenomicRanges::mcols(tssWin))<-"ID"
  idx<-match(matTable$region,tssWin$ID)
  tssWin<-tssWin[idx]
  newMatList<-getMatrices(matTable,regionName=c("WBGene00015955","WBGene00009621"))
  testthat::expect_equal(length(newMatList),dim(matTable)[1])
  testthat::expect_equal(dim(newMatList[[3]])[1],matTable$fewNAreads[3])
})



testthat::test_that("mergeSampleMats merges tables split by sample", {
  path="inst/extData"
  regionType="rawAmp"
  samples=c("dS03-N2","dS03-N2")
  allMats<-mergeSampleMats(path, regionType, samples, deleteSplitFiles=F)
  testthat::expect_equal(nrow(allMats),8)
})

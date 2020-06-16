context("Convert bam files to methylation matrices")


testthat::test_that("getting mpileup from bam works", {
  bamfile<-system.file("extdata", "aln/dS03-N2_20181119.noOL.bam",
                       package="methMatrix", mustWork=TRUE)
  genomefile<-system.file("extdata", "c_elegans.PRJNA13758.WS250.genomic.fa",
                          package="methMatrix", mustWork=TRUE)
  bedfile<-system.file("extdata", "c_elegans.PRJNA13758.WS250.genomic.CG.bed",
                       package="methMatrix", mustWork=TRUE)
  testthat::expect_equal()
})

---
title: "bwaMethDSMFanalysis"
author: "Jennifer Semple"
date: "2022-06-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{"Analysing bwaMeth dSMF data"} 
  %\VignetteEngine{knitr::knitr}
  %\usepackage[utf8]{inputenc}
---

# Dual enzyme single molecule footprinting (DSMF) analysis



## Using bwa-meth for bisulfite sequence analysis

Genome wide bisulfite data can be analysed relatively easily using QuasR, which can produce both average methylation frequency and single molecule methylation matrices (reads x positions with 1s/0s indicating methylation status). However QuasR uses Bowtie (v1?) for mapping, which is not the best mapper. Also in bisulfite mode it only considers reads in the FR orientation. Since some libraries seem to map extremely poorly with QuasR (only 10% of reads), we developed a new pipeline with [bwa-meth](https://github.com/brentp/bwa-meth) to align the reads and [MethylDackel](https://github.com/dpryan79/MethylDackel) to call methylation frequency at both CpG and GpC motifs. To get single molecule methylation matrices we used samtools mpileup to call variants. The methMatrix package contains scripts to parse the output of mpileup, to get methylation matrices for each motif. There are also functions to combine the CG and GC matrices, and then work with them.

## Prerequisites

Getting single read matrices requires: 
1) Bed files with the position of motifs.
2) Bam file with the reads. 
3) GRanges object with the region of interest.
4) GRanges object with unique non-overlapping CG and GC motifs in the genome

The bed files can be created from a file containing the genomic sequence. So first we read in the paths to the genome sequence, the bam file and regions of interest.

```r
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
```

To create the bed files with the CG or GC motifs, use the function **makeCGorGCbed**.

```r
makeCGorGCbed(genomeFile,"CG")
makeCGorGCbed(genomeFile,"GC")
```

These will save a bed file for that motif to the same directory as the genome file. Then in future one just needs to give the path to these files.


```r
#bedFileCG="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.CG.bed"
bedFileCG<-system.file("extdata",
                       "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CG.bed",
                       package="methMatrix",mustWork=TRUE)
#bedFileGC="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.GC.bed"
bedFileGC<-system.file("extdata",
                       "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.GC.bed",
                       package="methMatrix",mustWork=TRUE)
```

One also has to create a GRanges object with a list of unique, non-overlapping CG and GC motifs in the genome. This can be done with the **nanodsmf::findGenomeMotifs** function, and then saved as an RDS. Note that this function takes a long time to run (3.5h on the C. elegans genome, for example)!


```r
genomeMotifGR<-nanodsmf::findGenomeMotifs(genomeFile)
```

Once the RDS is saved it can be loaded for all future analyses.


```r
#genomeMotifGR<-readRDS("/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJ#NA13758.WS250.genomic.CGGC_motifs.RDS")
genomeMotifGR<-readRDS(system.file("extdata",
              "genome/c_elegans.PRJNA13758.WS250.genomic_Xchr.CGGC_motifs.RDS",
                      package="methMatrix", mustWork=TRUE))
```

Note that getting a unique, non-overlapping list of CG and GC motifs requires some manipulation, since triplet motifs such as GCG or CGC and longer runs CGCG.. or GCGC... will have multiple overlapping CG and GC motifs within them. To deal with this we treat isolated triplets as a single motif. Longer runs of CGCG.. or GCGC... are artifically split. If their length is even, they are split into doublets and if their length is odd, one triplet is split off and the rest is split into doublets. Unlike NOME-seq, where you might want to distinguish between CG and GC motifs because they convey different types of information (and therefore GCG motifs are discarded), in DSMF, both motifs give an indication of DNA accessibility, therefore this approach allows simplification of the data without loss of overlapping sites.

## Extracting methylation matrix from bam file

Use the **getReadMatrix** to extract methylation status in each read for a specific motif. Each motif is performed separately. The matrices will contain reads as rows and C positions as columns. 


```r
matCG<-getReadMatrix(bamFile=bamFile,genomeFile=genomeFile,bedFile=bedFileCG,region=regionGR, samtoolsPath="~/miniconda3/bin/")
matGC<-getReadMatrix(bamFile,genomeFile,bedFileGC,regionGR,samtoolsPath="~/miniconda3/bin/")

matCG[2:5,18:26]
```

```
##                                               8069559 8069564 8069565 8069599
## M02442:245:000000000-C5LKB:1:1101:11219:26350       1       0       1       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1       0       1       0
##                                               8069600 8069610 8069611 8069617
## M02442:245:000000000-C5LKB:1:1101:11219:26350       1       0       1       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1       0       1       0
##                                               8069618
## M02442:245:000000000-C5LKB:1:1101:11219:26350       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1
```

```r
matGC[2:5,8:16]
```

```
##                                               8069484 8069519 8069520 8069652
## M02442:245:000000000-C5LKB:1:1101:11219:26350       0       1       0       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682       0      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       0       1       0       0
##                                               8069653 8069666 8069667 8069694
## M02442:245:000000000-C5LKB:1:1101:11219:26350       0       0       0       1
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       1
## M02442:245:000000000-C5LKB:1:1101:13980:22059       0       0       0       1
##                                               8069695
## M02442:245:000000000-C5LKB:1:1101:11219:26350       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       0
```

This are very sparse matrices because forward and reverse reads will have C calls in different positions even within the same motif (C is shifted by 1 between the two strands). To simplify the data and get a more workable matrix we perform the following steps: 

1) The data is "destranded" by taking an average of the methylation call at both positions within a motif. To do this, NAs are ignored and average of the call at both positions in a motif (or all three in the case of a triplet) is calculated for each read. Motif positions that had only NAs are called as NAs.

2) To combine CG and GC matrices, positions that overlap between both matrices (GCGorCGC motifs will be called in both matrices) need to be dealt with. Again, the average of methylation calls for any read at a given motif is caclulated from both matrices while ignoring NAs. If all positions are NAs, they are scored as NAs. 

3) The averaged overlapping motifs are combined with unique motifs from each matrix into a single matrix

4) The motifs are reduced to a single bp as follows: CG motifs take the position of the first bp in the motif, GC motifs take the position of the second bp in the motif, and GCGorCGC motifs take the position of the third bp in the motif. 

The output of the **combineCGandGCmatrices** function is a single matrix with reads as rows and C positions within unique, non-overlapping CG, GC or GCGorCGC motifs as columns. The matrix contains values between 0 (non methylated) and 1 (methylated), as well as NAs for positions without a methylation call.


```r
methMat<-combineCGandGCmatrices(matCG,matGC,regionGR,genomeMotifGR)
methMat[2:5,10:18]
```

```
##                                               8069530 8069549 8069558 8069564
## M02442:245:000000000-C5LKB:1:1101:11219:26350       1       1       1       1
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1       1       1       1
##                                               8069599 8069610 8069617 8069631
## M02442:245:000000000-C5LKB:1:1101:11219:26350       1       1       0       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA      NA      NA      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0       0       0       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1       1       1       1
##                                               8069645
## M02442:245:000000000-C5LKB:1:1101:11219:26350       0
## M02442:245:000000000-C5LKB:1:1101:12625:13682      NA
## M02442:245:000000000-C5LKB:1:1101:13549:26394       0
## M02442:245:000000000-C5LKB:1:1101:13980:22059       1
```

```r
colSums(methMat,na.rm=T)
```

```
## 8069424 8069426 8069431 8069443 8069453 8069472 8069483 8069498 8069520 8069530 
##   241.0   191.5   240.0   261.0   406.0   395.0   217.0   440.0   345.0   454.0 
## 8069549 8069558 8069564 8069599 8069610 8069617 8069631 8069645 8069653 8069667 
##   440.0   478.0   460.0   462.0   465.0   452.0   448.0   373.0   299.0   258.0 
## 8069670 8069674 8069695 8069723 8069756 
##   420.0   433.0   259.0   491.0   476.0
```

This methylation matrix has more the look of a standard methylation matrix with consecutive sites having methylation scores, though the reads here are a bit short (2x100bp).

# Manipulating matrices to relative coordinates

The initial matrices can be converted from absolute genomic coordinates to coordinates relative to a reference point, like a transcription start site (TSS). First methylation data should be extracted from bam file with genomic ranges specifying the feature. These will be expanded to specify a window around the feature e.g +-250bp. We will load such matrices with genomic coordinates. The matrix files must all be listed in a csv file under the column heading "filename". 


```r
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
```

```
##     sample         region
## 1 dS02-182 WBGene00015947
##     sample         region
## 2 dS02-182 WBGene00015955
##     sample         region
## 3 dS02-182 WBGene00009620
##     sample         region
## 4 dS02-182 WBGene00009621
##     sample         region
## 5 dS02-182 WBGene00009801
##    sample         region
## 6 dS03-N2 WBGene00015947
##    sample         region
## 7 dS03-N2 WBGene00015955
##    sample         region
## 8 dS03-N2 WBGene00009620
##    sample         region
## 9 dS03-N2 WBGene00009621
##     sample         region
## 10 dS03-N2 WBGene00009801
```

```r
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
```

```
## [1] "plotting WBGene00015947"
```

```
## Warning: Removed 32025 rows containing missing values (geom_tile).
```

```
## Warning: Removed 13250 rows containing missing values (geom_tile).
```

```
## [1] "plotting WBGene00015955"
```

```
## Warning: Removed 6732 rows containing missing values (geom_tile).
```

```
## Warning: Removed 2109 rows containing missing values (geom_tile).
```

```
## [1] "plotting WBGene00009620"
```

```
## Warning: Removed 78182 rows containing missing values (geom_tile).
```

```
## Warning: Removed 26055 rows containing missing values (geom_tile).
```

```
## [1] "plotting WBGene00009621"
```

```
## Warning: Removed 575 rows containing missing values (geom_tile).
```

```
## Warning: Removed 2232 rows containing missing values (geom_tile).
```

```
## [1] "plotting WBGene00009801"
```

```
## Warning: Removed 21672 rows containing missing values (geom_tile).
```

```
## Warning: Removed 5775 rows containing missing values (geom_tile).
```


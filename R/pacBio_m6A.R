###### Functions to create matrices from fiber-seq m6A data #########
## Paper: https://www.science.org/doi/10.1126/science.aaz1646


#' Convert fiberseq bed file to bigwig
#'
#' Processes data from bed file that is the output of
#' Andrew Stergachis fiberseq scripts where each row is a
#' fiber and the final columns indicate methylated positions
#' in that fiber. Output is: 1) the original GRanges with all invalide fibers removed
#' and with additional column denoting fraction methylation, 2) A new GRanges
#' with all valid A/T postions and their fraction methylation.
#' @param bedGR GRanges object of fiberseq output (bed) imported into R
#' @param genome BSgenome object for your organism with UCSC style seqinfo (Default
#' is C. elegans ce11)
#' @param ATpositionGR GRanges object with position of all As and Ts in genome
#' @param minSubreadCov Minimum coverage of the read with subreads (default=10)
#' @param minReadCov Minimum read coverage per genomic position (default=5)
#' @return List of two GRanges: The first is the original bedGR but with additonal
#' metadata columns for fraction methylation per fiber. The second contains all the
#' A/T sites along the genome for which a valid fraction methylation could be calculated.
#' @export
fiberseqBedToBigwig<-function(bedGR,
                genome=BSgenome.Celegans.UCSC.ce11::Celegans,
                ATpositionGR=NULL,
                minSubreadCov=10,
                minReadCov=5){
  print("Extracting methylation position data...")
  bedGR<-bedGR[bedGR$score>=minSubreadCov]

  # use blck widths (0 or 1) to get length of vector of 1s per fiber
  blcks<-unlist(bedGR$blocks)
  blckrle<-rle(IRanges::width(blcks))
  blcklengths<-blckrle$lengths[blckrle$values==1]
  rm(list=c("blckrle"))
  # remove blcks with width 0 (these just denote the first and last position on the fiber, but not methylation info)
  blcks<-blcks[IRanges::width(blcks)==1]

  #make granges, adding seqnames replicated by blcklengths
  sn<-as.vector(GenomicRanges::seqnames(bedGR))
  starts<-rep(GenomicRanges::start(bedGR),blcklengths)

  gr<-GenomicRanges::GRanges(seqnames=rep(sn,blcklengths),
              ranges=IRanges::IRanges(start=IRanges::start(blcks)+starts-1,
                                      width=1))

  rm(list=c("starts","sn","blcks"))
  # get seqinfo data from genome object to be able to calculate coverage
  GenomeInfoDb::seqlevelsStyle(gr)<-"UCSC"
  GenomeInfoDb::seqlevels(gr)<-GenomeInfoDb::seqlevels(genome)
  GenomeInfoDb::seqinfo(gr)<-GenomeInfoDb::seqinfo(genome)
  gr<-GenomicRanges::trim(gr)

  # calculate coverage of methylation
  print("Calculating coverage...")
  mecov<-GenomicRanges::coverage(gr)
  rm(list=c("gr")) # remove big gr object

  GenomeInfoDb::seqlevelsStyle(bedGR)<-"UCSC"
  GenomeInfoDb::seqlevels(bedGR)<-GenomeInfoDb::seqlevels(genome)
  GenomeInfoDb::seqinfo(bedGR)<-GenomeInfoDb::seqinfo(genome)
  allcov<-GenomicRanges::coverage(bedGR)

  fracMe<-mecov/allcov

  rm(list=c("mecov"))
  print("limiting to AT positions only...")
  # get location of AT bases
  if(is.null(ATpositionGR)){
    ATpositionGR<-makeATgrObj(genome)
  }
  atfracMe<-GenomicRanges::binnedAverage(ATpositionGR,fracMe,varname="fracMe")
  rm(list=c("fracMe"))

  #look at percentage of positions methylated per fiber
  bedGR$totalAT<-GenomicRanges::countOverlaps(bedGR,ATpositionGR)
  bedGR$fracMePerFiber<-blcklengths/bedGR$totalAT
  # TODO: filter out fibers with > 20-30% methylation?

  # remove regions that have less the minReadCov reads
  allcovgr<-GenomicRanges::bindAsGRanges(allcov)
  removeLowCov<-allcovgr[allcovgr$V1<minReadCov]
  ol<-GenomicRanges::findOverlaps(atfracMe,removeLowCov)
  atfracMe<-atfracMe[-S4Vectors::queryHits(ol)]

  # remove regions that have NA values
  idx<-is.na(atfracMe$fracMe)
  print(paste0(sum(idx)," AT positions are NA"))
  atfracMe<-atfracMe[!is.na(atfracMe$fracMe)]

  #GenomeInfoDb::seqinfo(atfracMe)<-GenomeInfoDb::seqinfo(genome)
  return(list(bedGR,atfracMe))
}




#' Make GRanges object of all As and Ts in genome
#'
#' Make GRanges object of all As and Ts, with runs of A/Ts
#' unmerged.
#' @param genome BSgenome object for your organism (default is C.elegans)
#' @return GRanges object
#' @export
makeATgrObj<-function(genome=BSgenome.Celegans.UCSC.ce11::Celegans){
  dnass<-BSgenome::getSeq(genome)
  As<-grangesutils::mIdxToGR(Biostrings::vmatchPattern("A",dnass))
  Ts<-grangesutils::mIdxToGR(Biostrings::vmatchPattern("T",dnass))
  #ATranges<-GenomicRanges::reduce(c(As,Ts),ignore.strand=T)
  ATranges1<-sort(c(As,Ts))
  GenomeInfoDb::seqinfo(ATranges1)<-GenomeInfoDb::seqinfo(genome)
  return(ATranges1)
}


#' Convert fiberseq bed file to bigwig
#'
#' Processes data from bed file that is the output of
#' Andrew Stergachis fiberseq scripts where each row is a
#' fiber and the final columns indicate methylated positions
#' in that fiber. Output is: 1) the original GRanges with all invalide fibers removed
#' and with additional column denoting fraction methylation, 2) A new GRanges
#' with all valid A/T postions and their fraction methylation.
#' @param bedGR GRanges object of fiberseq output (bed) imported into R
#' @param genome BSgenome object for your organism with UCSC style seqinfo (Default
#' is C. elegans ce11)
#' @param ATpositionGR GRanges object with position of all As and Ts in genome
#' @param minSubreadCov Minimum coverage of the read with subreads (default=10)
#' @param minReadCov Minimum read coverage per genomic position (default=5)
#' @param regionGR Genomic ranges object for region of interest
#' @return List of two GRanges: The first is the original bedGR but with additonal
#' metadata columns for fraction methylation per fiber. The second contains all the
#' A/T sites along the genome for which a valid fraction methylation could be calculated.
#' @export
fiberseqBedToMatrix<-function(bedGR,
                              genome=BSgenome.Celegans.UCSC.ce11::Celegans,
                              ATpositionGR=NULL,
                              minSubreadCov=10,
                              minReadCov=5,
                              regionGR=NULL){
  print("Extracting methylation position data...")
  if(!is.null(regionGR)){
    ol<-GenomicRanges::findOverlaps(regionGR,bedGR,ignore.strand=T,type="within")
    bedGR<-bedGR[S4Vectors::subjectHits(ol)]
  }
  bedGR<-bedGR[bedGR$score>=minSubreadCov]

  # use blck widths (0 or 1) to get length of vector of 1s per fiber
  blcks<-unlist(bedGR$blocks)
  blckrle<-rle(IRanges::width(blcks))
  blcklengths<-blckrle$lengths[blckrle$values==1]
  rm(list=c("blckrle"))
  # remove blcks with width 0 (these just denote the first and last position on the fiber, but not methylation info)
  blcks<-blcks[IRanges::width(blcks)==1]

  #make granges, adding seqnames replicated by blcklengths
  sn<-as.vector(GenomicRanges::seqnames(bedGR))
  starts<-rep(GenomicRanges::start(bedGR),blcklengths)
  readNames<-rep(bedGR$name,blcklengths)

  gr<-GenomicRanges::GRanges(seqnames=rep(sn,blcklengths),
                             ranges=IRanges::IRanges(start=IRanges::start(blcks)+starts-1,
                                                     width=1))

  gr$readNames<-readNames
  rm(list=c("starts","sn","blcks"))
  # get seqinfo data from genome object to be able to calculate coverage
  GenomeInfoDb::seqlevelsStyle(gr)<-"UCSC"
  GenomeInfoDb::seqlevels(gr)<-GenomeInfoDb::seqlevels(genome)
  GenomeInfoDb::seqinfo(gr)<-GenomeInfoDb::seqinfo(genome)
  gr<-GenomicRanges::trim(gr)

  print("getting all AT positions...")
  # get location of all AT bases
  if(is.null(ATpositionGR)){
    ATpositionGR<-makeATgrObj(genome)
  }

  # make a matrix with all A/T positions in the region
  ol<-GenomicRanges::findOverlaps(ATpositionGR,regionGR)
  dfall<-data.frame(ATpositionGR[S4Vectors::queryHits(ol)])
  matall<-matrix(data=0,nrow=length(bedGR),ncol=nrow(dfall))
  rownames(matall)<-bedGR$name
  colnames(matall)<-dfall$start

  # remove any methylated positions that are not A/T
  ol<-GenomicRanges::findOverlaps(ATpositionGR,gr)
  gr<-gr[S4Vectors::subjectHits(ol)]
  # remove parts of the reads that are not within the
  # region of interest
  ol<-GenomicRanges::findOverlaps(regionGR,gr)
  gr<-gr[S4Vectors::subjectHits(ol)]
  if(length(bedGR)>=minReadCov){
    # convert to dataframe and reshape to wide format
    df<-data.frame(sort(gr))
    df$methylation<-1
    df1<-df %>% tidyr::pivot_wider(id_cols=readNames,names_from=start,values_from=methylation)
    # convert to matrix with rownames=reads and colnames=position
    methMat<-as.matrix(df1[,-1])
    row.names(methMat)<-df1$readNames
    methMat[is.na(methMat)]<-0
    #copy data into matrix with all A/T positions in the region
    matall[rownames(methMat),colnames(methMat)]<-methMat
  } else {
    matall<-NULL
  }
  print(dim(matall))
  return(matall)
}

# getFiberSeqMatrices<-function(sampleTable, bedGR,
#                               genome=BSgenome.Celegans.UCSC.ce11::Celegans,
#                               regionGRs,
#                               ATpositionGR=NULL,
#                               minSubreadCov=10,
#                               minReadCov=5,
#                               path=".",
#                               convRatePlots=FALSE, nThreads=1,
#                               overwriteMatrixLog=FALSE){
#
# }



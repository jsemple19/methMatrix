###### Functions to create matrices from fiber-seq m6A data #########
## Paper: https://www.science.org/doi/10.1126/science.aaz1646


#' Convert fiberseq bed file to bigwig
#'
#' Processes data from bed file that is the output of
#' Andrew Stergachis fiberseq scripts where each row is a
#' fiber and the final columns indicate methylated positions
#' in that fiber
#' @param bedGR GRanges object of fiberseq output (bed) imported into R
#' @param genome BSgenome object for your organism with UCSC style seqinfo
#' @param ATpositionGR GRanges object with position of all As and Ts in genome
#' @param minSubreadCov Minimum coverage of the read with subreads (default=10)
#' @param minReadCov Minimum read coverage per genomic position (default=5)
#' @return GRanges with score column containing % methylation
#' along genome
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
  hist(bedGR$fracMePerFiber,breaks=100,xlim=c(0,1))
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


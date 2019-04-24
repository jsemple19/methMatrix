
#' Create bed file of CpG or GpC motifs
#'
#' @param genomeFile String with path to fasta file with genome sequence
#' @param CpGorGpC String with motif to be found (either "CG" or "GC")
#' @return A bed file with locations of the C in all occurrences of this motif on both strands will be written to
#' the same directory as genomeFile with the extension .CpG.bed or .GpC.bed
#'
#' @export
makeCGorGCbed<-function(genomeFile,CGorGC) {
  if (! CGorGC %in% c("CG","GC")) {
    print("The variable 'CGorGC' must have the value 'CG' or 'GC'")
  }
  genome<-Biostrings::readDNAStringSet(genomeFile)
  Cmotifs<-Biostrings::vmatchPattern(CGorGC,genome)

  # create gr for Cs on positive strand
  gr<-GenomicRanges::GRanges(seqnames=S4Vectors::Rle(names(Cmotifs),sapply(Cmotifs,length)),
                             ranges=unlist(Cmotifs), strand="+")
  allGR<-gr
  # create gr for Cs on negative strand
  gr<-GenomicRanges::GRanges(seqnames=S4Vectors::Rle(names(Cmotifs),sapply(Cmotifs,length)),
                             ranges=unlist(Cmotifs), strand="-")
  allGR<-c(allGR,gr)
  # resize to cover only Cs
  allGR<-GenomicRanges::resize(allGR,width=1,fix=ifelse(CGorGC=="CG","start","end"))
  allGR<-sort(allGR,ignore.strand=T)

  #export as bed file to genomeFile directory
  outDir<-dirname(genomeFile)
  outFile<-gsub("fa",paste0(CGorGC,".bed"),basename(genomeFile))
  rtracklayer::export.bed(allGR,con=paste0(outDir,"/",outFile))
}


#' Extract methylation matrix from alignment file
#'
#' @param bamFile String with path to bam file with alignments of reads to genome
#' @param genomeFile String with path to fasta file with genome sequence
#' @param bedFile String with path to .bed file with locations of Cs to be evaluated
#' @param region Genomic range of region for which to extract reads. It can also be a string
#' denoting the region, e.g. "X:8069396-8069886"
#' @return A matrix with reads as row names and C positions as column names.
#' Every position in matrix has a value of 0 (not methylated), 1 (methylated) or
#' NA (undetermined)
#'
#' @export
getReadMatrix<-function(bamFile,genomeFile,bedFile,region) {
  if (class(region)=="GRanges") {
    region<-paste0(GenomeInfoDb::seqnames(region),":",IRanges::start(region),"-",IRanges::end(region))
  }
  # use samtools mpileup to call C methylation
  tab<-system(paste0("samtools mpileup -f ",genomeFile," -l ",bedFile," -r ",region,
                "  --output-QNAME ",bamFile),intern=T)
  #convert output to data frame
  tab<-lapply(tab,strsplit,"\t")
  tab<-lapply(1:length(tab),function(x){t(tab[[x]][[1]])})
  tab<-as.data.frame(do.call(rbind,tab),stringsAsFactors=F)
  colnames(tab)<-c("chr","start","ref","count","matches","BQ","reads")
  # make empty matrix of reads x C positions
  allReads<-unique(sort(unlist(sapply(tab$reads,strsplit,split=","),use.names=FALSE)))
  allPos<-tab$start
  mat<-matrix(data=NA,nrow=length(allReads),ncol=length(allPos))
  colnames(mat)<-allPos
  rownames(mat)<-allReads
  # scroll through tab to extract methylation values in strand aware way
  #plusStrand<-tab[tab$ref=="c" | tab$ref=="C",]
  #minusStrand<-tab[tab$ref=="g" | tab$ref=="G",]
  #matchList<-stringr::str_locate_all(tab[,"matches"],"\\.")
  for (line in 1:nrow(tab)) { #nrow(tab)
    df<-pileupToMethStatus(tab[line,])
    mat[df$reads,df$pos]<-df$meth
  }
  return(mat)
}

#' Extract methylation calls from single line of pileup file
#'
#' @param pileupLine One line data frame from pileup file. Columns have been named: c("chr","start","ref","count","matches","BQ","reads")
#' @return A data frame with C position, read name and methylation call  0 (not methylated), 1 (methylated)
pileupToMethStatus<-function(pileupLine) {
  # if read matches forward strand
  match<-parseMatchString(pileupLine$matches)
  if(nchar(match)!=pileupLine$count) {
    "problem with match string"
  }
  df=NULL
  posMeth=NULL
  posUnmeth=NULL
  if (pileupLine$ref=="c" | pileupLine$ref=="C") {
    matchList<-stringr::str_locate_all(match,"[\\.]")[[1]]
    if(dim(matchList)[1]>0) {
      posMeth<-rep(pileupLine$start,dim(matchList)[1])
      readsMeth<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
    matchList<-stringr::str_locate_all(match,"[T]")[[1]]
    if(dim(matchList)[1]>0) {
      posUnmeth<-rep(pileupLine$start,dim(matchList)[1])
      readsUnmeth<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
  }
  # if read matches reverse strand
  if (pileupLine$ref=="g" | pileupLine$ref=="G") {
    matchList<-stringr::str_locate_all(match,"[,]")[[1]]
    if(dim(matchList)[1]>0) {
      posMeth<-rep(pileupLine$start,dim(matchList)[1])
      readsMeth<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
    matchList<-stringr::str_locate_all(match,"[a]")[[1]]
    if(dim(matchList)[1]>0) {
      posUnmeth<-rep(pileupLine$start,dim(matchList)[1])
      readsUnmeth<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
  }
  #prepare df for export
  if (!is.null(posMeth)) {
    methdf<-data.frame(pos=posMeth,reads=readsMeth,meth=1,stringsAsFactors=F)
    df<-methdf
  }
  if (!is.null(posUnmeth)) {
    unmethdf<-data.frame(pos=posUnmeth,reads=readsUnmeth,meth=0,stringsAsFactors=F)
    if (exists("df")) {
      df<-rbind(df,unmethdf)
    } else {
      df<-unmethdf
    }
  }
  return(df)
}

#a dot stands for a match to the reference base on the forward strand,
#a comma for a match on the reverse strand, a '>' or '<' for a reference skip,
#`ACGTN' for a mismatch on the forward strand and
#`acgtn' for a mismatch on the reverse strand.
# `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position
# `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference

# * does not add character
# $ indicates end of read, it adds one character and should be removed
# ^ indicates a read start mark followed by mapping quality. two characters should be removed
# +[1-9]+[]


#' Clean up match string from pileup file
#'
#' @param match Match string from pileup file.
#' @return A cleaned up match string with no indels or start/end line info
parseMatchString<-function(match) {
  match<-gsub("\\$","",match) # remove end of read mark
  match<-gsub("\\^.?","",match) # remove start of read mark + read qual char
  indelCount<-stringr::str_extract_all(match,"[0-9]+")[[1]]
  for (i in indelCount) { # remove indel info around match
    m<-paste0("[\\+-]",i,".{",i,"}")
    match<-stringr::str_replace(match,m,"")
  }
  return(match)
}

#' Convert methylation matrix to genomic ranges
#'
#' @param mat Methylation matrix with reads x Cposition
#' @param gr GRanges object for which the methylation matrix was made
#' @return genomic ranges of Cpositions  with mcols containing the reads as columns
matToGR<-function(mat,gr) {
  matgr<-GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gr),
                         IRanges::IRanges(start=as.integer(colnames(mat)),width=1),
                         strand="*")
  GenomicRanges::mcols(matgr)<-t(mat)
  return(matgr)
}


#' Check if vector is only NAs
#'
#' @param vec Vector of values to be checked
#' @examples
#' allNAs(c(NA, NA, NA))
#' allNAs(c(NA,1,NA))
#' allNAs(c(1,2,3))
#' allNAs(c())
#' @return boolean TRUE or FALSE
allNAs<-function(vec) {
  if (sum(is.na(vec)) == length(vec))  {
    returnVal=TRUE
  } else {
    returnVal=FALSE
  }
  return(returnVal)
}


genomeFile="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"
bedFileCG="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.CG.bed"
bedFileGC="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.GC.bed"
bamFile="/Users/semple/Documents/MeisterLab/sequencingData/20190206_dSMFv016v020_N2gw_spike100pe/bwaMeth/test20190405/aln/dS16N2gw_20190206.noOL.bam"
amplicons=readRDS("/Users/semple/Documents/MeisterLab/dSMF/PromoterPrimerDesign/usefulObjects/ampliconGR.RDS")

#makeCGorGCbed(genomeFile,"CG")
#makeCGorGCbed(genomeFile,"GC")

region="I:10560-11095"
genomeFile="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"
region=nanodsmf::ucscToWbGR(amplicons[1])
names(GenomicRanges::mcols(region))<-"ID"
mCG<-getReadMatrix(bamFile,genomeFile,bedFileCG,region)
mGC<-getReadMatrix(bamFile,genomeFile,bedFileGC,region)

#gnmgr<-nanodsmf::findGenomeMotifs(genomeFile)
gnmgr<-readRDS("/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS260.genomicCGGC_motifs.RDS")

# convert mCG and mGC to genomic ranges with transposed matrix
mCGgr<-matToGR(mCG,region)
mGCgr<-matToGR(mGC,region)

#subset gnmgr by maxInterval
regGCCG<-IRanges::subsetByOverlaps(gnmgr,region)

# use gnmgr subset to apply sum over mCG and mGC
CGreads<-colnames(GenomicRanges::mcols(mCGgr))
GCreads<-colnames(GenomicRanges::mcols(mGCgr))
library(magrittr)

cg<-nanodsmf::applyGRonGR(regGCCG,mCGgr,CGreads,sum,na.rm=T)
gc<-nanodsmf::applyGRonGR(regGCCG,mGCgr,GCreads,sum,na.rm=T)

ol<-IRanges::findOverlaps(cg,gc)


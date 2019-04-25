
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
  # convert output to data frame
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
  for (line in 1:nrow(tab)) {
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

# a dot stands for a match to the reference base on the forward strand,
# a comma for a match on the reverse strand, a '>' or '<' for a reference skip,
# ACGTN for a mismatch on the forward strand and
# acgtn for a mismatch on the reverse strand.
# \\+[0-9]+[ACGTNacgtn]+ indicates there is an insertion between this reference position and the next reference position
# -[0-9]+[ACGTNacgtn]+ represents a deletion from the reference

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



#' Combine CG and GC methylation matrices
#'
#' Methylation matrices are tricky to merge because:
#'
#' 1. Forward and reverse reads will have methylation calls at different positions
#' even if they are part of the same motif because the C is shifted by 1 between
#' the two strands. Therefore we "destrand" the reads by averaging all methylation
#' calls within the same motif (and ingoring NAs). The genomic positions in the final merged
#' matrix will be on the position of the 1st bp of the motif for CG motifs and the genomic
#' position of the 2nd bp of the motif for GC motifs.
#'
#' 2. GCG/CGC motifs will have methylation calls from three Cs and it is not clear
#' whether the middle C should be part of one motif or the other. Therefore we consider
#' such triplets as a single motif and average the methylation calls from all three Cs.
#' The genomic position used in the final merged matrix will be that of the middle bp
#' of the motif.
#'
#' 3. Longer GCGC/CGCG runs have multiple overlapping motifs layered within them. To deal
#' with that we split these runs into 2-3bp motifs. If the length of the run is an even
#' number, it is split into doublets and considered as simple CG or GC motifs. If the
#' length of the run is an odd number, one triplet is created, and the rest is split into
#' doublets. This creates a set of unique non-overlapping motifs in the genome that are either
#' CG or GC 2bp motifs, or a GCG/CGC 2bp motif. This unique set is created before hand
#' with the nanodsmf::findGenomeMotifs function and can be stored in a RDS file. Doublet
#' and triplet motifs created from the run are treated the same as isolated doublet and triplet
#' motifs, as described in points 1. and 2. above.
#'
#' @param matCG Methylation matrix (reads x positons) of C positions within CG motifs
#' @param matGC Methylation matrix (reads x positons) of C positions within GC motifs
#' @param regionGR GRanges object of region used to make methylation matrices
#' @param gnmMotifGR GRanges object with all unique CG/GC/GCGorCGC motifs in genome
#' @return Merged methylation matrix
#' @export
combineCGandGCmatrices<-function(matCG,matGC,regionGR,gnmMotifGR){
  # convert matCG and matGC to genomic ranges with transposed matrix (positions x reads)
  matCGgr<-matToGR(matCG,regionGR)
  matGCgr<-matToGR(matGC,regionGR)

  #subset gnmMotifGR by regionGR to get motifs that should be present in the matrices
  regGCCG<-IRanges::subsetByOverlaps(gnmMotifGR,regionGR)

  # get vector of read names for each gr
  CGreads<-colnames(GenomicRanges::mcols(matCGgr))
  GCreads<-colnames(GenomicRanges::mcols(matGCgr))

  # use gnmMotifGR subset to "destrand" CG and GC calls by summing positions within motifs
  cg<-nanodsmf::applyGRonGR(regGCCG,matCGgr,CGreads,sum,na.rm=T)
  gc<-nanodsmf::applyGRonGR(regGCCG,matGCgr,GCreads,sum,na.rm=T)

  # find gr that overlap between cg and gc calls
  ol<-IRanges::findOverlaps(cg,gc)
  # find reads that are in both grs (GCGorCGC motifs)
  idxCGinGC<-CGreads %in% GCreads
  idxGCinCG<-GCreads %in% CGreads
  # average values from both matrices at these overlapping sites (NAs are ignored, but
  # kept if not real values is found)
  m1<-as.matrix(GenomicRanges::mcols(cg)[S4Vectors::queryHits(ol),CGreads[idxCGinGC]])
  m2<-as.matrix(GenomicRanges::mcols(gc)[S4Vectors::subjectHits(ol),GCreads[idxGCinCG]])
  mAvr<-ifelse(is.na(m1), ifelse(is.na(m2), NA, m2), ifelse(is.na(m2), m1, (m1 + m2)/2))
  GenomicRanges::mcols(cg)[S4Vectors::queryHits(ol),CGreads[idxCGinGC]]<-tibble::as_tibble(mAvr)

  # make nonoverlapping combined list
  allGR<-c(cg,gc[!idxGCinCG])
  # shrink gr to single bp
  CG1bp<-resize(allGR[allGR$context=="HCG"],width=1,fix="start")
  GC1bp<-resize(allGR[allGR$context=="GCH"],width=1,fix="end")
  GCGorCGC1bp<-resize(allGR[allGR$context=="GCGorCGC"],width=1,fix="center")
  # combine and sort
  allGR<-sort(c(CG1bp,GC1bp,GCGorCGC1bp))

  # convert back to matrix
  readNum<-dim(GenomicRanges::mcols(allGR))[2]-1 # don't include the "context" column
  methMat<-t(as.matrix(GenomicRanges::mcols(allGR)[,2:(readNum+1)]))
  colnames(methMat)<-GenomicRanges::start(allGR)
  # remove reads with aboslutly no methylation info
  methMat<-methMat[rowSums(is.na(methMat))!=dim(methMat)[2],]
  return(methMat)
}



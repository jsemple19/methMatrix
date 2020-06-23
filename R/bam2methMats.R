
#' Create bed file of CpG or GpC motifs
#'
#' @param genomeFile String with path to fasta file with genome sequence
#' @param CGorGC String with motif to be found (either "CG" or "GC")
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
  allGR<-GenomicRanges::sort(allGR,ignore.strand=T)

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
#' @param region Genomic range of region for which to extract reads. It can also
#' be a string denoting the region, e.g. "X:8069396-8069886"
#' @param samtoolsPath Path to samtools executable (default="") if not in unix $PATH
#' @param maxDepth Maximum number of reads to take for a given region (default=10,000)
#' @return A matrix with reads as row names and C positions as column names.
#' Every position in matrix has a value of 0 (non-converted C = methylated), 1 (converted C (T) = not methylated) or NA (undetermined)
#'
#' @export
getReadMatrix<-function(bamFile, genomeFile, bedFile, region, samtoolsPath="",
                        maxDepth=10000) {
  if (class(region)=="GRanges") {
    region<-paste0(GenomeInfoDb::seqnames(region), ":", IRanges::start(region),
                   "-", IRanges::end(region))
  }
  # use samtools mpileup to call C methylation
  tab<-system(paste0(samtoolsPath, "samtools mpileup -f ", genomeFile, " -l ",
                     bedFile, " -r ", region,
                "  --output-QNAME --max-depth ", maxDepth, " --min-BQ 8 ",
                " --ff UNMAP,QCFAIL ", bamFile),intern=T)
  if (length(tab)>0) {
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
      df<-pileupToConversionStatus(tab[line,])
      m1<-mat[df$reads,df$pos]
      m2<-df$Cconv
      mat[df$reads,df$pos]<-ifelse(is.na(m1), ifelse(is.na(m2), NA, m2),
                                   ifelse(is.na(m2), m1, (m1 + m2)))
    }
  } else {
    mat<-NULL
  }
  return(mat)
}



#' Extract methylation calls from single line of pileup file
#'
#' @param pileupLine One line data frame from pileup file. Columns have been named: c("chr","start","ref","count","matches","BQ","reads")
#' @return A data frame with C position, read name and C conversion call  0 (not converted (C) = methylated), 1 (converted (T) = not methylated)
#' @export
pileupToConversionStatus<-function(pileupLine) {
  # if read matches forward strand
  match<-parseMatchString(pileupLine$matches)
  if(nchar(match)!=pileupLine$count) {
    "problem with match string"
  }
  df=NULL
  posConv=NULL
  posUnconv=NULL
  if (pileupLine$ref=="c" | pileupLine$ref=="C") {
    matchList<-stringr::str_locate_all(match,"[\\.,]")[[1]]
    if(dim(matchList)[1]>0) {
      posConv<-rep(pileupLine$start,dim(matchList)[1])
      readsConv<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
    matchList<-stringr::str_locate_all(match,"[Tt]")[[1]]
    if(dim(matchList)[1]>0) {
      posUnconv<-rep(pileupLine$start,dim(matchList)[1])
      readsUnconv<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
  }
  # if read matches reverse strand
  if (pileupLine$ref=="g" | pileupLine$ref=="G") {
    matchList<-stringr::str_locate_all(match,"[\\.,]")[[1]]
    if(dim(matchList)[1]>0) {
      posConv<-rep(pileupLine$start,dim(matchList)[1])
      readsConv<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
    matchList<-stringr::str_locate_all(match,"[aA]")[[1]]
    if(dim(matchList)[1]>0) {
      posUnconv<-rep(pileupLine$start,dim(matchList)[1])
      readsUnconv<-unlist(strsplit(pileupLine$reads,","))[matchList[,1]]
    }
  }
  #prepare df for export
  first<-TRUE
  if (!is.null(posConv)) {
    convdf<-data.frame(pos=posConv,reads=readsConv,Cconv=0,stringsAsFactors=F)
    df<-convdf
    first<-FALSE
  }
  if (!is.null(posUnconv)) {
    unconvdf<-data.frame(pos=posUnconv,reads=readsUnconv,Cconv=1,stringsAsFactors=F)
    if (first==FALSE) {
      df<-rbind(df,unconvdf)
    } else {
      df<-unconvdf
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
#' @export
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
#' @export
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
#' @param genomeMotifGR GRanges object with all unique non-overlapping CG/GC/GCGorCGC motifs in genome
#' @return Merged methylation matrix
#' @export
combineCGandGCmatrices<-function(matCG,matGC,regionGR,genomeMotifGR){
  if (!is.null(matCG) & !is.null(matGC)) {
    # convert matCG and matGC to genomic ranges with transposed matrix (positions x reads)
    matCGgr<-matToGR(matCG,regionGR)
    matGCgr<-matToGR(matGC,regionGR)

    #subset genomeMotifGR by regionGR to get motifs that should be present in the matrices
    regGCCG<-IRanges::subsetByOverlaps(genomeMotifGR,regionGR,ignore.strand=TRUE)

    # get vector of read names for each gr
    CGreads<-colnames(GenomicRanges::mcols(matCGgr))
    GCreads<-colnames(GenomicRanges::mcols(matGCgr))

    # use genomeMotifGR subset to "destrand" CG and GC calls by summing positions within motifs
    cg<-grangesutils::applyGRonGR(regGCCG,matCGgr,CGreads,sum,na.rm=T)
    gc<-grangesutils::applyGRonGR(regGCCG,matGCgr,GCreads,sum,na.rm=T)

    # ensure no value in the matrix exceeds 1
    maxval1<-function(m1){
      ifelse(is.na(m1), NA, ifelse(m1>1, 1, m1))
    }
    GenomicRanges::mcols(cg)[,2:dim(GenomicRanges::mcols(cg))[2]]<-
      data.table::as.data.table(maxval1(as.matrix(GenomicRanges::mcols(cg)[,2:dim(GenomicRanges::mcols(cg))[2]])))
    GenomicRanges::mcols(gc)[,2:dim(GenomicRanges::mcols(gc))[2]]<-
      data.table::as.data.table(maxval1(as.matrix(GenomicRanges::mcols(gc)[,2:dim(GenomicRanges::mcols(gc))[2]])))

    # find gr that overlap between cg and gc calls
    ol<-IRanges::findOverlaps(cg,gc)

    if (length(ol)>=1) {
      # find reads that are in both grs (GCGorCGC motifs)
      idxCGinGC<-CGreads %in% GCreads
      idxGCinCG<-GCreads %in% CGreads
      # average values from both matrices at these overlapping sites (NAs are ignored, but
      # kept if no real values is found)
      m1<-as.matrix(GenomicRanges::mcols(cg)[S4Vectors::queryHits(ol),CGreads[idxCGinGC]])
      m2<-as.matrix(GenomicRanges::mcols(gc)[S4Vectors::subjectHits(ol),GCreads[idxGCinCG]])
      mAvr<-ifelse(is.na(m1), ifelse(is.na(m2), NA, m2), ifelse(is.na(m2), m1, (m1 + m2)/2))
      GenomicRanges::mcols(cg)[S4Vectors::queryHits(ol),CGreads[idxCGinGC]]<-tibble::as_tibble(mAvr)
      # make nonoverlapping combined list
      allGR<-c(cg,gc[-S4Vectors::subjectHits(ol)])
    } else {
      allGR<-c(cg,gc)
    }

    # shrink gr to single bp
    CG1bp<-GenomicRanges::resize(allGR[allGR$context=="HCG"],width=1,fix="start")
    GC1bp<-GenomicRanges::resize(allGR[allGR$context=="GCH"],width=1,fix="end")
    GCGorCGC1bp<-GenomicRanges::resize(allGR[allGR$context=="GCGorCGC"],width=1,fix="center")
    # combine and sort
    allGR<-GenomicRanges::sort(c(CG1bp,GC1bp,GCGorCGC1bp))

    # convert back to matrix
    readNum<-dim(GenomicRanges::mcols(allGR))[2]-1 # don't include the "context" column
    convMat<-t(as.matrix(GenomicRanges::mcols(allGR)[,2:(readNum+1)]))
    colnames(convMat)<-GenomicRanges::start(allGR)
    # remove reads with aboslutly no methylation info
    convMat<-convMat[rowSums(is.na(convMat))!=dim(convMat)[2],]
  } else {
    convMat<-NULL
  }
  return(convMat)
}




#' Get methylation frequency genomic ranges
#'
#' Takes MethylDackel methylation call bedgraph files for CpG CHH and CHG motifs and creates genomic ranges for CG, GC and all other Cs. Output is a list of these three genomic ranges objects.
#'
#' @param baseFileName base file name used in MethylDackel extract command (assumes _CpG.bedgraph, _CHH.bedgraph and _CHG.bedgraph extensions were later added by MethylDackel)
#' @param pathToMethCalls Path to folder in which MethylDackel output is found
#' @param motifFile Path to GRanges object with all unique CG/GC/GCGorCGC motifs in genome
#' @param minDepth Minimum read depth. Positions with fewer than this number are discarded (default=5)
#' @return Named list with CG, GC and C genomic ranges with metadata about methylation.
#' @export
getMethFreqGR<-function(baseFileName,pathToMethCalls,motifFile,minDepth=5) {
  # get genome motifs
  gnmMotifs<-readRDS(motifFile)
  GCmotifs<-gnmMotifs[gnmMotifs$context=="GCH"]
  CGmotifs<-gnmMotifs[gnmMotifs$context=="HCG"]
  GCGmotifs<-gnmMotifs[gnmMotifs$context=="GCGorCGC"]

  # make sure the path variable is set
  pathToMethCalls<-gsub("/$","",pathToMethCalls)
  if (! exists("pathToMethCalls")) {pathToMethCalls="."}
  # read in methyldackel output files
  methCG<-rtracklayer::import(paste0(pathToMethCalls,"/",baseFileName,"_CpG.bedGraph"),format="bedGraph")
  colnames(GenomicRanges::mcols(methCG))<-c("methPercent","methylated","nonMethylated")
  methCG$readDepth<-rowSums(cbind(methCG$methylated,methCG$nonMethylated))
  methCG<-methCG[methCG$readDepth>minDepth]
  methCHH<-rtracklayer::import(paste0(pathToMethCalls,"/",
                                      baseFileName,"_CHH.bedGraph"),format="bedGraph")
  colnames(GenomicRanges::mcols(methCHH))<-c("methPercent","methylated","nonMethylated")
  methCHG<-rtracklayer::import(paste0(pathToMethCalls,"/",
                                      baseFileName,"_CHG.bedGraph"),format="bedGraph")
  colnames(GenomicRanges::mcols(methCHG))<-c("methPercent","methylated","nonMethylated")
  methNonCG<-GenomicRanges::sort(c(methCHH,methCHG))

  # create methGC
  ol4<-IRanges::findOverlaps(methNonCG,GCmotifs)
  ol5<-IRanges::findOverlaps(methNonCG,GCGmotifs)
  GCidx<-c(S4Vectors::queryHits(ol4),S4Vectors::queryHits(ol5))
  methGC<-GenomicRanges::sort(methNonCG[GCidx])
  methGC$readDepth<-rowSums(cbind(methGC$methylated,methGC$nonMethylated))
  methGC<-methGC[methGC$readDepth>minDepth]

  # create methC
  methC<-methNonCG[-GCidx]
  methC$readDepth<-rowSums(cbind(methC$methylated,methC$nonMethylated))
  methC<-methC[methC$readDepth>minDepth]

  return(list(CG=methCG,GC=methGC,C=methC))
}



#' Convert Genomic ranges list to long data frame for plotting
#'
#' Methylation frequency data stored in  a list(by sample) of list (by C context) of genomic ranges is extracted into a data frame. By specifying a context type ("CG","GC" or "C") this type of data is extracted for all samples listed in "samples" can converted to a dataframe.
#'
#' @param methFreqGR A list (by sample) of a list (by C context) of genomic ranges for cytosine methylation frquency
#' @param samples a vector of sample names for which to extract the data
#' @param Ctype The type of sequence context for the cytosine _("CG","GC" or "C")
#' @return A long-form data frame with methylation frequency and counts at different sites in different samples
#' @export
grlToDf<-function(methFreqGR,samples,Ctype) {
  first<-TRUE
  for (sampleName in samples) {
    df<-as.data.frame(methFreqGR[[sampleName]][[Ctype]])
    df$sampleName<-sampleName
    if (first==FALSE) {
      alldf<-rbind(alldf,df)
    } else {
      alldf<-df
      first<-FALSE
    }
  }
  alldf$sampleName<-factor(alldf$sampleName)
  return(alldf)
}


#' Convert Genomic ranges list to single genomic range with sample methylation and total counts in metadata
#'
#' genomic ranges for CG and GC positions in each samples are combined and sorted. Then genomic ranges
#' from different samples are merged together putting methylation frequency (_M) and total read counts (_T) for each sample in the metadata.
#'
#' @param methFreqGR A list (by sample) of a list (by C context) of genomic ranges for cytosine methylation frquency
#' @param samples a vector of sample names for which to extract the data
#' @return A genomic ranges wtih all CG and GC positions found in the samples. The metadata columns contain methylation frequency (_M) and total read count (_T) for each sample
#' @export
combineCGGCgr<-function(methFreqGR,samples) {
  first<-TRUE
  for (sampleName in samples) {
    cg<-methFreqGR[[sampleName]][["CG"]]
    gc<-methFreqGR[[sampleName]][["GC"]]
    cggc<-GenomicRanges::sort(c(cg,gc))
    fractionMeth<-cggc$methylated/cggc$readDepth
    readDepth<-cggc$readDepth
    GenomicRanges::mcols(cggc)<-NULL
    GenomicRanges::mcols(cggc)[,paste0(sampleName,"_M")]<-fractionMeth
    GenomicRanges::mcols(cggc)[,paste0(sampleName,"_T")]<-readDepth

    if (first==FALSE) {
      allcggc<-GenomicRanges::merge(allcggc,cggc,all=T)
    } else {
      allcggc<-cggc
      first<-FALSE
    }
  }
  return(allcggc)
}



#' Find bisulfite conversion rate of Cs in non-methylated context
#'
#' Use bedfile with positons of Cs that are in non-methylated context to obtain stats
#' about the number of informative Cs and conversion status of Cs per read
#'
#' @param bamFile String with path to bam file with alignments of reads to genome
#' @param genomeFile String with path to fasta file with genome sequence
#' @param bedFileC String with path to .bed file with locations of Cs to be evaluated (for forward strand calls)
#' @param bedFileG String with path to .bed file with locations of Gs to be evaluated (for reverse strand calls)
#' @param regionGR Genomic range of region for which to extract reads. It can
#' also be a string denoting the region, e.g. "X:8069396-8069886"
#' @param samtoolsPath Path to samtools executable (default="") if not in unix $PATH
#' @return A data frame with the names of the reads, the count of the number of informative Cs
#' per region (not NAs), the maximum number of possible Cs in the regon, and
#' the fraction of the informative Cs which have been bisulfite converted
#' @export
poorBisulfiteConversion<-function(bamFile,genomeFile,bedFileC,bedFileG,
                                  regionGR,samtoolsPath="") {
  # calls on forward strand
  matC<-getReadMatrix(bamFile, genomeFile, bedFileC, regionGR, samtoolsPath)
  dfc<-data.frame(reads=row.names(matC),stringsAsFactors=F)
  dfc$informativeCs<-rowSums(!is.na(matC))
  dfc$totalCs<-dim(matC)[2]
  dfc$fractionConverted<-rowMeans(matC,na.rm=T)
  # calls on reverse strand
  matG<-getReadMatrix(bamFile, genomeFile, bedFileG, regionGR, samtoolsPath)
  dfg<-data.frame(reads=row.names(matG),stringsAsFactors=F)
  dfg$informativeCs<-rowSums(!is.na(matG))
  dfg$totalCs<-dim(matG)[2]
  dfg$fractionConverted<-rowMeans(matG,na.rm=T)
  #choose highest for each read (calls of 0, or close to 0 will be produced if reads on other strand from bedfile)
  df<-merge(dfc,dfg,by=c("reads"),all=TRUE)
  df[is.na(df)]<-0
  keepDFC<-df$fractionConverted.x>=df$fractionConverted.y
  df[!keepDFC,c("informativeCs.x","totalCs.x","fractionConverted.x")]<-
    df[!keepDFC,c("informativeCs.y","totalCs.y","fractionConverted.y")]
  df[,c("informativeCs.y","totalCs.y","fractionConverted.y")]<-NULL
  colnames(df)<-gsub("\\.x","",colnames(df))
  return(df)
}




#' Make directories
#'
#' @param path String with path to where the directories should be made
#' @param dirNameList Vector of strings with names of directories to create (can include multilevel directories)
#' @return Creates the directories listed in dirNameList
#' @examples
#' makeDirs(path=".",dirNameList=c("/txt","/rds/sample1"))
#' @export
makeDirs<-function(path,dirNameList=c()) {
  for (d in dirNameList) {
    if (!dir.exists(paste0(path,"/",d))){  # for alignments
      dir.create(paste0(path,"/",d), recursive=TRUE, showWarnings=FALSE)
    }
  }
}



#' Extract list of methylation matrices
#'
#' Input requires a sampleTable with two columns: FileName contains the path to
#' a bam file of aligned sequences. SampleName contains the name of the sample.
#' The function will return a table of all matrix-files for all samples and all
#' regions listed in regionGRs, together with some info about the matrices.
#' Matrices contain values between 0 (not methylated) and 1 (methylated), or NA (undefined)
#' @param sampleTable Table with FileName column listing the full path to bam files belonging to the samples listed in the SampleName column
#' @param genomeFile String with path to fasta file with genome sequence
#' @param regionGRs A genomic regions object with all regions for which matrices should be extracted. The metadata columns must contain a column called "ID" with a unique ID for that region.
#' @param regionType A collective name for this list of regions (e.g TSS or amplicons)
#' @param genomeMotifGR A GenomicRanges object with a unique set of non-overlapping CG, GC and GCGorCGC sites
#' @param minConversionRate Minimal fraction of Cs from a non-methylation context that must be converted to Ts for the read to be included in the final matrix (default=0.8)
#' @param maxNAfraction Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)
#' @param bedFilePrefix The full path and prefix of the bed file for C, G, CG and GC positions in the genome (i.e path and name of the file without the ".C.bed",".G.bed", ".CG.bed" or ".GC.bed" suffix). Defulat is NULL and assumes the bed file are in the same location as the genome sequence file.
#' @param path Path for output. "plots", "csv" and "rds" directories will be created here. Default is current directory.
#' @param convRatePlots Boolean value: should bisulfite conversion rate plots be created for each region? (default=FALSE)
#' @param nThreads number of threads for parallelisation
#' @param samtoolsPath Path to samtools executable (default="") if not in unix $PATH
#' @param overwriteMatrixLog Should matrixLog file be overwritten (in case of
#'  change in analysis or data), or should already computed matrices be used and
#'  script skips to next matrix (in case of premature termination of analysis)
#'  (default=FALSE)
#' @return A list (by sample) of lists (by regions) of methylation matrices
#' @export
getSingleMoleculeMatrices<-function(sampleTable, genomeFile, regionGRs,
                                    regionType, genomeMotifGR,
                                    minConversionRate=0.8, maxNAfraction=0.2,
                                    bedFilePrefix=NULL, path=".",
                                    convRatePlots=FALSE, nThreads=1,
                                    samtoolsPath="", overwriteMatrixLog=FALSE) {
  totalCs<-informativeCs<-fractionConverted<-i<-NULL
  #create pathnames to bedfiles
  if (is.null(bedFilePrefix)){
    bedFilePrefix=gsub("\\.fa","", genomeFile)
  }
  bedFileC=paste0(bedFilePrefix,".C.bed")
  bedFileG=paste0(bedFilePrefix,".G.bed")
  bedFileCG=paste0(bedFilePrefix,".CG.bed")
  bedFileGC=paste0(bedFilePrefix,".GC.bed")

  # make output directories
  makeDirs(path,c("csv", paste0("rds/methMats_",regionType)))
  if (convRatePlots==TRUE) {
    makeDirs(path,c("plots/informativeCsPlots",
                    "plots/conversionRatePlots"))
  }
  samples<-sampleTable$SampleName

  addSampleName<-ifelse(length(unique(samples))==1,paste0("_",samples[1]),"")

  matrixLog<-getMatrixLog(paste0(path,"/csv/MatrixLog_",regionType,
                                 addSampleName,"_log.csv"))

  #print(matrixLog)

  for (currentSample in samples) {
    print(currentSample)
    bamFile<-sampleTable$FileName[sampleTable$SampleName==currentSample]
    #lists to collect plots for all regions in a particular sample
    if (convRatePlots==TRUE) {
      informativeCsPlots=vector()
      conversionRatePlots=vector()
    }
    clst<-parallel::makeCluster(nThreads)
    doParallel::registerDoParallel(clst)
    pmatrixLog<-foreach::foreach(i=1:length(regionGRs),
                                 .combine=rbind,
                                 .packages=c("GenomicRanges",
                                             "S4Vectors")) %dopar% {
    #for (i in seq_along(regionGRs)) {
      # find appropriate line of matrixLog, and check if data already exists
      regionGR<-regionGRs[i]
      logLine<-data.frame(regionNum=i,filename=NA,sample=currentSample,
                          region=regionGR$ID, numCGpos=NA, numGCpos=NA,
                          numUniquePos=NA, CGreads=NA, GCreads=NA,
                          methMatReads=NA, goodConvReads=NA,
                          fewNAreads=NA,stringsAsFactors=F)
      alreadyDone<-(currentSample %in% matrixLog$sample &
                      regionGR$ID %in% matrixLog$region)
      #j<-which(matrixLog$sample==currentSample & matrixLog$region==regionGR$ID)
      if(!alreadyDone | overwriteMatrixLog==T){
        # get C conversion matrices
        matCG<-getReadMatrix(bamFile, genomeFile, bedFileCG, regionGR,
                           samtoolsPath)
        matGC<-getReadMatrix(bamFile, genomeFile, bedFileGC, regionGR,
                           samtoolsPath)
        convMat<-combineCGandGCmatrices(matCG,matGC,regionGR,genomeMotifGR)
        if(! is.null(dim(convMat))) {
          # combine CG and GC matrices and change conversion=1 to methylation=1
          methMat<-1-convMat
          # record number of reads in the matrices
          logLine[1,"numCGpos"]<-ifelse(!is.null(dim(matCG)[2]),
                                          dim(matCG)[2], 0)
          logLine[1,"numGCpos"]<-ifelse(!is.null(dim(matGC)[2]),
                                          dim(matGC)[2], 0)
          logLine[1,"numUniquePos"]<-ifelse(!is.null(dim(methMat)[2]),
                                            dim(methMat)[2], 0)
          logLine[1,"CGreads"]<-ifelse(!is.null(dim(matCG)[1]),
                                         dim(matCG)[1], 0)
          logLine[1,"GCreads"]<-ifelse(!is.null(dim(matGC)[1]),
                                         dim(matGC)[1], 0)
          logLine[1,"methMatReads"]<-ifelse(!is.null(dim(methMat)[1]),
                                            dim(methMat)[1], 0)

          # get bisulfite conversion stats for Cs in non-methylated context
          df<-poorBisulfiteConversion(bamFile, genomeFile, bedFileC, bedFileG,
                                    regionGR, samtoolsPath)
          removeReads<-df[df$fractionConverted<minConversionRate,"reads"]
          methMat<-methMat[!(rownames(methMat) %in% removeReads),]
          logLine[1,"goodConvReads"]<-ifelse(!is.null(dim(methMat)[1]),
                                             dim(methMat)[1],0)

          if (is.null(dim(methMat))) {
            next
          }
          if (convRatePlots==TRUE) {
          ## plot histogram of number informative Cs per read
            p<-ggplot2::ggplot(df,ggplot2::aes(x=informativeCs/totalCs)) +
              ggplot2::geom_histogram() +
              ggplot2::ggtitle(paste0(regionGR$ID,
                                    " Informative Cs per read (totalCs: ",
                                    df$totalCs[1]," )"))
            # save to file
            plotName<-paste0(path,"/plots/informativeCsPlots/infC_", regionType,
                             "_", currentSample,"_", regionGR$ID, ".pdf")
            ggplot2::ggsave(plotName, plot=p, device="pdf", width=7.25, height=5,
                          units="cm")
            informativeCsPlots<-c(informativeCsPlots,plotName)
            ## plot histogram of number of bisulfite converted Cs per read as
            ## fraction of informative Cs
            p<-ggplot2::ggplot(df,ggplot2::aes(x=fractionConverted)) +
              ggplot2::geom_histogram() +
              ggplot2::xlim(c(0,1)) +
              ggplot2::ggtitle(paste0(regionGR$ID,
                    " Bisulfite converted Cs per read out of informative Cs")) +
              ggplot2::geom_vline(xintercept=minConversionRate, col="red",
                                linetype="dashed")
            # save to file
            plotName<-paste0(path,"/plots/conversionRatePlots/convR_",
                             regionType, "_", currentSample,"_",regionGR$ID,"
                             .pdf")
            ggplot2::ggsave(plotName, plot=p, device="pdf", width=7.25, height=5,
                          units="cm")
            conversionRatePlots<-c(conversionRatePlots,plotName)
          }
          # count reads that do not cover at least 1-maxNAfraction of the cytosines
          logLine[1,"fewNAreads"]<-sum(rowMeans(is.na(methMat))<maxNAfraction)
          sink(type="message")
          print(paste(i,regionType,currentSample,regionGR$ID,sep=" "))
          sink()
          matName<-paste0(path,"/rds/methMats_", regionType,"/", currentSample,
                        "_", regionGR$ID, ".rds")
          saveRDS(methMat,file=matName)
          logLine[1,"filename"]<-matName
          # write intermediate data to file so that if it crashes one can restart
          sink(file=paste0(path, "/csv/MatrixLog_", regionType,
                           addSampleName, "_log.csv"), append=TRUE,
                           type="output")
          cat(paste(logLine,collapse=","), sep="\n")
          sink()
          #utils::write.csv(logLine,paste0(path, "/csv/MatrixLog_", regionType,
          #                                addSampleName, ".csv"),
          #               quote=F, row.names=F)
        }
        #print(matrixLog[j,])
        logLine
      }
    }
    print(pmatrixLog)
    if (convRatePlots==TRUE) {
      #TODO:combine PDF function
    }
  }
  #tidy and sort pmatrixLog
  pmatrixLog<-rbind(matrixLog,pmatrixLog)
  pmatrixLog<-pmatrixLog[rowSums(!is.na(pmatrixLog))>0,]
  pmatrixLog<-pmatrixLog[order(pmatrixLog$regionNum),]
  # write file with final data
  utils::write.csv(pmatrixLog,paste0(path, "/csv/MatrixLog_", regionType,
                            addSampleName, ".csv"), quote=F, row.names=F)
  #file.remove(paste0(path, "/csv/MatrixLog_", regionType,
  #                             addSampleName, "_log.csv"))
  return(pmatrixLog)
}


#' Read in an existing matrixLog, or create a new one
#'
#' Input requires matrixLog file name with correct path. If file exists it will be read
#' in, otherwise an empty data.frame will be created with the correct number of fields
#' for the samples (samples) and the regions (regionGRs). Teh function returns this
#' matrixLog data.frame
#' @param matrixLogFile Name of matrix log file complete with relative or absolut path
#' @param samples Vector with names of samples to be processed
#' @param regionGRs A genomic regions object with all regions for which matrices should be extracted. The metadata columns must contain a column called "ID" with a unique ID for that region.
#' @return A data frame with columns to record data about the single molecule matrices
getMatrixLog<-function(matrixLogFile){
  if (file.exists(matrixLogFile) & file.info(matrixLogFile)$size>0) {
    # this allows restarting
    matrixLog<-utils::read.csv(matrixLogFile, stringsAsFactors=F, header=T,
                               row.names=NULL)
  } else {
    writeLines(paste(c("regionNum","filename","sample","region","numCGpos",
                       "numGCpos","numUniquePos","CGreads","GCreads",
                       "methMatReads","goodConvReads", "fewNAreads"),
                     collapse=","), con=matrixLogFile)
    #log table to record number of reads in matrix at various steps
    matrixLog<-data.frame(regionNum=NA, filename=NA,sample=NA, region=NA,
                          numCGpos=NA, numGCpos=NA, numUniquePos=NA,
                          CGreads=NA, GCreads=NA, methMatReads=NA,
                          goodConvReads=NA, fewNAreads=NA,
                          stringsAsFactors=F)
  }
  return(matrixLog)
}




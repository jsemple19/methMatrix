############ functions for working with methylation matrix lists ##################



# the methylation matrices for individual genes by TSS have the following structure:
# a list (by sample) of lists of matrices (by TSS)
# e.g
# [1] sample1
#     [1] TSS1 matrix of reads x Cpositions
#     [2] TSS2 matrix of reads x Cpositions
# [2] sample2
#     [1] TSS1 matrix of reads x Cpositions
#     [2] TSS2 matrix of reads x Cpositions
#
# the matrices contain METHYLATION values: 0 = not methylated, 1 = methylated.
# to avoid confusion, keep them this way and only convert to dSMF (1-methylation) for plotting
# or correlations

## necessary libraries
library(rtracklayer)
#library("RColorBrewer")
#library(GGally)
#library(gridExtra)
library(dplyr)
library(tidyr)
library(ggpubr)

# load example methMatrices
#methMats<-readRDS("/Users/semple/Documents/MeisterLab/sequencingData/20181119_dSMFv002-4amp_N2-182/dSMFamplicon/methylation_calls/allSampleRelCoordMats_TSS.rds")
#methMats<-readRDS("/Users/semple/Documents/MeisterLab/sequencingData/20181119_dSMFv002-4amp_N2-182/dSMFamplicon/methylation_calls/allSampleMergedMats_TSS.rds")
#save(methMats, file="data/methMats.RData")
#TSS<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/PromoterPrimerDesign/usefulObjects/ampliconMaxTSSgr.RDS")
#names(mcols(TSS))<-"ID"
#save(TSS, file="data/TSSgr.RData")

#' Extract matrices into a simple list
#'
#' @param methMats A table of filepaths to methylation matrices
#' (usually named by the genomic region they come from)
#' @param regionName A vector containing the names of regions of interest (default is all regions)
#' @param sampleName A vector containing the names of samples of interest (default is all samples)
#' @return A simple, one level list of methylation matrices for the samples and regions of interest.
#' Matrices will be named with both the sample and the region name joined by "__"
#' @examples
#' getMatrices(methMats)
#' #will return a simple list of all regions in all samples, changing the name of the matrices
#' getMatrices(methMats,regionName="WBGene00009234")
#' #will return a simple list of the matrices for that gene from all samples
#' getMatrices(methMats,regionName=c("WBGene00009234","WBGene00009621"))
#' #will return a simple list of the matrices for both those genes from all samples
#' getMatrices(methMats,sampleName="182_dSMFv002amp")
#' #will return a simple list of the matrices for all gene from this sample only
#' getMatrices(methMats,regionName="WBGene00009234",sampleName="182_dSMFv002amp")
#' #will return a single item list of a matrix for that gene in that sample.
#' @export
getMatrices<-function(methMats,regionName=c(),sampleName=c()) {
  # extracts a simple list of matrices. It is possible to also subset by regionName or sampleName
  # each matrix will have a name composed of sampleName__regionName
  newMats<-list()
  if(length(sampleName)==0) {
    sampleName<-unique(methMats$sample)
  }
  if(length(regionName)==0) {
    regionName<-unique(methMats$region)
  }
  s<-methMats$sample %in% sampleName
  r<-methMats$region %in% regionName
  idx<-r*s
  for (i in 1:nrow(methMats[idx,])) {
    mat<-readRDS(methMats$filename[i])
    newName<-paste0(methMats[idx,"sample"][i], "__", methMats[idx,"region"][i])
    newMats[newName]<-mat
  }
  return(newMats)
}


#' Merge multiple matrices into a single matrix
#'
#' This function takes a list of matrices with the same number of columns and merges
# them to a single matrix. row names will now contain sampleName__regionName__readName
#'
#' @param matList A table of paths to matrices which have the same columns
#' @return A single methylation matrix
#' @export
rbindMatrixList<-function(matList) {
  naRows<-is.na(matList$filename)
  matList<-matList[!naRows,]
  first<-TRUE
  for (i in 1:nrow(matList)) {
    mat<-readRDS(matList$filename[i])
    row.names(mat)<-paste0(matList$sample[i],"__",matList$region[i],"__",row.names(mat))
    if (first==TRUE) {
      mergedMat<-mat
      first<-FALSE
    } else {
      mergedMat<-rbind(mergedMat,mat)
    }
  #row.names(mergedMat)<-paste0(rep(names(matList),sapply(matList,nrow)),"__",row.names(mergedMat))
  }
  return(mergedMat)
}


#' Convert C position numbering from genomic to relative coordinates
#'
#' @param mat A methylation matrix
#' @param regionGR A genomicRanges object of the region relative to which the new coordinates are caclulated
#' @param invert  A logical variable to indicate if the region should be inverted (e.g. if it is on the negative strand). Default: FALSE
#' @return A methylation matrix in which the column names have been changed from absolute genomic positions to relative
#' positions within the genomicRange regionGR
#' @export
getRelativeCoord<-function(mat,regionGR,invert=F){
  # converts matrix from absolute genome coordinates to
  # relative coordinates within a genomic Range
  pos<-as.numeric(colnames(mat))
  regionStart<-GenomicRanges::start(regionGR)
  regionEnd<-GenomicRanges::end(regionGR)
  if (invert==F) {
    newPos<-pos-regionStart
    colnames(mat)<-as.character(newPos)
  } else {
    newPos<-regionEnd-pos
    colnames(mat)<-newPos
    mat<-mat[,order(as.numeric(colnames(mat))),drop=F]
    colnames(mat)<-as.character(colnames(mat))
  }
  return(mat)
}


#' Change the anchor coordinate
#'
#' @param mat A methylation matrix
#' @param anchorCoord The coordinate which will be set as the 0 position for relative coordinates
#' (default=0)
#' @return A methylation matrix in which the column names have been changed to indicate relative position
#' with reference to the anchor coordinate. e.g. a 500bp matrix centered around the TSS can have its column
#' names changed from 0 to 500 range to -250 to 250 range by setting anchorCoord=250.
#' @export
changeAnchorCoord<-function(mat,anchorCoord=0) {
  # changes the 0 coordinate position (anchorCoord) of
  # a matrix. e.g. sets position 250 to 0 in a 500 region around TSS
  # to get +-250 bp around TSS
  pos<-as.numeric(colnames(mat))
  newPos<-pos-anchorCoord
  colnames(mat)<-as.character(newPos)
  return(mat)
}

#' Get full matrix of all positions in a window even if no C
#'
#' In order to compare different promoters, we need to create a padded
#' matrix with NAs in positions in between Cs.
#'
#' @param matList A table of paths to matrices which have the same columns
#' @param regionType A name for the type of region the matrices desribe
#' @param winSize The size (in bp) of the window containing the matrices
#' @param workDir The path to working directory
#' @return A table of file paths to padded methylation matrices
#' @export
getFullMatrices<-function(matList,regionType,winSize=500, workDir=".") {
  naRows<-is.na(matList$filename)
  matList<-matList[!naRows,]
  makeDirs(workDir,paste0("rds/paddedMats_",regionType))
  matrixLog<-matList[,c("filename","sample","region")]
  matrixLog$filename<-NA
  for (i in 1:nrow(matList)) {
    mat<-readRDS(matList$filename[i])
    # create a matrix with winSize columns and one row per seq read
    Cpos<-colnames(mat)
    withinRange<- -winSize/2<=as.numeric(Cpos) & winSize/2>=as.numeric(Cpos)
    fullMat<-matrix(data=NaN,nrow=dim(mat)[1],ncol=winSize)
    colnames(fullMat)<-c(seq(-winSize/2,-1),seq(1,winSize/2))
    fullMat[,Cpos[withinRange]]<-mat[,withinRange]
    matName<-paste0(workDir,"/rds/paddedMats_",regionType,"/",currentSample,"_",regionGR$ID,".rds")
    saveRDS(fullMat,file=matName)
    matrixLog[i,"filename"]<-matName
  }
  utils::write.csv(matrixLog,paste0(workDir,"/csv/MatrixLog_paddedMats_",regionType,".csv"), quote=F, row.names=F)
  return(matrixLog)
}




#' Convert C position numbering from genomic to relative coordinates for a list of matrices
#'
#' @param matList A table of paths to methylation matrices with names that match the regionGRs object
#' @param regionGRs A genomicRanges object of the regions relative to which the new coordinates are caclulated
#' with a metadata column called "ID" containing names that match the methylation matrices in matList
#' @param anchorCoord  The coordinate which will be set as the 0 position for relative
#' coordinates (default=0)
#' @param workDir path to working directory
#' @return A list of methylation matrices that have been converted from abslute genomic coordinates
#'  to relativepositions within the genomicRanges regionGRs. regionGRs on the negative strand will be flipped to be
#'  in the forward orientation.
#' @export
getRelativeCoordMats<-function(matList, regionGRs, regionType, anchorCoord=0,workDir=".") {
  naRows<-is.na(matList$filename)
  matList<-matList[!naRows,]
  makeDirs(workDir,paste0("/rds/relCoord_",regionType))
  matrixLog<-matList[,c("filename","sample","region")]
  matrixLog$filename<-NA
  matrixLog$reads<-NA
  matrixLog$motifs<-NA
  for (i in 1:nrow(matList)) {
    print(matList[i,c("sample","region")])
    mat<-readRDS(matList$filename[i])
    if(sum(dim(mat)==c(0,0))<1) {
       regionID<-matList$region[i]
       regionGR<-regionGRs[regionGRs$ID==regionID]
       newMat<-getRelativeCoord(mat, regionGR,
                                invert=ifelse(GenomicRanges::strand(regionGR)=="+",F,T))
       newMat<-changeAnchorCoord(mat=newMat,anchorCoord=anchorCoord)
    } else {
      newMat<-mat
    }
    matName<-paste0(workDir,"/rds/relCoord_",regionType,"/",matList$sample[i],"_",regionGR$ID,".rds")
    saveRDS(newMat,file=matName)
    matrixLog[i,"filename"]<-matName
    matrixLog[i,"reads"]<-dim(newMat)[1]
    matrixLog[i,"motifs"]<-dim(newMat)[2]
  }
  if (! dir.exists(paste0(workDir,"/csv/"))) {
    dir.create(paste0(workDir,"/csv/"))
  }
  utils::write.csv(matrixLog,paste0(workDir,"/csv/MatrixLog_relCoord_",regionType,".csv"), quote=F, row.names=F)
  return(matrixLog)
}


#' Calculate methylation frequency at each position for metagene plots
#'
#' @param matList A table of filepaths of methylation matrices with names that match the regionsGRs object
#' @param regionGRs A genomicRanges object of the regions for which aggregate methylation frequency
#' will be calculated. the object must contain a metadata column called "ID" containing names that
#' match the methylation matrices in matList
#' @param minReads  The minimal number of reads a matrix must have in order to be used (default=50)
#' @return A long form  dataframe with four columns: "position" is the C position within the genomic Ranges,
#' "methFreq" is the frequency of methylation at that position, "ID" is the name of the region, "chr"
#' is the chromosome on which that region is present.
#' @export
getMetaMethFreq<-function(matList,regionGRs,minReads=50) {
  naRows<-is.na(matList$filename)
  matList<-matList[!naRows,]
  first=TRUE
  for (i in 1:nrow(matList)) {
    mat<-readRDS(matList$filename[i])
    if (dim(mat)[1]>minReads) {
      vecSummary<-colMeans(mat,na.rm=T)
      df<-data.frame("position"=names(vecSummary),"methFreq"=vecSummary,stringsAsFactors=F)
      df$ID<-matList$region[i]
      df$chr<-as.character(GenomicRanges::seqnames(regionGRs)[match(matList$region[i],regionGRs$ID)])
      if (first==TRUE){
        methFreqDF<-df
        first=FALSE
      } else {
        methFreqDF<-rbind(methFreqDF,df)
      }
    }
  }
  return(methFreqDF)
}





#' Single molecule plot of a methylation matrix
#'
#' @param mat A methylation matrix
#' @param regionName The name of the region the matrix is taken from should match on of the IDs in regionGRs
#' @param regionGRs A genomicRanges object which includes the region which the mat matrix provides the data for.
#' The object must contain a metadata column called "ID" containing names that match the methylation
#' matrices in mat
#' @param featureGRs  A genomicRanges object denoting features to be plotted such as the TSS
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS
#' (default="TSS)
#' @param drawArrow Boolean: should the feature be drawn as an arrow or just a line? (default=TRUE)
#' @param title A title for the plot (default will be the name of the region, the chr and strand on which
#' the region is present)
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @param maxNAfraction Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)
#' @return A ggplot2 plot object
#' @export
plotSingleMolecules<-function(mat,regionName, regionGRs, featureGRs="",
                              myXlab="CpG/GpC position",
                              featureLabel="TSS", drawArrow=TRUE, title=NULL, baseFontSize=12,
                              maxNAfraction=0.2) {
  ### single molecule plot. mat is matrix containing methylation values at different postions
  # (columns) in individual reads (rows). regionName is the ID of the amplicon or genomic
  # region being plotted. regionGRs is a genomicRanges object containing the region being
  # plotted. one of its mcols must have a name "ID" in which the same ID as in regionName
  # appears. featureGRs is genomic ranges object for plotting location of some feature in
  # the region, such as the TSS. myXlab is the X axis label. featureLabel is the label for
  # the type of feature that will be plotted underneath the feature
  tooManyNAs<-rowMeans(is.na(mat))>maxNAfraction
  mat<-mat[!tooManyNAs,]
  if (!is.null(dim(mat)) & any(dim(mat)[1]>10)) {
    regionGR<-regionGRs[match(regionName,regionGRs$ID)]
    if (length(featureGRs)>1) {
      featGR<-featureGRs[match(regionName,featureGRs$ID)]
    }
    na.matrix<-is.na(mat)
    mat[na.matrix]<- -1
    # try to perform heirarchical clustering
    hc <- try(
      stats::hclust(stats::dist(apply(mat,2,as.numeric))),
      silent = TRUE)
    mat[na.matrix]<-NA
    if (class(hc) == "try-error") {
      df<-as.data.frame(mat,stringsAsFactors=F)
      print("hclust failed. Matrix dim: ")
      print(dim(mat))
    } else {
      df<-as.data.frame(mat[hc$order,], stringsAsFactors=F)
    }

    reads<-row.names(df)
    d<-tidyr::gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    #d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      strandInfo<-ifelse(GenomicRanges::strand(featGR)!=GenomicRanges::strand(regionGR),
                         paste0("reg: ",GenomicRanges::strand(regionGR),"ve, ",
                                featureLabel,": ",GenomicRanges::strand(featGR),"ve strand"),
                         paste0(GenomicRanges::strand(regionGR),"ve strand"))
      title=paste0(regionName, ": ",GenomicRanges::seqnames(regionGR)," ",strandInfo)
    }
    p<-ggplot2::ggplot(d,ggplot2::aes(x=position,y=molecules,width=2)) +
      ggplot2::geom_tile(ggplot2::aes(width=3,fill=methylation),alpha=0.8) +
      ggplot2::scale_fill_gradient(low="blue", high="red", na.value="transparent",
                                   breaks=c(0,1), labels=c("protected","accessible"),
                                   limits=c(0,1), name="dSMF\n\n") +
      #ggplot2::scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white", labels=c("protected","accessible"),name="dSMF") +
      ggplot2::theme_light(base_size=baseFontSize) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold",hjust = 0.5),
            legend.position="bottom", legend.key.height = ggplot2::unit(0.2, "cm"),
            legend.key.width=ggplot2::unit(0.5,"cm")) +
      ggplot2::ggtitle(title) +
      ggplot2::xlab(myXlab) +
      ggplot2::ylab("Single molecules") +
      ggplot2::xlim(GenomicRanges::start(regionGR),GenomicRanges::end(regionGR)+10)
    if(length(featureGRs)>0) {
      p<-p+ggplot2::geom_linerange(ggplot2::aes(x=GenomicRanges::start(featGR), y=NULL, ymin=0,
                                       ymax=length(reads)+max(3,0.04*length(reads))),col="black")+
        ggplot2::annotate(geom="text", x=GenomicRanges::start(featGR), y=-max(2,0.03*length(reads)),
                          label=featureLabel,color="black")
      if (drawArrow==TRUE) {
        p<-p+ggplot2::annotate("segment", x = GenomicRanges::start(featGR),
                          xend = GenomicRanges::start(featGR)+
                            20*ifelse(GenomicRanges::strand(featGR)=="-",-1,1),
                          y = length(reads)+max(3,0.04*length(reads)),
                          yend =length(reads)+max(3,0.04*length(reads)),
                          colour = "black", arrow=ggplot2::arrow(length = ggplot2::unit(0.3, "cm")), size=0.7)
      }

    }
  } else {
    p<-NULL
  }
  return(p)
}



#' Single molecule plot of a methylation matrix with average methylation frequency
#'
#' @param mat A methylation matrix
#' @param regionName The name of the region the matrix is taken from should match on of the IDs in regionGRs
#' @param regionGRs A genomicRanges object which includes the region which the mat matrix provides the data for.
#' The object must contain a metadata column called "ID" containing names that match the methylation
#' matrices in mat
#' @param featureGRs  A genomicRanges object denoting features to be plotted such as the TSS
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS
#' (default="TSS)
#' @param drawArrow Boolean: should the feature be drawn as an arrow or just a line? (default=TRUE)
#' @param title A title for the plot (default will be the name of the region, the chr and strand on which
#' the region is present)
#' @param baseFontSize The base font for the plotting theme (default=11 works well for 4x plots per A4 page)
#' @param maxNAfraction Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)
#' @return A ggplot2 plot object
#' @export
plotSingleMoleculesWithAvr<-function(mat, regionName, regionGRs, featureGRs,
                                     myXlab="CpG/GpC position", featureLabel="TSS", drawArrow=TRUE,
                                     title=NULL, baseFontSize=11, maxNAfraction=0.2) {
  # remove reads with more than maxNAfraction positions with NAs
  tooManyNAs<-rowMeans(is.na(mat))>maxNAfraction
  mat<-mat[!tooManyNAs,]
  if(!is.null(dim(mat)) & any(dim(mat)[1]>10)) {
    regionGR<-regionGRs[match(regionName,regionGRs$ID)]
    if (length(featureGRs)>0) {
      featGR<-featureGRs[match(regionName,featureGRs$ID)]
    }
    na.matrix<-is.na(mat)
    mat[na.matrix]<--1
    # try to perform heirarchical clustering
    hc <- try(
      stats::hclust(stats::dist(apply(mat,2,as.numeric))),
      silent = TRUE)
    mat[na.matrix]<-NA
    if (class(hc) == "try-error") {
      df<-as.data.frame(mat,stringsAsFactors=F)
      print("hclust failed. Matrix dim: ")
      print(dim(mat))
    } else {
      df<-as.data.frame(mat[hc$order,],stringsAsFactors=F)
    }

    reads<-row.names(df)
    d<-tidyr::gather(df,key=position,value=methylation)
    d$molecules<-seq_along(reads)
    #d$methylation<-as.character(d$methylation)
    d$position<-as.numeric(d$position)
    if (is.null(title)) {
      strandInfo<-ifelse(GenomicRanges::strand(featGR)!=GenomicRanges::strand(regionGR),
                         paste0("reg: ",GenomicRanges::strand(regionGR),"ve, ",
                                featureLabel,": ",GenomicRanges::strand(featGR),"ve strand"),
                         paste0(GenomicRanges::strand(regionGR),"ve strand"))
      title=paste0(regionName, ": ",GenomicRanges::seqnames(regionGR)," ",strandInfo)
    }
    dAvr<-data.frame(position=as.numeric(colnames(df)),dSMF=1-colSums(df,na.rm=T)/length(reads))
    # average plot
    p1<-ggplot2::ggplot(dAvr,ggplot2::aes(x=position,y=dSMF,group=1)) +
      ggplot2::geom_point()+
      ggplot2::geom_line(size=1,show.legend=F) +
      ggplot2::guides(fill=FALSE, color=FALSE) +
      ggplot2::theme_light(base_size=baseFontSize) +
      ggplot2::ylab("Mean dSMF") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank()) +
      ggplot2::ylim(0,1) +
      ggplot2::xlim(GenomicRanges::start(regionGR),GenomicRanges::end(regionGR)+20)
    if (length(featureGRs)>0) { # plot feature if present
      p1<-p1 + ggplot2::geom_linerange(ggplot2::aes(x=GenomicRanges::start(featGR),
                                                    y=NULL, ymin=0, ymax=1),
                                       col="black",size=0.7)
      if (drawArrow==TRUE) {
        p1<-p1+ggplot2::annotate("segment", x = GenomicRanges::start(featGR),
                              xend = GenomicRanges::start(featGR)+
                                20*ifelse(GenomicRanges::strand(featGR)=="-",-1,1),
                              y = 1, yend = 1, colour = "black", size=0.7,
                              arrow=ggplot2::arrow(length = ggplot2::unit(0.2, "cm")))
      }
    }
    #single molecule plot
    p2<-ggplot2::ggplot(d,ggplot2::aes(x=position,y=molecules,width=2)) +
      ggplot2::geom_tile(ggplot2::aes(width=3,fill=methylation),alpha=0.8) +
      ggplot2::scale_fill_gradient(low="blue", high="red", na.value="transparent",
                                   breaks=c(0,1), labels=c("protected","accessible"),
                                   limits=c(0,1), name="dSMF\n\n") +
      ggplot2::theme_light(base_size=baseFontSize) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_blank(), legend.position="bottom",
                     legend.key.height = ggplot2::unit(0.2, "cm"),
                     legend.key.width=ggplot2::unit(0.5,"cm")) +
      ggplot2::xlab(myXlab) +
      ggplot2::ylab("Single molecules") +
      ggplot2::xlim(GenomicRanges::start(regionGR),GenomicRanges::end(regionGR)+20)
    if (length(featureGRs)>0) { # plot feature if present
      p2<-p2+ggplot2::geom_linerange(ggplot2::aes(x=GenomicRanges::start(featGR), y=NULL, ymin=0,
                                         ymax=length(reads)+max(3,0.04*length(reads))), col="black") +
             ggplot2::annotate(geom="text", x=GenomicRanges::start(featGR),
                               y=-max(2,0.03*length(reads)),
                          label=featureLabel, color="black")

        if (drawArrow==TRUE) {
          p2<-p2+ggplot2::annotate("segment", x = GenomicRanges::start(featGR),
                                 xend = GenomicRanges::start(featGR)+
                                   20*ifelse(GenomicRanges::strand(featGR)=="-",-1,1),
                                 y = length(reads)+max(3,0.04*length(reads)),
                                 yend =length(reads)+max(3,0.04*length(reads)),
                                 colour = "black", size=0.7,
                                 arrow=ggplot2::arrow(length = ggplot2::unit(0.3, "cm")))
        }
    }
    figure<-ggpubr::ggarrange(p1, p2, heights = c(0.5, 2), ncol = 1, nrow = 2, align = "v")
    figure<-ggpubr::annotate_figure(figure, top = ggpubr::text_grob(title, face = "bold"))
  } else {
    figure=NULL
  }
  return(figure)
}


#' Plot a list of list of single molecule matrices
#'
#' This function takes a list (by sample) of a list (by genomic region) of
#' methylation matrices and produces single molecule plots for each amplicon with
#' four samples per page.
#'
#' @param allSampleMats A list (by sample) of lists (by regions) of methylation matrices
#' @param samples A list of samples to plot (same as sample names in allSampleMats)
#' @param regionGRs A genomic regions object with all regions for which matrices should be extracted (same as in allSampleMats). The metadata columns must contain a column called "ID" with a unique ID for each region.
#' @param featureGRs A genomic regions object for features (such as TSS) to be plotted. Feature must be identified with the same ID as the regionGRs
#' @param regionType A collective name for this list of regions (e.g TSS or amplicons). It will be used in naming the output directories
#' @param featureLabel A string with a label for the feature to be added to the plot (default="TSS")
#' @param maxNAfraction Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)
#' @param withAvr Boolean value: should single molecule plots be plotted together with the average profile (default=FALSE)
#' @param includeInFileName String to be included at the end of the plot file name, e.g. experiment name (default="")
#' @param drawArrow Boolean: should the feature be drawn as an arrow or just a line? (default=TRUE)
#' @param workDir Path to working directory
#' @return Plots are written to plots directory
#' @export
plotAllMatrices<-function(allSampleMats, samples, regionGRs, featureGRs, regionType,
                          featureLabel="TSS", maxNAfraction=0.2,withAvr=FALSE,
                          includeInFileName="", drawArrow=TRUE, workDir=".") {
  # convert any factor variables to character
  f <- sapply(allSampleMats, is.factor)
  allSampleMats[f] <- lapply(allSampleMats[f], as.character)
  # remove any regions with no matrix
  naRows<-is.na(allSampleMats$filename)
  allSampleMats<-allSampleMats[!naRows,]
  # get list of all regions in the object
  allAmp2plot<-unique(allSampleMats$region)
  # plot single molecule matrices on their own
  for (i in allAmp2plot) {
    if (withAvr==TRUE) {
      makeDirs(workDir,paste0("plots/singleMoleculePlotsAvr_",regionType))
    } else {
      makeDirs(workDir,paste0("plots/singleMoleculePlots_",regionType))
    }
    plotList=list()
    print(paste0("plotting ", i))
    currentRegion<-allSampleMats[allSampleMats$region==i,]
    for (j in seq_along(currentRegion$sample)) {
      mat<-readRDS(allSampleMats[allSampleMats$region==i &
                                   allSampleMats$sample==currentRegion$sample[j],"filename"])
      maxReads=10000
      if (!is.null(dim(mat))) {
        if (dim(mat)[1]>maxReads) { # if matrix contains more than 10000 reads, do a random subsample
          set.seed(1)
          chooseRows<-sample(1:dim(mat)[1],maxReads)
          mat<-mat[chooseRows,]
        }
        if (withAvr==TRUE) {
          p<-plotSingleMoleculesWithAvr(mat=mat, regionName=i, regionGRs=regionGRs,
                                        featureGRs=featureGRs, myXlab="CpG/GpC position",
                                        featureLabel=featureLabel, drawArrow=drawArrow,
                                        title=currentRegion$sample[j], baseFontSize=11,
                                        maxNAfraction=maxNAfraction)
        } else {
          p<-plotSingleMolecules(mat=mat, regionName=i, regionGRs=regionGRs,
                                 featureGRs=featureGRs, myXlab="CpG/GpC position",
                                 featureLabel=featureLabel, drawArrow=drawArrow,
                                 title=currentRegion$sample[j], baseFontSize=12,
                                 maxNAfraction=maxNAfraction)
        }
        if (!is.null(p)) {
          plotList[[currentRegion$sample[j]]]<-p
        }
      }
    }
    if (length(plotList)>0) {
      numPages=ceiling(length(plotList)/4)
      for (page in 1:numPages) {
        regionGR<-regionGRs[match(i,regionGRs$ID)]
        featGR<-featureGRs[match(i,featureGRs$ID)]
        chr<-GenomicRanges::seqnames(regionGR)
        strandInfo<-ifelse(GenomicRanges::strand(featGR)!=GenomicRanges::strand(regionGR),
                         paste0("reg: ",GenomicRanges::strand(regionGR),"ve, ",
                                featureLabel,": ",GenomicRanges::strand(featGR),
                                "ve strand"),
                         paste0(GenomicRanges::strand(regionGR),"ve strand"))
        title<-paste0(i, ": ",chr," ",strandInfo)
        spacer<-ifelse(length(includeInFileName)>0,"_","")
        toPlot<-c(1:4)+4*(page-1) #get plots for this page
        toPlot<-toPlot[toPlot<=length(plotList)] #make sure only valid plot numbers used
        mp<-gridExtra::marrangeGrob(grobs=plotList[toPlot],nrow=2,ncol=2,top=title)
        if (withAvr==TRUE) {
          ggplot2::ggsave(paste0(workDir,"/plots/singleMoleculePlotsAvr_", regionType,
                                 "/", chr,"_",i,spacer,includeInFileName,"_",page,".png"),
                        plot=mp, device="png", width=20, height=29, units="cm")
        } else {
          ggplot2::ggsave(paste0(workDir,"/plots/singleMoleculePlots_",regionType,
                                 "/",chr,"_",i, spacer, includeInFileName,"_",page,".png"),
                        plot=mp, device="png", width=29, height=20, units="cm")
        }
      }
    }
  }
}


#' Convert Genomic Ranges to relative cooridnates
#'
#' Convert a normal genomic ranges to one where the start and end are relative to some
#' anchor point - either the middle or the start of the genomic ranges (e.g. -250 to 250, or
#' 0 to 500). The original start end and strand are stored in the metadata.
#' @param grs A GenomicRanges object to be converted to relative coordinates
#' @param winSize The size of the window you wish to create
#' @param anchorPoint One of "middle" or "start": the position from which numbering starts
#' @return A GenomicRanges object with relative coordinate numbering
#' @export
convertGRtoRelCoord<-function(grs,winSize,anchorPoint="middle") {
  grsRelCoord<-grs
  GenomicRanges::mcols(grsRelCoord)$gnmStart<-GenomicRanges::start(grs)
  GenomicRanges::mcols(grsRelCoord)$gnmEnd<-GenomicRanges::end(grs)
  GenomicRanges::mcols(grsRelCoord)$gnmStrand<-GenomicRanges::strand(grs)
  grsRelCoord<-GenomicRanges::resize(grsRelCoord,width=winSize,fix="center")
  GenomicRanges::strand(grsRelCoord)<-"*"
  if (anchorPoint=="middle") {
    GenomicRanges::start(grsRelCoord)<- -winSize/2
    GenomicRanges::end(grsRelCoord)<- winSize/2
  } else if (anchorPoint=="start") {
    GenomicRanges::start(grsRelCoord)<- 1
    GenomicRanges::end(grsRelCoord)<- winSize
  } else {
    print("anchorPoint must be one of 'middle' or 'start'")
  }
  return(grsRelCoord)
}


#' get average methylation frequency from all matrices
#'
#' In order to do a metagene plot from matrices, the average methylation frequency from all
#' matrices with more reads than minReads is collected into a data frame
#' @param relCoordMats A table of paths to methylation matrices which have been converted to relative coordinates
#' @param samples A list of samples to plot (same as sample names in relCoordMats)
#' @param regionGRs A genomic regions object with all regions for which matrices should be extracted (same as in relCoordMats). The metadata columns must contain a column called "ID" with a unique ID for each region.
#' @param minReads The minimal number of reads required in a matrix for average frequency to be calculated.
#' @return Data frame with average methylation at relative coordinates extracted from all matrices for all samples
#' @export
getAllSampleMetaMethFreq<-function(relCoordMats,samples,regionGRs,minReads=10) {
  naRows<-is.na(relCoordMats$filename)
  relCoordMats<-relCoordMats[!naRows,]
  first<-TRUE
  for (i in seq_along(samples)) {
    idx<-relCoordMats$sample==samples[i]
    metaMethFreqDF<-getMetaMethFreq(matList=relCoordMats[idx,],
                                    regionGRs=regionGRs, minReads=minReads)
    print(samples[i])
    metaMethFreqDF$sample<-samples[i]
    if(first==TRUE) {
      allSampleMetaMethFreqDF<-metaMethFreqDF
      first<-FALSE
    } else {
      allSampleMetaMethFreqDF<-rbind(allSampleMetaMethFreqDF,metaMethFreqDF)
    }
  }
  # convert position from factor to numeric
  allSampleMetaMethFreqDF$position<-as.numeric(as.character(allSampleMetaMethFreqDF$position))
  return(allSampleMetaMethFreqDF)
}




#' Plot metagene by sample
#'
#' Plots metagene methylation frequency from dataframe extracted from matrices.
#' @param metageneDF ata frame with average methylation at relative coordinates extracted from all matrices for all samples
#' @param maxPoints The maximum number of points to plot per sample. To avoid large files with too much overplotting, the defualt limit is set to 10000. Larger dataframes will be randomly sub-sampled
#' @return A ggplot2 plot object
#' @export
plotDSMFmetageneDF<-function(metageneDF,maxPoints=10000) {
  # subsample if too many points
  if (nrow(metageneDF)>maxPoints) {
    set.seed(1)
    idx<-sample(1:nrow(metageneDF),maxPoints)
  } else {
    idx<-1:nrow(metageneDF)
  }

  p1<-ggplot2::ggplot(metageneDF[idx,],ggplot2::aes(x=position,y=1-methFreq)) +
    ggplot2::theme_light(base_size=16) + ggplot2::ylim(0,1) +
    ggplot2::xlab("Position relative to TSS") + ggplot2::ylab("dSMF (1-%methylation)") +
    ggplot2::geom_linerange(ggplot2::aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
    ggplot2::geom_point(alpha=0.1) +
    ggplot2::geom_smooth(colour="red",fill="red") +
    ggplot2::facet_wrap(~sample)

  p2<-ggplot2::ggplot(metageneDF,ggplot2::aes(x=position,y=1-methFreq,colour=sample)) +
    ggplot2::theme_light(base_size=16) + ggplot2::ylim(0,1) +
    ggplot2::xlab("Position relative to TSS") + ggplot2::ylab("dSMF (1-%methylation)") +
    ggplot2::geom_linerange(ggplot2::aes(x=0, y=NULL, ymin=0, ymax=1),color="steelblue",size=1) +
    ggplot2::geom_smooth(se=FALSE)

  ml <- gridExtra::marrangeGrob(list(p1,p2), nrow=1, ncol=1)
  return(ml)
}

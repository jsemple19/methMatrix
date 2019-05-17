#' Apply a function on a data from gr2 with windows provided by gr1
#'
#' Apply a function to summarise data from multiple metadata columns of GRanges in gr2 that
#' fall within GRanges provided by gr1
#' @param gr1 A GenomicRanges object with windows of interest
#' @param gr2 A GenomicRanges object with data of interest
#' @param applyTo Vector with names of columns in gr2 on which to apply the function
#' @param fun Function to apply (e.g. sum, mean, paste0)
#' @return GRanges from gr1 with summed values of metadata columns from gr2
#' @examples
#' # single gr
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1,3,7), c(5,6,10),names=paste0("win", letters[1:3])), score=4:6)
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 3, 8), c(1, 3, 8),names=paste0("dataID:", letters[1:3])), score=c(10,20,30))
#' applyGRonGR(gr1,gr2,applyTo="score",fun=sum)
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1,3,7), c(5,6,10),names=paste0("win", letters[1:3])), score=4:6, other=1:3)
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 3, 8), c(1, 3, 8),names=paste0("dataID:", letters[1:3])), score=c(10,20,30), other=c(100,NA,300))
#' applyGRonGR(gr1,gr2,c("score","other"),fun=sum,na.rm=T)
#' @export
applyGRonGR<-function(gr1,gr2,applyTo,fun,...) {
  newGR<-IRanges::subsetByOverlaps(gr1,gr2)
  ol<-IRanges::findOverlaps(gr1,gr2)
  grps<-S4Vectors::queryHits(ol)
  dataCols<-tibble::as_tibble(cbind(S4Vectors::DataFrame(grps),
                                    tibble::as_tibble(GenomicRanges::mcols(gr2)[S4Vectors::subjectHits(ol),applyTo])))
  colnames(dataCols)<-c("grps",applyTo)
  newData<-dataCols %>% dplyr::group_by(grps) %>% dplyr::summarise_all(fun,na.rm=T) %>% dplyr::select(applyTo)
  # mark groups where all values are NA
  NAidx<-dataCols %>% dplyr::group_by(grps) %>% dplyr::summarise_all(allNAs) %>% dplyr::select(applyTo)
  newData[as.matrix(NAidx)]<-NA
  # either replace columns or add columns
  if (sum(!(applyTo %in% colnames(GenomicRanges::mcols(newGR))))==0) {
    GenomicRanges::mcols(newGR)[,applyTo]<-newData[,applyTo]
  } else {
    GenomicRanges::mcols(newGR)<-cbind(GenomicRanges::mcols(newGR),newData)
  }
  return(newGR)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bam2methMats.R
\name{combineCGandGCmatrices}
\alias{combineCGandGCmatrices}
\title{Combine CG and GC methylation matrices}
\usage{
combineCGandGCmatrices(matCG, matGC, regionGR, genomeMotifGR)
}
\arguments{
\item{matCG}{Methylation matrix (reads x positons) of C positions within CG motifs}

\item{matGC}{Methylation matrix (reads x positons) of C positions within GC motifs}

\item{regionGR}{GRanges object of region used to make methylation matrices}

\item{genomeMotifGR}{GRanges object with all unique non-overlapping CG/GC/GCGorCGC motifs in genome}
}
\value{
Merged methylation matrix
}
\description{
Methylation matrices are tricky to merge because:
}
\details{
1. Forward and reverse reads will have methylation calls at different positions
even if they are part of the same motif because the C is shifted by 1 between
the two strands. Therefore we "destrand" the reads by averaging all methylation
calls within the same motif (and ingoring NAs). The genomic positions in the final merged
matrix will be on the position of the 1st bp of the motif for CG motifs and the genomic
position of the 2nd bp of the motif for GC motifs.

2. GCG/CGC motifs will have methylation calls from three Cs and it is not clear
whether the middle C should be part of one motif or the other. Therefore we consider
such triplets as a single motif and average the methylation calls from all three Cs.
The genomic position used in the final merged matrix will be that of the middle bp
of the motif.

3. Longer GCGC/CGCG runs have multiple overlapping motifs layered within them. To deal
with that we split these runs into 2-3bp motifs. If the length of the run is an even
number, it is split into doublets and considered as simple CG or GC motifs. If the
length of the run is an odd number, one triplet is created, and the rest is split into
doublets. This creates a set of unique non-overlapping motifs in the genome that are either
CG or GC 2bp motifs, or a GCG/CGC 2bp motif. This unique set is created before hand
with the nanodsmf::findGenomeMotifs function and can be stored in a RDS file. Doublet
and triplet motifs created from the run are treated the same as isolated doublet and triplet
motifs, as described in points 1. and 2. above.
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/Utilities.R
\name{compareChipInRNASeq}
\alias{compareChipInRNASeq}
\title{Compare selected ncRNAs from the chip with RNA-Seq results}
\usage{
compareChipInRNASeq(oligoTableIDs, rnaSeqResultSeqs, type = "local-global",
  ...)
}
\arguments{
\item{oligoTableIDs}{e-numbers of the microarray results that need to be investigated}

\item{...}{further arguments to pairwiseAlignment method}

\item{character}{vector or DNAStringSet of the RNA-Seq contigs/ncRNAs}

\item{character}{vector of mapping type for pairwiseAlignment method (default: local-global means local in the RNA-Seq and global oligo)}
}
\value{
data.frame with the closest match to each sequence
}
\description{
compares sequences from the microarray with sequences of the RNA-Seq run and returns the best matching sequences
}



#' @title Compare sequences to chip oligos
#' @description compares sequences to the oligonucleotide sequence (reverseComplement of oligo)
#' @param character vector of sequences to be analysed
#' @return data.frame with the closest match to each sequence
#' @export
compareRNASeqInChip = function( rnaSeqCandidateSequences ){
  require(CustomMicroarrayPackage)
  require(sqldf)
  require(Biostrings)
  rnaSeqCandidateSequenceL = as.list( as.character(rnaSeqCandidateSequences ))
  db <- dbConnect(SQLite(), dbname=system.file("data/ncrna.chip2.annotation.db",package="CustomMicroarrayPackage"))
  oligoTable = dbGetQuery(db, "SELECT * FROM oligoIDMapping oidm JOIN oligo o ON oidm.oligo_oligoTableID = o.oligoTableID")
  oligoTable$targetSequence = as.character(reverseComplement(DNAStringSet(oligoTable$sequence)))  
  
  resultDFL = lapply(rnaSeqCandidateSequenceL, function(x){
    ali = pairwiseAlignment(pattern=DNAStringSet(oligoTable$targetSequence), subject = DNAString(x))
    maxScoreIdx = which.max( score(ali) )
    curr = oligoTable[maxScoreIdx,]
    currAlign = ali[maxScoreIdx]
    curr$align_pattern = as.character(pattern(currAlign))
    curr$align_subject = as.character(subject(currAlign))
    curr$align_score = as.numeric(score(currAlign))
    return( curr )
  })
  
  resultDF = do.call(rbind, resultDFL)
  return(resultDF)  
}

#' @title Compare query to subject sequences
#' @description This method compares query sequences to subject sequences and returns a data.frane with the closest match as rows of subject sequences
#' @param character vector of querySeqs
#' @param character vector of subjectSeqs
#' @return data.frame with the closest match to each sequence
#' @export
compareQueryWithSubjectSequences = function(querySeqs, subjectSeqs){
  require(Biostrings)
  subjectSeqsL = as.list( as.character(subjectSeqs ))
  
  resultDFL = lapply(subjectSeqsL, function(x){
    ali = pairwiseAlignment(pattern=DNAStringSet(querySeqs), subject = DNAString(x))
    maxScoreIdx = which.max( score(ali) )    
    currAlign = ali[maxScoreIdx]
    curr = data.frame("align_pattern"= as.character(pattern(currAlign)), "align_subject"= as.character(subject(currAlign)), "align_score"= as.numeric(score(currAlign)))
    return( curr )
  })
  
  resultDF = do.call(rbind, resultDFL)
  return(resultDF)  
  
}



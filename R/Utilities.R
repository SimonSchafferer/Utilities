
#' @title Wrapper for draw.pairwise.venn function in VennDiagram
#' @description This function is a wrapper for draw.pairwise.venn.
#' @param ids of sample A
#' @param ids of sample B
#' @param category (see draw.pairwise.venn)
#' @param fill (see draw.pairwise.venn)
#' @param cat.col (see draw.pairwise.venn)
#' @param lty (see draw.pairwise.venn)
#' @param cex (see draw.pairwise.venn)
#' @param cat.cex (see draw.pairwise.venn)
#' @return result object of draw.pairwise.venn
#' @export
createDoubleVenn = function( idsA, idsB, idsC, category = c("First", "Second"), 
                             fill = c("blue", "red"), cat.col=c("blue", "red"),
                             lty = "blank", cex = 2, cat.cex = 2, ... ){
  
  require(VennDiagram)
  # A more complicated diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(idsA),
    area2 = length(idsB),
    cross.area = length(intersect(idsA,idsB)),
    category = category,
    fill = fill,
    lty = lty,
    cex = cex,
    cat.cex = cat.cex,
    cat.col = cat.col
    , ...);
  return(venn.plot)
}

#' @title Wrapper for draw.triple.venn function in VennDiagram
#' @description This function is a wrapper for draw.triple.venn
#' @param ids of sample A
#' @param ids of sample B
#' @param ids of sample C
#' @param category (see draw.triple.venn)
#' @param fill (see draw.triple.venn)
#' @param cat.col (see draw.triple.venn)
#' @param lty (see draw.triple.venn)
#' @param cex (see draw.triple.venn)
#' @param cat.cex (see draw.triple.venn)
#' @return result object of (draw.triple.venn)
#' @export
createTripleVenn = function( idsA, idsB, idsC, category = c("First", "Second", "Third"), 
                             fill = c("blue", "red", "green"), cat.col=c("blue", "red", "green"),
                             lty = "blank", cex = 2, cat.cex = 2, ... ){
  require(VennDiagram)
  # A more complicated diagram
  venn.plot <- draw.triple.venn(
    area1 = length(idsA),
    area2 = length(idsB),
    area3 = length(idsC),
    n12 = length(intersect(idsA,idsB)),
    n23 = length(intersect(idsB,idsC)),
    n13 = length(intersect(idsA,idsC)),
    n123 = length(intersect( intersect(idsA,idsB), idsC) ),
    category = category,
    fill = fill,
    lty = lty,
    cex = cex,
    cat.cex = cat.cex,
    cat.col = cat.col
    , ...);
  return(venn.plot)
}

#' @title Wrapper for draw.quad.venn function in VennDiagram
#' @description This function is a wrapper for draw.quad.venn
#' @param ids of sample A
#' @param ids of sample B
#' @param ids of sample C
#' @param ids of sample D
#' @param category (see draw.quad.venn)
#' @param fill (see draw.quad.venn)
#' @param cat.col (see draw.quad.venn)
#' @param lty (see draw.quad.venn)
#' @param cex (see draw.quad.venn)
#' @param cat.cex (see draw.quad.venn)
#' @return result object of (draw.quad.venn)
#' @export
createQuadVenn = function( idsA, idsB, idsC, idsD, category = c("First", "Second", "Third","Fourth"), 
                           fill = c("orange", "red", "green", "blue"), cat.col=c("orange", "red", "green", "blue"),
                           lty = "blank", cex = 2, cat.cex = 2, ... ){
  
  require(VennDiagram)
  # A more complicated diagram
  venn.plot <- draw.quad.venn(
    area1 = length(idsA),
    area2 = length(idsB),
    area3 = length(idsC),
    area4 = length(idsD),
    n12 = length(intersect(idsA,idsB)),
    n13 = length(intersect(idsA,idsC)),
    n14 = length(intersect(idsA,idsD)),
    n23 = length(intersect(idsB,idsC)),
    n24 = length(intersect(idsB,idsD)),
    n34 = length(intersect(idsC,idsD)),
    n123 = length(intersect( intersect(idsA,idsB), idsC) ),
    n124 = length(intersect( intersect(idsA,idsB), idsD) ),
    n134 = length(intersect( intersect(idsA,idsC), idsD) ),
    n234 = length(intersect( intersect(idsB,idsC), idsD) ),
    n1234 = length(intersect( intersect( intersect(idsA,idsB), idsC), idsD ) ),
    category = category,
    fill = fill,
    lty = lty,
    cex = cex,
    cat.cex = cat.cex,
    cat.col = cat.col
    , ...);
  return(venn.plot)
}

#' @title Calls HyperGTest
#' @description This function is a wrapper for GO-Term HyperG test. It will allow to modify all parameters for creation of GOHyperGParams object
#' @param bg.genes background gene list (ENTREZ IDs)
#' @param geneL genes to test (May be significant genes from Analysis) (ENTREZ IDs)
#' @param conditional (see GOHyperGParams)
#' @param ontology (see GOHyperGParams)
#' @param annotation (see GOHyperGParams)
#' @param testDirection (see GOHyperGParams)
#' @param ... (see GOHyperGParams)
#' @return List containing the summarized result and the result from the hyperGTest
#' @export
performGOTermAnalysis = function( bg.genes, geneL,  conditional=FALSE, pvalueCutoff=0.1, pValCorMethod="hochberg", ontology="BP", annotation="org.Mm.eg.db",testDirection="over",... ){
  require(GOstats)
  message("Performing a Hypergeometric Test on the gene list for GO term association")
  geneL = na.omit(unique(as.numeric(geneL)))
  bg.genes = na.omit(unique(as.numeric(bg.genes)))
  
  params = new("GOHyperGParams", conditional=conditional, geneIds=geneL,
                universeGeneIds=bg.genes, annotation=annotation,
                ontology=ontology, pvalueCutoff=pvalueCutoff, testDirection=testDirection )
  ## performing the analysis
  GO.over = GOstats::hyperGTest(params)
  GO.res = GOstats::summary(GO.over)
  rownames( GO.res ) = GO.res$GOBPID
  
  GO.res$padj = p.adjust(p=GO.res$Pvalue, method=pValCorMethod)
  
  return(list( GORes = GO.res, hyperGTestRes = GO.over))
}

#' @title Get Genes from GO-IDs
#' @description This function returns the genes present in the given dataset for each GO-ID
#' @param geneDF dataframe containing the input genes (ensembl_gene_ids)
#' @param goIDs vector of GO-IDs
#' @param geneDF_ensembleGeneColumnName (The column name in the geneDF that represents the ensembl_gene_ids for merging)
#' @param dataset for biomart (human/mouse etc.) - default: mmusculus_gene_ensembl
#' @param annotDBIdb databse default: org.Mm.eg.db
#' @return List of Data frames (subset of geneDF plus GOID) splitted by GO-ID
#' @export
getEnrichedGenesFromDataset = function( geneDF, goIDs, geneDF_ensembleGeneColumnName = "UID",dataset="mmusculus_gene_ensembl", annotDBIdb = "org.Mm.eg.db"){
  require(biomaRt)
  require(AnnotationDbi)
  require( annotDBIdb, character.only=TRUE)
  if(!is.data.frame(geneDF)){stop("geneDF has to be a data.frame!")}
  
  gGOcatL <- list()
  go2all = eval(parse(text=paste0(annotDBIdb,"::", sub(".db","GO2ALLEGS",annotDBIdb) ) ))
  for(i in goIDs){
    gGOcatL[[i]] = as.vector(unlist(AnnotationDbi::mget(i, go2all)))
  }
  ensembl = useMart("ensembl",dataset=dataset) #uses human ensembl annotations
  
  datL = lapply( gGOcatL, function(gGOcat){
    gene.data <- getBM(attributes=c("ensembl_gene_id",'entrezgene'),
                       filters = 'entrezgene', values = gGOcat, mart = ensembl)
    return( merge( geneDF, gene.data, by.x=geneDF_ensembleGeneColumnName, by.y="ensembl_gene_id" ) )
  } )
  
  return(datL)
}


#' @title Prepares miRNA cluster
#' @description loads clusters of miRNAs according to Wen-Ching Chan et al 2012
#' @param character default human, other possibility: mouse
#' @return list that contains some annotation for each cluster and a named list (names miRNA cluster ID) and miRNA names as vector
#' @export
getMiRNAClusterL = function( species="human"){
  #http://www.sciencedirect.com/science/article/pii/S0888754312001140
  
  #   Table was obtained from the material and methods section of the publication: 
  #     MetaMirClust: Discovery of miRNA cluster patterns using a data-mining approach
  #   Wen-Ching Chan et al 2012
  #   miClust = read.table(file.path(resourcesDir,"miRNA_clusters.csv"), sep="\t", header=TRUE )
  #   idType = miClust[ !miClust$MirClustID_Species %in% c("hg19","mm9") ,c("MirClustID_Species","MirClustType")]
  #   colnames(miClust) = sub("\\..*$","",colnames(miClust))
  #   save( miClust, file="/home/simon/RWorkspace/Utilities/data/miClust.rda" )
  load( system.file( file.path("data","miClust.rda"), package="Utilities") )
  if( species == "human" ){
    species = "hg19"
  } else if(species == "mouse") {
    species = "mm9"
  } else{
    stop("Only mouse and human allowed!")
  }
  
  miClust = miClust[which(miClust$MirClustID_Species == species), ]
  miClust$ID = paste0("Clust",1:dim(miClust)[1])
  miClustL = apply( miClust, 1, function(x){
    return( unlist( strsplit( x["MirGenes"], split=" " ) ) )
  })
  names(miClustL) = miClust$ID
  
  return( list(miClustL, miClust) )  
}


#' @title Download miRNA annotations from mirbase
#' @description downloads the latest mirbase miRNA fasta file and extracts the names and identifier as data.frame
#' @param species optional species identifier to subset the data (mmu for example)
#' @return data.frame containing Name, Id, Name_wo_Species, Species
#' @export
downloadImportMiRNAAnnotation_mirBase = function(speciesFilter=""){
  require(Biostrings)
  speciesColumn = 3
  currTmpFile = tempdir()
  dir.create(currTmpFile)
  download.file("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", file.path(currTmpFile,"mature.fa.gz") )
  system( paste0( "gunzip ", file.path(currTmpFile,"mature.fa.gz") ) )
  matureMiRNAs = readRNAStringSet(file.path(currTmpFile,"mature.fa"))
  if( speciesFilter != "" ){
    mmuMirNames = names(matureMiRNAs)[grep(speciesFilter,names(matureMiRNAs))]  
  } else{
    mmuMirNames = names(matureMiRNAs)
  }
  
  mmuMirNamesL = strsplit(mmuMirNames, " ")
  mmuMirNamesDF = as.data.frame( do.call(rbind, mmuMirNamesL) )
  mmuMirNamesDF$Species = paste0( mmuMirNamesDF[,speciesColumn], "_", mmuMirNamesDF[,(dim(mmuMirNamesDF)[2]-1) ])
  mmuMirNamesDF = mmuMirNamesDF[,-c(speciesColumn,(dim(mmuMirNamesDF)[2]-2))] #now -2 since we added a column (Species) in the last step
  colnames(mmuMirNamesDF) = c("Name","ID","Name_wo_Species","Species")
  return(mmuMirNamesDF)
}

#' @title Download and import MiRNAFamilies from mirbase
#' @description downloads the latest mirbase miRNA family file and imports it into R as list
#' @return list named list containing miRNAfam ID and miRNAfam name as name and a data.frame for each familiy containing miRNA name and ID
#' @export
downloadImportMiRNAFamilies_mirBase = function(){
  currTmpFile = tempdir()
  dir.create(currTmpFile)
  download.file("ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.gz", file.path(currTmpFile,"miFam.dat.gz") )
  system( paste0( "gunzip ", file.path(currTmpFile,"miFam.dat.gz") ) )
  miFamMirBase = read.delim(file.path(currTmpFile,"miFam.dat"), sep="", head=FALSE )
  miFamMirBase = miFamMirBase[which(miFamMirBase$V1 != "//"),]
  miFamMirBase$V2[which(miFamMirBase$V1 == "AC")] = paste0( miFamMirBase$V2[which(miFamMirBase$V1 == "AC")],"_",miFamMirBase$V2[which(miFamMirBase$V1 == "ID")])
  miFamMirBase$V1[ miFamMirBase$V1 == "AC"] = "AC_ID"
  miFamMirBase = miFamMirBase[which( miFamMirBase$V1 != "ID"),]
  miFamMirBase$V4 = miFamMirBase$V1
  miFamMirBase$V4[ miFamMirBase$V4 == "AC_ID"] = paste0( miFamMirBase$V4[miFamMirBase$V4 == "AC_ID"], 1:length(which(miFamMirBase$V4== "AC_ID")) ) 
  miFamMirBase$V5 = miFamMirBase$V4
  for (i in 1:dim(miFamMirBase)[1] ){
    if(miFamMirBase$V1[i]== "AC_ID"){
      curr = miFamMirBase$V4[i]  
    } else if( miFamMirBase$V4[i] == "MI" ){
      miFamMirBase$V5[i] = curr
    }  
  }
  miFamMirBaseL = split(miFamMirBase, miFamMirBase$V5)
  names(miFamMirBaseL) = lapply( miFamMirBaseL, function(x){x[1,2]}  )
  miFamMirBaseL = lapply(miFamMirBaseL, function(x){ x = x[-1,]; x = x[,-c(1,4,5)]; colnames(x) = c("ID","Name"); return(x) })
  
  return(miFamMirBaseL)
}


#' @title Compare selected ncRNAs from the chip with RNA-Seq results
#' @description compares sequences from the microarray with sequences of the RNA-Seq run and returns the best matching sequences
#' @param oligoTableIDs e-numbers of the microarray results that need to be investigated
#' @param character vector or DNAStringSet of the RNA-Seq contigs/ncRNAs
#' @param character vector of mapping type for pairwiseAlignment method (default: local-global means local in the RNA-Seq and global oligo)
#' @param ... further arguments to pairwiseAlignment method
#' @return data.frame with the closest match to each sequence
#' @export
compareChipInRNASeq = function( oligoTableIDs, rnaSeqResultSeqs, type="local-global", ... ){
  require(CustomMicroarrayPackage)
  require(sqldf)
  require(Biostrings)
  db <- dbConnect(SQLite(), dbname=system.file("data/ncrna.chip2.annotation.db",package="CustomMicroarrayPackage"))
  
  oligoTable = dbGetQuery(db, paste0("SELECT * FROM oligoIDMapping oidm JOIN oligo o ON 
                                     oidm.oligo_oligoTableID = o.oligoTableID WHERE o.oligoTableID IN( ",
                                     paste(  sub("e","",oligoTableIDs) , collapse=","),");",sep=""))
  
  oligoTable$targetSequence = as.character(reverseComplement(DNAStringSet(oligoTable$sequence)))  
  
  resultDFL = vector("list", dim(oligoTable)[1] )
  for(i in 1:dim(oligoTable)[1]){
    x = oligoTable[i,]
    ali = pairwiseAlignment(pattern=DNAStringSet(rnaSeqResultSeqs), subject = DNAString(x$targetSequence), type=type, ...)
    maxScoreIdx = which.max( score(ali) )
    curr = x
    currAlign = ali[maxScoreIdx]
    curr$align_pattern = as.character(pattern(currAlign))
    curr$align_subject = as.character(subject(currAlign))
    curr$subjectSeq = as.character(rnaSeqResultSeqs)[maxScoreIdx]
    curr$align_score = as.numeric(score(currAlign))
    curr$subjectID = names(rnaSeqResultSeqs)[maxScoreIdx]
    resultDFL[[i]] = curr
  }
  resultDF = do.call(rbind, resultDFL)
  return(resultDF)  
}

#' @title Compare sequences to chip oligos
#' @description compares sequences to the oligonucleotide sequence (reverseComplement of oligo)
#' @param character vector of sequences to be analysed
#' @param character vector of mapping type for pairwiseAlignment method (default: global-local means global oligo and local in the RNA-Seq)
#' @param ... further arguments to pairwiseAlignment method
#' @return data.frame with the closest match to each sequence
#' @export
compareRNASeqInChip = function( rnaSeqCandidateSequences, type="global-local" ,... ){
  require(CustomMicroarrayPackage)
  require(sqldf)
  require(Biostrings)
  rnaSeqCandidateSequenceL = as.list( as.character(rnaSeqCandidateSequences ))
  db <- dbConnect(SQLite(), dbname=system.file("data/ncrna.chip2.annotation.db",package="CustomMicroarrayPackage"))
  oligoTable = dbGetQuery(db, "SELECT * FROM oligoIDMapping oidm JOIN oligo o ON oidm.oligo_oligoTableID = o.oligoTableID")
  oligoTable$targetSequence = as.character(reverseComplement(DNAStringSet(oligoTable$sequence)))  
  
  resultDFL = lapply(rnaSeqCandidateSequenceL, function(x){
    ali = pairwiseAlignment(pattern=DNAStringSet(oligoTable$targetSequence), subject = DNAString(x), type=type, ...)
    maxScoreIdx = which.max( score(ali) )
    curr = oligoTable[maxScoreIdx,]
    currAlign = ali[maxScoreIdx]
    curr$align_pattern = as.character(pattern(currAlign))
    curr$align_subject = as.character(subject(currAlign))
    curr$subjectSeq = x
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



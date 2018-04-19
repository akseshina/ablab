library(Biostrings)
library(FactoMineR)
library(ggplot2)
library(gridExtra)
library(mclust)
library(cluster)
library(e1071)
library(elasticnet)
library(kernlab)
library(scatterplot3d)


# get and filter information from FASTA file of assembly
assembly_data <- function(path_to_data, max_length=800){
  seqs <- readDNAStringSet(path_to_data, format="fasta", use.names=TRUE)
  
  total_length1 <- sum(as.numeric(seqs@ranges@width))
  seqs <- seqs[as.numeric(seqs@ranges@width) > max_length]
  total_length2  <- sum(as.numeric(seqs@ranges@width))
  cat(paste(round(100*(total_length1-total_length2)/total_length1, 2), 
            "% of total length was removed \n", sep=""))
  
  # find GC content
  a_fr = alphabetFrequency(seqs, baseOnly=TRUE, as.prob=TRUE)
  GC_content = a_fr[,"C"] + a_fr[,"G"]
  
  # find coverage
  seq_names <- strsplit(paste(seqs@ranges@NAMES, '_NA_NA_NA_NA_NA_NA_NA_NA_NA', sep=''), "_", fixed=TRUE)
  coverage <- as.numeric(lapply(seq_names,FUN=function(x){x[6]}))
  
  organism <- unlist(lapply(seq_names, '[[', 9))
  
  # find tetranucleotide frequency
  tetra_nucl_fr <- as.data.frame(oligonucleotideFrequency(seqs, width=4, as.prob=TRUE))
  s <- colnames(tetra_nucl_fr)
  flag <- rep(TRUE, 256)
  for (i in 1:256) {  # join reverse complement columns together
    if (flag[i]) {
      flag[i] <- FALSE
      s1 <- s[i]
      s2 <- toString(reverseComplement(DNAString(s1)))
      num_s2 <- match(s2, s)
      if (flag[num_s2]) {
        flag[num_s2] <- FALSE
        tetra_nucl_fr[,s1] <- tetra_nucl_fr[,s1] + tetra_nucl_fr[,s2]
        tetra_nucl_fr[,s2] <- NULL
      }
    }
  }
  
  # find length of contigs
  seq_length <- as.numeric(seqs@ranges@width)
  
  # turn tetranucleotide frequency into primary components
  cur_pca <- PCA(tetra_nucl_fr, scale.unit=TRUE, ncp=30, graph = F)
  coord_pca<- as.data.frame(cur_pca$ind$coord[])
  
  object <- list(tetra_nucl=tetra_nucl_fr,
                 coord_pca=coord_pca,
                 coverage=coverage, 
                 GC_content=GC_content, 
                 seq_length=seq_length,
                 seq_name=seqs@ranges@NAMES,
                 organism=organism)
  class(object) <- "assembly_data"
  return (object)
}

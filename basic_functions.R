library(Biostrings)
library(FactoMineR)
library(ggplot2)
library(gridExtra)
library(mclust)
library(cluster)
library(e1071)


# get information from FASTA file of assembly
assembly_data <- function(path_to_data, max_length){
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
  seq_names <- strsplit(seqs@ranges@NAMES, "_", fixed=TRUE)
  coverage <- as.numeric(lapply(seq_names,FUN=function(x){x[6]}))
  
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
                 seq_name=seqs@ranges@NAMES)
  class(object) <- "assembly_data"
  return (object)
}


# draw plot which shows weighted coverage density
plot_coverage_density <- function(cur_data){
  sum_length = sum(cur_data$seq_length)
  png(filename="coverage_density.png", width=500, height=450)
  plot(density(cur_data$coverage, weights = cur_data$seq_length / sum_length, bw=0.05),
       xlim=c(0, 120), main = "density(coverage, weights=seq_length/sum_length)")
  dev.off()
}


# draw plot which shows coverage on first two principal components
plot_rainbow_coverage <- function(cur_data) {
  
  coord <- cur_data$coord_pca
  coord$cover <- cur_data$coverage
  
  coord <- coord[order(coord$cover), ]
  assign("pl1", qplot(Dim.1, Dim.2, data=coord, colour = log(cover)) 
         + scale_colour_gradientn(colours = rainbow(7)) 
         + theme(panel.background=element_rect(fill="black")))
  
  coord <- coord[order(coord$cover, decreasing=TRUE), ]
  assign("pl2", qplot(Dim.1, Dim.2, data=coord, colour = log(cover)) 
         + scale_colour_gradientn(colours = rainbow(7)) 
         + theme(panel.background=element_rect(fill="black")))
  
  png(filename="coverage_on_pc.png", width=1000, height=450)
  grid.arrange(pl1, pl2, ncol=2)
  dev.off()
}


# count number of clusters by mclust on GC-content and plot results
gc_content_mclust <- function(cur_data) {
  png(filename="GC_mclust.png", width=500, height=450)
  coord <- cur_data$GC_content
  fit <- Mclust(coord)
  cluster <- fit$classification
  n <- fit$G
  
  colors <- sample(rainbow(n))  # different colors for clusters
  
  cur_hist <- hist(coord, plot=FALSE, breaks=100)
  plot(cur_hist, xlab="GC-content")
  
  tu <- par('usr')
  par(xpd=FALSE)
  temp_df <- as.data.frame(coord)
  temp_df$cl <- cluster
  
  for (i in c(1:n)) {
    tu1 <- min(subset(temp_df, cl==i)$coord)
    tu2 <- max(subset(temp_df, cl==i)$coord)
    clip(tu1, tu2, tu[3], tu[4])
    plot(cur_hist, col=colors[i], add=TRUE)  # add color for cluster
  }
  dev.off()
}


# return column with clusterization results
add_clusterization <- function(cur_data, method, num_of_cl, num_of_pc) {
  
  features <- data.frame(cur_data$coord_pca[,"Dim.1"])  # WORK FOR RITA
  for (i in c(2: num_of_pc)) {
    s <- paste("Dim.", toString(i), sep="")
    features[,s] <- cur_data$coord_pca[,s]
  }
  
  #features[,"coverage"] <- cur_data$coverage
  #features <- scale(features)
  
  if (method == "mclust") {
    fit <- Mclust(features, G=num_of_cl)
    return(fit$classification)
    
  } else if (method == "hclust") {
    d <- dist(features, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward")
    cluster <- cutree(fit, k=num_of_cl) # cut tree into n clusters
    print(cluster)
    return(cluster)
    
  } else if (method == "pam") {
    fit <- pam(features, num_of_cl)
    return(fit$clustering)
  }
}


# draw plot on primary components with clusters colored differently
plot_cluster_solution <- function(cur_data, name="") {
  
  pca_with_clusters <- data.frame(Dim.1=cur_data$coord_pca$Dim.1,
                                  Dim.2=cur_data$coord_pca$Dim.2,
                                  cluster=cur_data$cluster)
  
  qplot(Dim.1, Dim.2, data=pca_with_clusters, colour = cluster) +
    scale_colour_gradientn(colours = rainbow(7)) +
    theme(panel.background=element_rect(fill="black"))
  
  ggsave(filename=paste("clusters_on_pc", name, ".png", sep=""), 
         width=5, height=4.5)
}


# write clusters statistics to file
write_cluster_information <- function(cur_data, name="") {
  
  clusters <- data.frame(cluster=cur_data$cluster, 
                         seq_length=cur_data$seq_length, 
                         GC_content=cur_data$GC_content, 
                         coverage=cur_data$coverage)  
  
  sink(paste("stats_", name, ".txt", sep=""))
  
  for (i in c(1:max(cur_data$cluster))) {
    cat("cluster", i, "\n")
    
    cur_cluster <- subset(clusters, cluster==i)
    
    cat(sum(cur_cluster$seq_length), "base pairs \n")
    
    cat(mean(cur_cluster$GC_content), "average GC-content \n")
    
    cat(mean(cur_cluster$coverage), "average coverage \n")
    cat(median(cur_cluster$coverage), "median coverage \n")
    cat("coverage is from", min(cur_cluster$coverage), "to", max(cur_cluster$coverage),"\n")
    
    cat("\n\n")
  }
  
  sink()
}


# carry out a clusterization analysis with given parameters
do_clusterization <- function(cur_data, method, num_of_cl, num_of_pc) {
  cur_data$cluster <- add_clusterization(cur_data, method, num_of_cl, num_of_pc)
  cur_name <- paste("(", method, "_", num_of_cl, "_", num_of_pc, ")", sep="")
  plot_cluster_solution(cur_data, cur_name)
  write_cluster_information(cur_data, cur_name)
}


make_data_analysis <- function(cur_data) {
  plot_coverage_density(cur_data)
  plot_rainbow_coverage(cur_data)
  gc_content_mclust(cur_data)
}
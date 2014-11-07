library(Biostrings)
library(FactoMineR)
library(ggplot2)
library(gridExtra)
library(mclust)
library(cluster)
library(e1071)


# create data frame with tetranucleotide frequency and other features
read_data <- function(path_to_data){
  seqs <- readDNAStringSet(path_to_data, format="fasta", use.names=TRUE)
  
  # find GC content
  a_fr = alphabetFrequency(seqs, baseOnly=TRUE, as.prob=TRUE)
  GC_content = a_fr[,"C"] + a_fr[,"G"]
  
  # find coverage
  seq_names <- strsplit(seqs@ranges@NAMES, "_", fixed=TRUE)
  coverage <- as.numeric(lapply(seq_names,FUN=function(x){x[6]}))
  
  # find tetranucleotide frequency
  tetra_nucl_fr <- oligonucleotideFrequency(seqs, width=4, as.prob=TRUE)
  tetra_nucl_fr <- as.data.frame(tetra_nucl_fr)
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
  
  dimensions <- tetra_nucl_fr
  dimensions$coverage <- coverage
  dimensions$GC_content <- GC_content
  dimensions$seq_length <- seq_length
  dimensions$seq_name <- seqs@ranges@NAMES
  return(dimensions)
}


# remove short contigs from data
length_filter <- function(cur_data, max_length) {
  total_length1 <- sum(cur_data$seq_length)
  cur_data <- subset(cur_data, seq_length > max_length)
  total_length2  <- sum(cur_data$seq_length)
  print(paste(round(100*(total_length1-total_length2)/total_length1, 2), "% of total length were removed", sep=""))
  return(cur_data)
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
  
  cur_pca <- PCA(cur_data[,1:136], scale.unit=TRUE, ncp=2, graph = F)
  
  coord <- cur_pca$ind$coord[,c(1,2)]
  coord <- as.data.frame(coord)
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
mclust_gc <- function(cur_data) {
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


# turn tetranucleotide frequency into primary components
read_pc <- function(data) {
  my_pca <- PCA(data[,1:136], scale.unit=TRUE, ncp=30, graph = F)
  coord_temp <- my_pca$ind$coord[]
  coord_temp <- as.data.frame(coord_temp)
  
  coord_temp$coverage <- data$coverage
  coord_temp$GC_content <- data$GC_content
  coord_temp$seq_length <- data$seq_length
  coord_temp$name <- data$seq_name
  
  return(coord_temp)
}


# add column with clusterization results
add_clusterization <- function(data, method, num_of_cl, num_of_pc) {
  features <- data.frame(data[,"Dim.1"])  # WORK FOR RITA
  for (i in c(2: num_of_pc)) {
    s <- paste("Dim.", toString(i), sep="")
    features[,s] <- data[,s]
  }
  
  #features[,"coverage"] <- data[,"coverage"] 
  #features <- scale(features)
  
  if (method == "mclust") {
    fit <- Mclust(features, G=num_of_cl)
    return(fit$classification)
    
  } else if (method == "hclust") {
    d <- dist(features, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward")
    cluster <- cutree(fit, k=num_of_cl) # cut tree into n clusters
    return(cluster)
    
  } else if (method == "pam") {
    fit <- pam(features, num_of_cl)
    return(fit$clustering)
  }
}


plot_cluster_solution <- function(cur_data, name) {
  
  assign("pl1", qplot(Dim.1, Dim.2, data=cur_data, colour = cluster) 
         + scale_colour_gradientn(colours = rainbow(7)) 
         + theme(panel.background=element_rect(fill="black")))
  
  png(filename=paste("clusters_on_pc", name, ".png", sep=""), width=500, height=450)
  grid.arrange(pl1, ncol=1)
  dev.off()
}


write_cluster_information <- function(cur_data, name) {
  
  n = max(cur_data$cluster)
  
  sink(paste(name, ".txt", sep=""))
  
  for (i in c(1:n)) {
    cat("cluster", i, "\n")
    
    cur_cluster <- subset(cur_data, cluster==i)
    
    cat(sum(cur_cluster$seq_length), "base pairs \n")
    
    cat(mean(cur_cluster$GC_content), "average GC-content \n")
    
    cat(mean(cur_cluster$coverage), "average coverage \n")
    cat(median(cur_cluster$coverage), "median coverage \n")
    cat("coverage is from", min(cur_cluster$coverage), "to", max(cur_cluster$coverage),"\n")
    
    cat("\n\n")
  }
  
  sink()
}


do_clusterization <- function(cur_data, method, num_of_cl, num_of_pc) {
  cur_data <- read_pc(cur_data)
  cur_data$cluster <- add_clusterization(cur_data, method, num_of_cl, num_of_pc)
  plot_cluster_solution(cur_data, paste("(", method, "_", num_of_cl, "_", num_of_pc, ")", sep=""))
  
  write_cluster_information(cur_data, paste(method, "_", num_of_cl, "_", num_of_pc, sep=""))
}
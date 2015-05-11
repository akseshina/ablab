# return GC content of windows and draw simple histogram 
gc_content <- function(path_to_data, max_length=800, name=""){
  seqs <- readDNAStringSet(path_to_data, format="fasta", use.names=TRUE)
  total_length1 <- sum(as.numeric(seqs@ranges@width))
  seqs <- seqs[as.numeric(seqs@ranges@width) > max_length]
  
  all_windows <- function(curseq) {
    n <- as.integer(length(seq(curseq)) / max_length)
    res <- list(subseq(curseq, 1, max_length))
    if (n==1)
      return(res)
    for (i in c(1:(n-1))) {
      res <- c(res, subseq(curseq, start=(max_length*i+1), width=max_length))
    }
    return(res)
  }
  
  subseqs <- unlist(sapply(seqs,  all_windows))
  
  # find GC content
  GC_content = as.numeric(lapply(subseqs, function(x) {
    alf <- alphabetFrequency(x, baseOnly=TRUE, as.prob=TRUE);
    sum(alf[c("G","C")])
  }))
  
  png(filename=paste(name,"_GC_content.png",sep=""), width=500, height=450)
  hist(GC_content, breaks=100, xlab="GC-content", ylab=paste("number of windows with length ", max_length, sep=""), main=name, col="grey")
  dev.off()
  
  return(GC_content)
}


# count number of clusters by mclust on GC-content and plot results
gc_content_mclust <- function(cur_data, name="") {
  
  png(filename=paste(name, "_GC_mclust.png", sep=""), width=500, height=450)
  coord <- cur_data$GC_windows
  fit <- Mclust(coord)
  cluster <- fit$classification
  n <- fit$G
  
  max_length <- 800
  
  colors <- sample(rainbow(n))  # different colors for clusters
  
  cur_hist <- hist(coord, plot=FALSE, breaks=100)
  plot(cur_hist, xlab="GC-content", ylab=paste("number of windows with length ", max_length, sep=""), main="")
  
  tu <- par('usr')
  par(xpd=FALSE)
  temp_df <- as.data.frame(coord)
  temp_df$cl <- cluster
  
  for (i in unique(cluster)) {
    tu1 <- min(subset(temp_df, cl==i)$coord)
    tu2 <- max(subset(temp_df, cl==i)$coord)
    clip(tu1, tu2, tu[3], tu[4])
    plot(cur_hist, col=colors[i], add=TRUE)  # add color for cluster
  }
  dev.off()
}


# draw plot which shows GC-content on first two principal components
plot_rainbow_GC <- function(cur_data, dir="") {
  
  coord <- cur_data$coord_pca
  coord$GC_content <- cur_data$GC_content
  
  qplot(Dim.1, Dim.2, data=coord, colour = GC_content) +
    scale_colour_gradientn(colours = rainbow(7)) +
    theme(panel.background=element_rect(fill="black"))
  
  ggsave(filename=paste(dir, "GC-content_on_pc.png", sep=""),
         width=5, height=4.5)
}

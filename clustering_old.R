
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

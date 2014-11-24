plot_rainbow_coverage_with_point <- function(cur_data) {
  
  cur_pca <- PCA(new[,1:136], scale.unit=TRUE, ncp=2, graph = F, ind.sup=c(27804))
  
  coord <- temp_pca$ind$coord[,c(1,2)]
  coord <- as.data.frame(coord)
  
  png(filename="contig_on_pc.png", width=600, height=450)
  plot(coord$Dim.1, coord$Dim.2, pch=20)
  points(temp_pca$ind.sup$coord[1], temp_pca$ind.sup$coord[2], col="red", pch=20, cex=2)
  dev.off()
}



my_data$cluster <- add_clusterization(my_data, "hclust", 5, 5)

clusters <- data.frame(cluster=my_data$cluster, 
                       seq_length=my_data$seq_length, 
                       coverage=my_data$coverage,
                       seq_name=my_data$seq_name, 
                       coord_pca=my_data$coord_pca)
cluster_1 <- subset(clusters, cluster==1)
cluster_others <- subset(clusters, cluster!=1)


cluster_1$coverage2 <- (cluster_1$coverage * 200)/(193 - 128 + 1)
sum(cluster_1$coverage2 * cluster_1$seq_length)
cluster_others$coverage2 <- (cluster_others$coverage * 200)/(193 - 128 + 1)
sum(cluster_others$coverage2 * cluster_others$seq_length)
sum(cluster_1$coverage2 * cluster_1$seq_length) / sum(cluster_others$coverage2 * cluster_others$seq_length)

temp_cl1 <- my_cluster$name
temp <- as.matrix(temp)

write(as.matrix(cluster_others$seq_name), file="list_not_cyano.txt")


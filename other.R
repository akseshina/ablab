# draw scree plot about PCA
scree_plot <- function(cur_data, dir="") {
  png(filename=paste(dir, "scree_plot.png", sep = ""), res=72,
      width=33, height=20, units="cm")
  cur_pca <- PCA(cur_data$tetra_nucl, scale.unit=TRUE, ncp=2, graph = F)
  barplot(cur_pca$eig[c(1:30),2], names.arg=c(1:30), 
          main="Scree plot", xlab="number of primary components", ylab="% of variance")
  dev.off()
  
}


# draw all general pictures about dataset
make_data_analysis <- function(cur_data) {
  plot_coverage_density(cur_data)
  plot_rainbow_coverage(cur_data)
  gc_content_mclust(cur_data)
  scree_plot(cur_data)
  plot_rainbow_GC(cur_data)
}
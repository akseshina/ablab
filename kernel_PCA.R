# draw Kernel PCA 2D
plot_kernel_pca <- function(cur_data, name="", kernel_type="rbfdot", params=list()) {
  
  kpc <- kpca(~.,data=cur_data$tetra_nucl, kernel=kernel_type, kpar=params, features=2)
  
  df_to_draw <- data.frame(component_1=pcv(kpc)[,1],
                           component_2=pcv(kpc)[,2],
                           organism=cur_data$organism)
  
  qplot(component_1, component_2, data = df_to_draw, colour=organism)
  
  ggsave(filename=paste("kPCA_", name, ".png", sep=""), 
         width=5, height=4.5)
}


# draw Kernel PCA 2D
plot_kPCA_3D <- function(cur_data, name="", angle=40) {
  
  kpc <- kpca(~.,data=as.data.frame(scale(cur_data$tetra_nucl)), kernel='laplacedot', kpar=list(sigma=1), features=3)
  
  df_to_draw <- data.frame(coordinate_1=pcv(kpc)[,1],
                           coordinate_2=pcv(kpc)[,2],
                           coordinate_3=pcv(kpc)[,3],
                           organism=cur_data$organism)
  
  # "pairs" let to manupulate normaly only this colors vector, so let's make it 
  n <- max(as.numeric(df_to_draw$organism))
  colors_set <- rainbow(n)
  cols <- character(nrow(df_to_draw))  
  for (i in c(1:n)) {
    cols[as.numeric(df_to_draw$organism) == i] <- colors_set[i]
    
  }  
  
  png(filename=paste("kPCA_3D_", name, ".png", sep=""), 
      width=600, height=450)
  scatterplot3d(df_to_draw$coordinate_1, df_to_draw$coordinate_2, df_to_draw$coordinate_3,
                pch=20, color=cols,
                xlab="PC 1", ylab="PC 2", zlab="PC 3", angle=angle)
  dev.off()
}


# draw Kernel PCA scatterplot matrices (3D)
plot_kPCA_pairs <- function(cur_data, name="") {
  
  kpc <- kpca(~.,data=as.data.frame(scale(cur_data$tetra_nucl)), kernel='laplacedot', kpar=list(sigma=1), features=3)
  
  df_to_draw <- data.frame(coordinate_1=pcv(kpc)[,1],
                           coordinate_2=pcv(kpc)[,2],
                           coordinate_3=pcv(kpc)[,3],
                           organism=cur_data$organism)
  
  # "pairs" let to manupulate normaly only this colors vector, so let's make it 
  n <- max(as.numeric(df_to_draw$organism))
  colors_set <- rainbow(n)
  cols <- character(nrow(df_to_draw))  
  for (i in c(1:n)) {
    cols[as.numeric(df_to_draw$organism) == i] <- colors_set[i]
    
  }
  
  png(filename=paste("kPCA_pairs_", name, ".png", sep=""), 
      width=600, height=450)
  pairs(~coordinate_1+coordinate_2+coordinate_3, data=df_to_draw,
        pch = 20, col=cols)
  #par(xpd=TRUE)
  #legend(0, 1, as.vector(unique(df_to_draw$organism)),  
  #       fill=colors_set)
  dev.off()
}
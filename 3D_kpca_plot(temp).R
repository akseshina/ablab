kpc <- kpca(~.,data=my_data$my_metagenome$tetra_nucl, kernel='laplacedot', kpar=list(sigma=1), features=3)





plot3d(df_to_draw$component_1, df_to_draw$component_2, df_to_draw$component_3, col=as.numeric(factor(df_to_draw$organism)), pch=20, size=10, xlab='PC 1', ylab='PC 2', zlab='PC 3')


plot3d(my_data$my_metagenome$coord_pca$Dim.1, my_data$my_metagenome$coord_pca$Dim.2, my_data$my_metagenome$coord_pca$Dim.3, col=as.numeric(factor(df_to_draw$organism)), pch=20, size=10)



#____________________
plot_kPCA_3D <- function(cur_data, name="") {
  
  kpc <- kpca(~.,data=scale(cur_data$tetra_nucl), kernel='laplacedot', kpar=list(sigma=1), features=3)
  
  df_to_draw <- data.frame(component_1=pcv(kpc)[,1],
                           component_2=pcv(kpc)[,2],
                           component_3=pcv(kpc)[,3],
                           organism=my_data$my_metagenome$organism)
  
  png(filename=paste("kPCA_", name, ".png", sep=""), 
      width=600, height=450)
  scatterplot3d(df_to_draw$coordinate_1, df_to_draw$coordinate_2, df_to_draw$coordinate_3,
                pch=20, color=as.numeric(df_to_draw$organism))
  dev.off()
}
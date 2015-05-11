


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


---------------------------------------------------------------------
  
  my_spca <- list() 
for (i in c(1:5)) {
  my_spca[[i]] <- spca(my_data[[i]]$tetra_nucl,K=1,type="predictor",
                       sparse="varnum",trace=TRUE,para=c(20))
}
names(my_spca) <- names(my_data)


for (i in c(1:5)) {
  my_spca[[i]]$non_zero_var <- names(which(my_spca[[i]]$loadings!=0,arr.ind = T)[,1])
}

non_zero_var <- list()
for (i in c(1:5)) {
  non_zero_var[[i]] <- names(which(my_spca[[i]]$loadings!=0,arr.ind = T)[,1])
}

Reduce(intersect, non_zero_var)

_________________________________________________________________________

x <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2)
y <- c()
for (i in x) {
  my_spca <- spca(my_data$utex$tetra_nucl,K=4,type="predictor", 
                  sparse="penalty",trace=TRUE,para=c(i, 0.01, 0.01, 0.01))
  y <- c(y, my_spca$pev[1])
}
my_data$utex$spca_penalties$x <- x
my_data$utex$spca_penalties$l1 <- y

cur_pca <- PCA(my_data$bastimolide$tetra_nucl, scale.unit=TRUE, ncp=30, graph = F)
my_data$bastimolide$PEV <- cur_pca$eig[[2]]

plot_spca_penalties <- function(cur_data, save_flag=F) {
  ggplot(as.data.frame(cur_data$spca_penalties), aes(x = x, y=l1_2)) + 
    geom_point() + geom_line() + 
    geom_hline(yintercept=cur_data$PEV[1]/100) + 
    ylab(label="PEV for PC 1") + 
    xlab("lambda_1") 
  #if (save_flag==T) {
  #  ggsave(filename=paste("l1.png", sep=""), 
  #         width=5, height=4.5)  
  #}
}

my_data$bastimolide$spca_penalties$l1_2 <- c(0.60905297, 0.51026312, 0.46072545, 0.15225548, 0.06465068, 0.17959345, 0.18049806, 0.17855573)

______________________________________________
# Multidimensional scaling
plot_MDS <- function(cur_data, name="", dist_type="euclidean") {
  
  d <- dist(cur_data$tetra_nucl, method=dist_type)
  fit <- cmdscale(d, k=2)
  
  df_to_draw <- data.frame(coordinate_1=fit[,1],
                           coordinate_2=fit[,2],
                           organism=cur_data$organism)
  
  qplot(coordinate_1, coordinate_2, data = df_to_draw, colour=organism)
  
  ggsave(filename=paste("MDS_", name, ".png", sep=""), 
         width=5, height=4.5)
}

# Multidimensional scaling in 3D
plot_MDS_3D <- function(cur_data, name="", dist_type="euclidean") {
  
  d <- dist(cur_data$tetra_nucl, method=dist_type)
  fit <- cmdscale(d, k=3)
  
  df_to_draw <- data.frame(coordinate_1=fit[,1],
                           coordinate_2=fit[,2],
                           coordinate_3=fit[,3],          
                           organism=cur_data$organism)
  
  png(filename=paste("MDS_", name, ".png", sep=""), 
      width=600, height=450)
  scatterplot3d(df_to_draw$coordinate_1, df_to_draw$coordinate_2, df_to_draw$coordinate_3,
                pch=20, color=as.numeric(df_to_draw$organism))
  dev.off()
}

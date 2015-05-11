# draw plot which shows weighted coverage density
plot_coverage_density <- function(cur_data, name="", max_x=120){
  sum_length = sum(cur_data$seq_length)
  png(filename=paste(name, "_coverage_density.png", sep=""), width=500, height=450)
  plot(density(cur_data$coverage, weights = cur_data$seq_length / sum_length, bw=0.05),
       xlim=c(0, max_x), main = name)
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
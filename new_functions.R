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
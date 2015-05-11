x=c(120,140,40,100,120,150,120)
n=names(my_data$r)

png(filename=paste("", "_coverage_density.png", sep=""), width=775, height=1000)
par(mfrow=c(4,2))
for (i in c(1:length(n))) {  plot_coverage_density(my_data$r[[n[i]]], name=n[i], max_x=x[i])}
dev.off()

#for (ns in names(my_data$s)) {  plot_coverage_density(my_data$s[[ns]], name=ns, max_x=40)}

#___________________________________


for (ns in names(my_data$s)) { gc_content_mclust(my_data$s[[ns]], name=ns)}
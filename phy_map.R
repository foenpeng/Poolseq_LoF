setwd("/Users/pengfoen/OneDrive - University of Connecticut/poolseq_LoF")

library(data.table)

site_pos<-fread("site_info.csv")

Fst_matrix<-read.csv("Fst_pairwise_pop_WithIndel.csv")
Fst_matrix<-as.matrix(Fst_matrix[,2:9],rownames.value=Fst_matrix[,1])

library(ape)
tree<-nj(Fst_matrix)
nj_tree<-root(tree,"SAY")
colors2 <- rep("black", Nedge(nj_tree))
colors2[which(nj_tree$edge[,2] %in% 1:8)] <- topo.colors(8)
is_tip <- nj_tree$edge[,2] <= length(nj_tree$tip.label)
ordered_tips <- nj_tree$edge[is_tip, 2]

png("./Figures/fig_phylogeny.png",res=300, width=1500,height = 2000)
plot(nj_tree, edge.color=colors2,edge.width=1.5, font=0.8, label.offset = 0.01)
tiplabels(pch = 21, bg=topo.colors(8)[match(c(1:8),ordered_tips)], cex = 1,adj=0.5)
dev.off()

png("./Figures/fig_map.png",res=300, width=1000,height = 1000)
long <- c(-126.2,-125)
lat <- c(48.8,50.4)
center <- c(mean(lat), mean(long))
zoom <- 8
terrmap2 <- GetMap(center=center, zoom=zoom, maptype= "terrain", destfile = "terrain2.png")
setorder(site_pos,site_ID)
site_pos<-site_pos[ordered_tips]
PlotOnStaticMap(terrmap2, lat = site_pos$lat, lon = site_pos$long, pch = 21, cex = 0.8, col="black",bg=topo.colors(8))
PlotOnStaticMap(terrmap2, lat = site_pos$lat, lon = site_pos$long+0.12, FUN= text,labels = site_pos$site_ID,cex=0.4, add=TRUE) 
dev.off()
dev.off()

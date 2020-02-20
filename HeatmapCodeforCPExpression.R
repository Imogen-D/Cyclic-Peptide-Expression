#script to heatmap rank of cylic peptides (and TPM and log(TPM))
#load packages
library(stringr)
library(tidyr)
library(XML)
library("rentrez")
library(RColorBrewer)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(rgl)
library(caTools)
library(magick)
library("viridis") 


#already made: all_read_trees, all_twentyfour_tips, strain_tree, species_tree in EpiAllInclCpur.R
#using allstrainframe from RNAquant.R script

unroot_strain_tree <- drop.tip(strain_rooted_tree, "Cpur")
CPtree <- keep.tip(unroot_strain_tree, c("C2857", "Eam702", "Eel732", "Efe2368", "Ety8"))
x <- plot(CPtree)
ggstraintree <- ggtree(CPtree)
joined <- ggstraintree %<+% allstrainframe
y <- joined +
  geom_tiplab(aes()) +
  theme(legend.position = "right") + # legend.title = element_text()) +
  scale_color_viridis()


allstrainframe <- rbind(Eammerged, Eelmerged, Fl1merged, F2368merged, E8merged)
names(allstrainframe)[6] <- "Proportion"


framespread <- spread(allstrainframe, V1, Proportion)
framespread <- framespread[,6:20]


compress <- function(x) c(na.omit(x), NA)[1]
compressed <- aggregate(framespread[2:15], framespread[1], compress)
rownames(compressed) <- compressed$V3
compressed$V3 <- NULL

gheatmap(y, compressed, offset = 4, colnames_angle = 90,low = "grey90", high = "slateblue4", colnames_offset_y = -0.5, width = 1.5) + theme(legend.title = element_text()) + labs(fill = "Distribution")  

#for TPM
allstrainframe <- rbind(Eammerged, Eelmerged, Fl1merged, F2368merged, E8merged)

framespread <- spread(allstrainframe, V1, TPM)
framespread <- framespread[,6:20]
compress <- function(x) c(na.omit(x), NA)[1]
compressed <- aggregate(framespread[2:15], framespread[1], compress)
rownames(compressed) <- compressed$V3
compressed$V3 <- NULL
gheatmap(y, compressed, offset = 4, colnames_angle = 90,low = "grey90", high = "slateblue4", colnames_offset_y = -0.5, width = 1.5) + theme(legend.title = element_text()) + labs(fill = "Distribution")  

#for log(TPM)

allstrainframe <- rbind(Eammerged, Eelmerged, Fl1merged, F2368merged, E8merged)

allstrainframe$logTPM <- log(allstrainframe$TPM)
framespread <- spread(allstrainframe, V1, logTPM)
framespread <- framespread[,7:21]
compress <- function(x) c(na.omit(x), NA)[1]
compressed <- aggregate(framespread[2:15], framespread[1], compress)
rownames(compressed) <- compressed$V3
compressed$V3 <- NULL
gheatmap(y, compressed, offset = 4, colnames_angle = 90,low = "grey90", high = "slateblue4", colnames_offset_y = -0.5, width = 1.5) + theme(legend.title = element_text()) + labs(fill = "Distribution")  

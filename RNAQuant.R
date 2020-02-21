library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)


Eamquant <- read.delim("<Eamquant.sf>", stringsAsFactors=FALSE)
Eelquant <- read.delim("~/CoxExtension/RNA/Eelquant.sf", stringsAsFactors=FALSE)
Fest2368quant <- read.delim("~/CoxExtension/RNA/Fest2368quant.sf", stringsAsFactors=FALSE)
FestFl1quant <- read.delim("~/CoxExtension/RNA/FestFl1quant.sf", stringsAsFactors=FALSE)
TyphE8quant <- read.delim("~/CoxExtension/RNA/TyphE8quant.sf", stringsAsFactors=FALSE)

ggplot(Eelquant, aes(TPM)) + geom_density() +scale_x_log10()
ggplot(Eamquant, aes(TPM)) + geom_density() +scale_x_log10()
ggplot(TyphE8quant, aes(TPM)) + geom_density() +scale_x_log10()
ggplot(Fest2368quant, aes(TPM)) + geom_density() +scale_x_log10()
ggplot(FestFl1quant, aes(TPM)) + geom_density() +scale_x_log10()


Eelquant <- Eelquant %>% arrange(desc(TPM))
Eamquant <- Eamquant %>% arrange(desc(TPM))
TyphE8quant <- TyphE8quant %>% arrange(desc(TPM))
Fest2368quant <- Fest2368quant %>% arrange(desc(TPM))
FestFl1quant <- FestFl1quant %>% arrange(desc(TPM))

Eamleng <- nrow(Eamquant)
Eamrows <- as.integer(attr(Eamquant, "row.names"))
Eamquant$rank <- Eamrows/Eamleng

Eeleng <- nrow(Eelquant)
Eelrows <- as.integer(attr(Eelquant, "row.names"))
Eelquant$rank <- Eelrows/Eeleng

Fl1leng <- nrow(FestFl1quant)
Fl1rows <- as.integer(attr(FestFl1quant, "row.names"))
FestFl1quant$rank <- Fl1rows/Fl1leng

F2368leng <- nrow(Fest2368quant)
F2368rows <- as.integer(attr(Fest2368quant, "row.names"))
Fest2368quant$rank <- F2368rows/F2368leng

E8leng <- nrow(TyphE8quant)
E8rows <- as.integer(attr(TyphE8quant, "row.names"))
TyphE8quant$rank <- E8rows/E8leng

subsetCP <- read.table("~/CoxExtension/subsetCP.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE)


EamsubsetCP <- subsetCP[c(which(subsetCP$V2 %in% Eamquant$Name)),]
EelsubsetCP <- subsetCP[c(which(subsetCP$V2 %in% Eelquant$Name)),]
E8subsetCP <- subsetCP[c(which(subsetCP$V2 %in% TyphE8quant$Name)),]
Fl1subsetCP <- subsetCP[c(which(subsetCP$V2 %in% FestFl1quant$Name)),]
F2368subsetCP <- subsetCP[c(which(subsetCP$V2 %in% Fest2368quant$Name)),]


Eammerged <- merge(Eamquant, EamsubsetCP, by.x = "Name", by.y = "V2")
Eelmerged <- merge(Eelquant, EelsubsetCP, by.x = "Name", by.y = "V2")
Fl1merged <- merge(FestFl1quant, Fl1subsetCP, by.x = "Name", by.y = "V2")
F2368merged <- merge(Fest2368quant, F2368subsetCP, by.x = "Name", by.y = "V2")
E8merged <- merge(TyphE8quant, E8subsetCP, by.x = "Name", by.y = "V2")

allstrainframe <- rbind(Eammerged, Eelmerged, Fl1merged, F2368merged, E8merged)

names(allstrainframe)[6] <- "Proportion"

plot <- ggplot(allstrainframe, aes(V1, Proportion, fill=V3))
plot <- plot + geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single"))
plot +  facet_grid(~V1, scales = "free_x", space = "free_x", switch = "x") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())


###all genes

ortho_long <- read.delim("~/CoxExtension/ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)

Eamsubset <- ortho_long[c(which(ortho_long$V2 %in% Eamquant$Name)),]
Eelsubset <- ortho_long[c(which(ortho_long$V2 %in% Eelquant$Name)),]
E8subset <- ortho_long[c(which(ortho_long$V2 %in% TyphE8quant$Name)),]
Fl1subset <- ortho_long[c(which(ortho_long$V2 %in% FestFl1quant$Name)),]
F2368subset <- ortho_long[c(which(ortho_long$V2 %in% Fest2368quant$Name)),]


Eamallmerged <- merge(Eamquant, Eamsubset, by.x = "Name", by.y = "V2")
Eelallmerged <- merge(Eelquant, Eelsubset, by.x = "Name", by.y = "V2")
Fl1allmerged <- merge(FestFl1quant, Fl1subset, by.x = "Name", by.y = "V2")
F2368allmerged <- merge(Fest2368quant, F2368subset, by.x = "Name", by.y = "V2")
E8allmerged <- merge(TyphE8quant, E8subset, by.x = "Name", by.y = "V2")


fullstrainframe <- rbind(Eamallmerged, Eelallmerged, Fl1allmerged, F2368allmerged, E8allmerged)
names(fullstrainframe)[6] <- "Proportion"


row.names(Eamallmerged) <- Eamallmerged$V1 #og_0719 and og_3334 repeated #easier to do it with this function as throws up warning 
which(Eamallmerged$V1 == "og_0719") #2286, 5362
which(Eamallmerged$V1 == "og_3334") #1195, 1750
Eamallmerged <- Eamallmerged[-c(1195, 2286, 5362, 1750)]

row.names(Eelallmerged) <- Eelallmerged$V1 #og_0724 and og_3334 repeated #easier to do it with this function as throws up warning
which(Eelallmerged$V1 == "og_0724") #164, 3617
which(Eelallmerged$V1 == "og_3334") #1195, 1750
Eelallmerged <- Eelallmerged[-c(911, 164, 3617, 1743)]

row.names(E8allmerged) <- E8allmerged$V1 #og_3334 repeated #easier to do it with this function as throws up warning
which(E8allmerged$V1 == "og_3334") #1905, 2918
E8allmerged <- E8allmerged[-c(1905, 2918)]

row.names(F2368allmerged) <- F2368allmerged$V1 #og_3334 repeated #easier to do it with this function as throws up warning
which(F2368allmerged$V1 == "og_3334") #4259, 4530
F2368allmerged <- F2368allmerged[-c(4259, 4530)]

row.names(Fl1allmerged) <- Fl1allmerged$V1 #og_3334 repeated #easier to do it with this function as throws up warning
which(Fl1allmerged$V1 == "og_3334") #1280, 4304
Fl1allmerged <- Fl1allmerged[-c(1280, 4304)]

fullstrainframeedited <- rbind(Eamallmerged, Eelallmerged, Fl1allmerged, F2368allmerged, E8allmerged)
fullstrainframeedited <- fullstrainframeedited[,-c(1,2,3,5)]

mydata <- melt(fullstrainframeedited, TPM=c("TPM","rank"))
write.table(mydata, "TMPRankallorthologs.tsv", quote=FALSE, sep='\t', row.names = FALSE)
nrow(mydata)


#making PCAs
TPM <- subset(mydata, variable == "TPM")
TPM_wider <- pivot_wider(TPM, names_from="V3", values_from="value", values_fn = list(value = mean))


TPM_wider <- TPM_wider[complete.cases(TPM_wider),]
row.names(TPM_wider) <- TPM_wider$V1

PCA = prcomp(TPM_wider[3:7], scale = TRUE, center = TRUE)
PCA_df = as.data.frame(PCA$rotation)
PCA_df$strain <- row.names(PCA_df)
ggplot(PCA_df, aes(PC1, PC2, label = strain)) + geom_text() + coord_equal()

PCA_values <- PCA$x
row.names(PCA_values) <- TPM_wider$V1

PCA_values <- PCA_values[,c(1,2)]

pythagorean <- function(a, b){
hypotenuse <- sqrt(a^2 + b^2)
return(hypotenuse)
}
hypo <- pythagorean(PCA_values[,1], PCA_values[,2])
hypo <- as.data.frame(hypo)
distance <- hypo$hypo
PCA_values <- as.data.frame(PCA_values)
PCA_values$distance <- distance

morethan10 <- which(PCA_values$distance >= 10)
subsetPCA_values <- PCA_values[c(morethan10),] #56

dist_orthologs <- row.names(subsetPCA_values)

distrows <- lapply(dist_orthologs, function(x) which(ortho_long$V1 == (x)))
rownums <- unlist(ortho_long_rows)
subsetorthos <- ortho_long[(rownums),]
write.table(subsetorthos, file = "distancegenesmorethan10.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Clustering project
# Author: Eduardo Arrieta
# LCG UNAM fourth semester
# Bioinformatics
# Clustering module

# Loading the libraries
library(cluster)
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))

# Retrieve the BLAST results
blastFile = 'C:/Users/Eduardo/Documents/IntensivoBioinfo/Clustering/ABCvABC_blast.out'
blastOut <- read.table(blastFile, sep = '\t')
# Taging the columns
colnames(blastOut) <- c('query_acc', 'subject_acc', 'identity', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score')

# Creat a 100*100 matrix to save the bitscores
bitscores_matrix <- matrix(data = rep(0, 100*100), nrow = 100, ncol = 100)
colnames(bitscores_matrix) <- unique(blastOut$query_acc)
rownames(bitscores_matrix) <- unique(blastOut$query_acc)

# from the BLAST results, get the only the bitscores for the specific protein pair
for (result in 1:nrow(blastOut)) {
  bitscores_matrix[blastOut[result, 'query_acc'], blastOut[result, 'subject_acc']] <- blastOut[result, 'bit_score']
}

# Get the dissimilarity for every record
disimilitud <- function(scores){
  d = 1 - scores/max(scores)
  return(d)
}
dis_matrix <- apply(bitscores_matrix, 2, disimilitud)

#build and save the dendograms for every method seen in class
csin <- hclust(dist(dis_matrix, method = "euclidean"), method = "single")
write.tree(phy = as.phylo(csin), file = "csin.tree")
coef_csin <- coef(csin)

cave <- hclust(dist(dis_matrix, method = "euclidean"), method = "average")
write.tree(phy = as.phylo(cave), file = "cave.tree")
coef_cave <- coef(cave)

ccom <- hclust(dist(dis_matrix, method = "euclidean"), method = "complete")
write.tree(phy = as.phylo(ccom), file = "ccom.tree")
coef_ccom <- coef(ccom)

cwar <- hclust(dist(dis_matrix, method = "euclidean"), method = "ward.D2")
write.tree(phy = as.phylo(cwar), file = "cwar.tree")
coef_cwar <- coef(cwar)

# Get the agglomerative coefficient for every tree
print(coef_csin)
print(coef_cave)
print(coef_ccom)
print(coef_cwar)
# [1] 0.5845887
# [1] 0.6850752
# [1] 0.7503544
# [1] 0.9306215
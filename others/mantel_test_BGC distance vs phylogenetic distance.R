#Calculating distance matrix from tree file for 1395 isolates
library(ape)
library(readxl)
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyverse)
tree<-read.tree("gtdbtk.bac120.rooted_new.tree")
DistMatrix<-cophenetic(tree)
#Calculate distance matrix based on BGC profile
data<-read_excel("/Users/arijitmukherjee/Downloads/Mantel_isolates.xlsx",col_names = T,skip = 0,sheet = "Sheet1")
head(data)
data$ANI_99.95<-as.factor(data$ANI_99.95)
data$ANI_99<-as.factor(data$ANI_99)
data$ANI_98<-as.factor(data$ANI_98)

data_98<-data%>%filter(ANI_98=="Y")
dim(data_98)
names(data_98)
data_98<-data_98[,-c(2,3,4)]
data_98$sum<-rowSums(data_98[,-1])
sum(data_98$sum==0)#there are 5 genomes without any BGC detected
#find those genome IDs
data_98_filt<-data_98%>%filter(sum!=0)
dim(data_98_filt)
#this is the dataset for no BGC filtered genomes
length(tree$tip.label)#1390 tips
dim(data_98_filt)

isolates_IDs<-tree$tip.label
length(isolates_IDs)
data_98_filt$Genome_ID
#filter the data based on the ids from tip labels of the tree file
data_98_filt<-data%>%filter(Genome_ID%in%isolates_IDs)
dim(data_98_filt)
data_98_filt$Genome_ID
#ordering assemblies based on tip labels
ordered_ids<-tree$tip.label
length(ordered_ids)
#names(ordered_data)
ordered_data<-left_join(data.frame(Genome_ID=ordered_ids),data_98_filt,by="Genome_ID")
dim(ordered_data)

names(ordered_data)
ordered_data_sub<-ordered_data[,5:15]
ordered_data_sub$Assembly<-ordered_data$Genome_ID
rownames(ordered_data_sub)<-ordered_data_sub$Assembly
names(ordered_data_sub)
ordered_data_sub<-ordered_data_sub[,-12]
sum(is.na(ordered_data_sub))
#There are missing va
#perfect 
rownames(ordered_data_sub)
#Now calculate the distance matrix based on the BGC profile on ordered dataframe
dist_isolates<-vegdist(as.matrix(ordered_data_sub),method = "bray")
dim(dist_isolates)
dist_isolates
library(ade4)
mentel_isolates<-mantel.rtest(as.dist(dist_isolates),as.dist(DistMatrix),nrepet = 1000)
saveRDS(mentel_isolates,"mantel_isolates_1390.rds")
#Mantel 1000 premutations p-value: 0.0009 mantel's r=0.31
mantel<-readRDS("mantel_isolates.rds")
mentel_isolates

#Now convert the data matrix into a dataframe for individual pairwise distance matric
library(reshape2)
#check if the rownames of the distmatrices do match, then merge them 
#mantel<-readRDS("mantel_isolates_1390.rds")

BGC_dist<-melt(as.matrix(dist_isolates))
phylo_dist<-melt(as.matrix(DistMatrix))
class(BGC_dist)
class(phylo_dist)
head(BGC_dist)


BGC_dist$phylo_distance<-phylo_dist$phylogenetic_distance
head(BGC_dist)


colnames(BGC_dist)<-c("genome_1","genome_2","BGC_distance")
colnames(phylo_dist)<-c("genome_1","genome_2","phylogenetic_distance")

write.table(BGC_dist,"BGC_dist.tsv",sep = "\t")
write.table(phylo_dist,"Phylogenetic_distance.tsv",sep = "\t")

#can't write this large file for generating figure
library(WriteXLS)
WriteXLS("BGC_dist1",ExcelFileName = "BGC_dist1.xlsx",row.names = F,col.names = T)
WriteXLS("BGC_dist2",ExcelFileName = "BGC_dist2.xlsx",row.names = F,col.names = T)

#this is beyond the dimension of excel, so cut the rows and two files and then merge

BGC_dist1<-BGC_dist[1:966050,]
BGC_dist2<-BGC_dist[1:966051,]








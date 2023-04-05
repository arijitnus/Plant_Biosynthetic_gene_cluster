library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)


library(readxl)
library(dplyr)
library(ggplot2)

plant_isolates<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet="plant_isolates",col_names=T,skip = 0)
head(plant_isolates)
plant_isolates<-as.data.frame(plant_isolates)
rownames(plant_isolates)<-plant_isolates$Genome_ID
plant_isolates<-plant_isolates[,-1]
head(plant_isolates)
plant_isolates_tree<-read.tree("plantisolate_iq.tre")
set.seed(1231)
#Reorder the dataframe based on the tip labels of the tree
plant_isolates$genomes<-rownames(plant_isolates)
plant_isolates_ord<-plant_isolates[match(plant_isolates_tree$tip.label,plant_isolates$genomes),]
#arylpolyene
res_aryl<-phylosig(plant_isolates_tree,plant_isolates_ord$arylpolyene, method="lambda", test=TRUE, nsim=999)
res_aryl

#betalactone
res_betalactone<-phylosig(plant_isolates_tree,plant_isolates$betalactone, method="lambda", test=TRUE, nsim=999)
res_betalactone

#hserlactone
res_hserlactone<-phylosig(plant_isolates_tree,plant_isolates$hserlactone, method="lambda", test=TRUE, nsim=999)
res_hserlactone0.852
#NRPS
res_NRPS<-phylosig(plant_isolates_tree,plant_isolates$NRPS, method="lambda", test=TRUE, nsim=999)
res_NRPS
#RiPPs
res_RiPPs<-phylosig(plant_isolates_tree,plant_isolates$RiPPs, method="lambda", test=TRUE, nsim=999)
res_RiPPs

#siderophore
res_siderophore<-phylosig(plant_isolates_tree,plant_isolates$siderophore, method="lambda", test=TRUE, nsim=999)
res_siderophore

#terpene
res_terpene<-phylosig(plant_isolates_tree,plant_isolates$Terpene, method="lambda", test=TRUE, nsim=999)
res_terpene

#PKS_NRP_hybrids
res_PKS_NRP_hybrids<-phylosig(plant_isolates_tree,plant_isolates$PKS.NRP_Hybrids, method="lambda", test=TRUE, nsim=999)
res_PKS_NRP_hybrids

#PKS_Other
res_PKS_Other<-phylosig(plant_isolates_tree,plant_isolates$PKSother, method="lambda", test=TRUE, nsim=999)
res_PKS_Other

#PKSI
res_PKSI<-phylosig(plant_isolates_tree,plant_isolates$PKSI, method="lambda", test=TRUE, nsim=999)
res_PKSI

#Others
res_Others<-phylosig(plant_isolates_tree,plant_isolates$Others, method="lambda", test=TRUE, nsim=999)
res_Others



#=========================Soil isolates=============
##================================================
library(readxl)
soil_isolates<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet="soil_isolates",col_names=T,skip = 0)
head(soil_isolates)
soil_isolates<-as.data.frame(soil_isolates)
sum(is.na(soil_isolates))
#rownames(soil_isolates)<-soil_isolates$genomeID
#soil_isolates<-soil_isolates[,-1]
#head(soil_isolates)
soil_isolates_tree<-read.tree("soil_isolate_iq2.tre")
set.seed(1231)

soil_isolates_tree$tip.label<-labels_df$lable1

labels_df<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet="tips_change",col_names=T,skip = 0)
head(labels_df)#change the labels based on this
#Reorder the dataframe based on the tip labels of the tree
#soil_isolates$genomes<-rownames(soil_isolates)
soil_isolates_ord<-soil_isolates[match(soil_isolates_tree$tip.label,soil_isolates$genomeID),]
sum(soil_isolates_ord$genomeID==soil_isolates_tree$tip.label)
sum(is.na(soil_isolates_ord$genomeID))#74 genomes here are NA
write.tree(soil_isolates_tree,"soil_isolates_iqtree_new.tre")

#arylpolyene
res_aryl<-phylosig(soil_isolates_tree,soil_isolates_ord$arylpolyene, method="lambda", test=TRUE, nsim=999)
res_aryl

#betalactone
res_betalactone<-phylosig(soil_isolates_tree,soil_isolates_ord$betalactone, method="lambda", test=TRUE, nsim=999)
res_betalactone

#hserlactone
res_hserlactone<-phylosig(soil_isolates_tree,soil_isolates_ord$hserlactone, method="lambda", test=TRUE, nsim=999)
res_hserlactone

#NRPS
res_NRPS<-phylosig(soil_isolates_tree,soil_isolates_ord$NRPS, method="lambda", test=TRUE, nsim=999)
res_NRPS

#RiPPs
res_RiPPs<-phylosig(soil_isolates_tree,soil_isolates_ord$RiPPs, method="lambda", test=TRUE, nsim=999)
res_RiPPs

#siderophore
res_siderophore<-phylosig(soil_isolates_tree,soil_isolates_ord$siderophore, method="lambda", test=TRUE, nsim=999)
res_siderophore

#terpene
res_terpene<-phylosig(soil_isolates_tree,soil_isolates_ord$Terpene, method="lambda", test=TRUE, nsim=999)
res_terpene

#PKS-NRP_hybrids
res_PKS_NRP_hybrids<-phylosig(soil_isolates_tree,soil_isolates_ord$`PKS-NRP_Hybrids`, method="lambda", test=TRUE, nsim=999)
res_PKS_NRP_hybrids

#PKSOther
res_PKSother<-phylosig(soil_isolates_tree,soil_isolates_ord$PKSother, method="lambda", test=TRUE, nsim=999)
res_PKSother

#PKSI
res_PKSI<-phylosig(soil_isolates_tree,soil_isolates_ord$PKSI, method="lambda", test=TRUE, nsim=999)
res_PKSI

#Others
res_Others<-phylosig(soil_isolates_tree,soil_isolates_ord$Others, method="lambda", test=TRUE, nsim=999)
res_Others


#=========================##==============
##=========================================
#=============================================
#Plant MAGs
plant_MAGs<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet="plant_MAGs",col_names=T,skip = 0)
head(plant_MAGs)
plant_MAGs<-as.data.frame(plant_MAGs)
rownames(plant_MAGs)<-plant_MAGs$genome_id
plant_MAGs<-plant_MAGs[,-1]
head(plant_MAGs)
plant_MAGs_tree<-read.tree("plant_MAGsiq.tre")
set.seed(1231)
#Reorder the dataframe based on the tip labels of the tree
plant_MAGs$genomes<-rownames(plant_MAGs)
plant_MAGs_ord<-plant_MAGs[match(plant_MAGs_tree$tip.label,plant_MAGs$genomes),]
#arylpolyene
res_aryl<-phylosig(plant_MAGs_tree,plant_MAGs_ord$arylpolyene, method="lambda", test=TRUE, nsim=999)
res_aryl

#betalactone
res_betalactone<-phylosig(plant_MAGs_tree,plant_MAGs_ord$betalactone, method="lambda", test=TRUE, nsim=999)
res_betalactone

#hserlactone
res_hserlactone<-phylosig(plant_MAGs_tree,plant_MAGs_ord$hserlactone, method="lambda", test=TRUE, nsim=999)
res_hserlactone
#NRPS
res_NRPS<-phylosig(plant_MAGs_tree,plant_MAGs_ord$NRPS, method="lambda", test=TRUE, nsim=999)
res_NRPS
#RiPPs
res_RiPPs<-phylosig(plant_MAGs_tree,plant_MAGs_ord$RiPPs, method="lambda", test=TRUE, nsim=999)
res_RiPPs

#siderophore
res_siderophore<-phylosig(plant_MAGs_tree,plant_MAGs_ord$siderophore, method="lambda", test=TRUE, nsim=999)
res_siderophore

#terpene
res_terpene<-phylosig(plant_MAGs_tree,plant_MAGs_ord$Terpene, method="lambda", test=TRUE, nsim=999)
res_terpene

#PKS_NRP_hybrids
res_PKS_NRP_hybrids<-phylosig(plant_MAGs_tree,plant_MAGs_ord$`PKS-NRP_Hybrids`, method="lambda", test=TRUE, nsim=999)
res_PKS_NRP_hybrids

#PKS_Other
res_PKS_Other<-phylosig(plant_MAGs_tree,plant_MAGs_ord$PKSother, method="lambda", test=TRUE, nsim=999)
res_PKS_Other

#PKSI
res_PKSI<-phylosig(plant_MAGs_tree,plant_MAGs_ord$PKSI, method="lambda", test=TRUE, nsim=999)
res_PKSI

#Others
res_Others<-phylosig(plant_MAGs_tree,plant_MAGs_ord$Others, method="lambda", test=TRUE, nsim=999)
res_Others

#========================###==========
#======#================Soil MAGs=====
#=====================Soil MAGs=======
soil_MAGs<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet="Soil_MAGs",col_names=T,skip = 0)
head(soil_MAGs)
soil_MAGs<-as.data.frame(soil_MAGs)
rownames(soil_MAGs)<-soil_MAGs$genome_id
soil_MAGs<-soil_MAGs[,-1]
head(soil_MAGs)
soil_MAGs_tree<-read.tree("soil_MAGs_final_iq.tre")
set.seed(1231)
#Reorder the dataframe based on the tip labels of the tree
soil_MAGs$genomes<-rownames(soil_MAGs)
soil_MAGs_ord<-soil_MAGs[match(soil_MAGs_tree$tip.label,soil_MAGs$genomes),]
#arylpolyene
res_aryl<-phylosig(soil_MAGs_tree,soil_MAGs_ord$arylpolyene, method="lambda", test=TRUE, nsim=999)
res_aryl

#betalactone
res_betalactone<-phylosig(soil_MAGs_tree,soil_MAGs_ord$betalactone, method="lambda", test=TRUE, nsim=999)
res_betalactone

#hserlactone
res_hserlactone<-phylosig(soil_MAGs_tree,soil_MAGs_ord$hserlactone, method="lambda", test=TRUE, nsim=999)
res_hserlactone
#NRPS
res_NRPS<-phylosig(soil_MAGs_tree,soil_MAGs_ord$NRPS, method="lambda", test=TRUE, nsim=999)
res_NRPS
#RiPPs
res_RiPPs<-phylosig(soil_MAGs_tree,soil_MAGs_ord$RiPPs, method="lambda", test=TRUE, nsim=999)
res_RiPPs

#siderophore
res_siderophore<-phylosig(soil_MAGs_tree,soil_MAGs_ord$siderophore, method="lambda", test=TRUE, nsim=999)
res_siderophore

#terpene
res_terpene<-phylosig(soil_MAGs_tree,soil_MAGs_ord$Terpene, method="lambda", test=TRUE, nsim=999)
res_terpene

#PKS_NRP_hybrids
res_PKS_NRP_hybrids<-phylosig(soil_MAGs_tree,soil_MAGs_ord$PKS.NRP_Hybrids, method="lambda", test=TRUE, nsim=999)
res_PKS_NRP_hybrids

#PKS_Other
res_PKS_Other<-phylosig(soil_MAGs_tree,soil_MAGs_ord$PKSother, method="lambda", test=TRUE, nsim=999)
res_PKS_Other

#PKSI
res_PKSI<-phylosig(soil_MAGs_tree,soil_MAGs_ord$PKSI, method="lambda", test=TRUE, nsim=999)
res_PKSI

#Others
res_Others<-phylosig(soil_MAGs_tree,soil_MAGs_ord$Others, method="lambda", test=TRUE, nsim=999)
res_Others













#Calculate D-statistics for all datasets based on tree from iqtree
library(readxl)
library(dplyr)
library(caper)
soil_isolates<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",
                col_names = T,skip = 0,sheet = "soil_isolates")
head(soil_isolates)

soil_isolate_iqtree<-read.tree("soilisolate_iq.tre")
length(soil_isolate_iqtree$tip.label)

#make the genomes as rownames for soil isolate data
soil_isolates<-as.data.frame(soil_isolates)
rownames(soil_isolates)<-soil_isolates$Genome_ID
soil_isolates<-soil_isolates[,-1]
head(soil_isolates)
#Make it a binary data
soil_isolates_bin<-soil_isolates%>%mutate_if(is.numeric,~1*(.>0))
head(soil_isolates_bin)
soil_isolates_bin$genomes<-rownames(soil_isolates_bin)
class(soil_isolates_bin)

#calculate D statistics
soil_isolates_phy_signal<-comparative.data(soil_isolate_iqtree,soil_isolates_bin,genomes)

soil_isolates_aryl<-phylo.d(soil_isolates_phy_signal,binvar = arylpolyene)
soil_isolates_aryl#0.11; random= 0, Brownian = 0.04

soil_isolates_betalactone<-phylo.d(soil_isolates_phy_signal,binvar = betalactone)
soil_isolates_betalactone#0.17; random= 0, Brownian = 0

soil_isolates_hserlactone<-phylo.d(soil_isolates_phy_signal,binvar = hserlactone)
soil_isolates_hserlactone#0.14 ; random=0, Brownian = 0.015

soil_isolates_NRPS<-phylo.d(soil_isolates_phy_signal,binvar = NRPS)
soil_isolates_NRPS#0.35; 0; 0

soil_isolates_RiPPs<-phylo.d(soil_isolates_phy_signal,binvar = RiPPs)
soil_isolates_RiPPs#0.35; 0; 0 

soil_isolates_siderophore<-phylo.d(soil_isolates_phy_signal,binvar = siderophore)
soil_isolates_siderophore#0.25; 0; 0

soil_isolates_Terpene<-phylo.d(soil_isolates_phy_signal,binvar = Terpene)
soil_isolates_Terpene#0.06; 0; 0.183

soil_isolates_PKS.NRP_Hybrids<-phylo.d(soil_isolates_phy_signal,binvar = PKS.NRP_Hybrids)
soil_isolates_PKS.NRP_Hybrids#0.46;0;0

soil_isolates_PKSother<-phylo.d(soil_isolates_phy_signal,binvar = PKSother)
soil_isolates_PKSother#0.13;0;0.009


soil_isolates_PKSI<-phylo.d(soil_isolates_phy_signal,binvar = PKSI)
soil_isolates_PKSI#0.13; 0; 0.038

soil_isolates_Others<-phylo.d(soil_isolates_phy_signal,binvar = Others)
soil_isolates_Others#0.32;0;0




####Doing it for soil isolates data
library(readxl)
library(dplyr)
library(caper)
soil_isolates<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",
                           col_names = T,skip = 0,sheet = "soil_isolates")
head(soil_isolates)

soil_isolate_iqtree<-read.tree("soil_isolate_iq.tre")
length(soil_isolate_iqtree$tip.label)
#make the genomes as rownames for soil isolate data
soil_isolates<-as.data.frame(soil_isolates)
rownames(soil_isolates)<-soil_isolates$genomeID
soil_isolates<-soil_isolates[,-1]
head(soil_isolates)
#Make it a binary data
soil_isolates_bin<-soil_isolates%>%mutate_if(is.numeric,~1*(.>0))
head(soil_isolates_bin)
soil_isolates_bin$genomes<-rownames(soil_isolates_bin)
class(soil_isolates_bin)
head(soil_isolates)
soil_isolates$genomes<-rownames(soil_isolates)


soil_isolates_phy_signal<-comparative.data(soil_isolate_iqtree,soil_isolates_bin,genomes)


#arylpolyenes
soil_isolates_aryl<-phylo.d(soil_isolates_phy_signal,binvar = arylpolyene)
soil_isolates_aryl#0.25; random= 0, Brownian = 0.002

soil_isolates_betalactone<-phylo.d(soil_isolates_phy_signal,binvar = betalactone)
soil_isolates_betalactone#0.45; random= 0, Brownian = 0

soil_isolates_hserlactone<-phylo.d(soil_isolates_phy_signal,binvar = hserlactone)
soil_isolates_hserlactone#0.11 ; random=0, Brownian = 0.121

soil_isolates_NRPS<-phylo.d(soil_isolates_phy_signal,binvar = NRPS)
soil_isolates_NRPS#0.42; 0; 0

soil_isolates_RiPPs<-phylo.d(soil_isolates_phy_signal,binvar = RiPPs)
soil_isolates_RiPPs#0.51; 0; 0 

soil_isolates_siderophore<-phylo.d(soil_isolates_phy_signal,binvar = siderophore)
soil_isolates_siderophore#0.27; 0; 0

soil_isolates_Terpene<-phylo.d(soil_isolates_phy_signal,binvar = Terpene)
soil_isolates_Terpene#0.20; 0; 0.005
soil_isolates$`PKS-NRP_Hybrids`

soil_isolates_PKS.NRP_Hybrids<-phylo.d(soil_isolates_phy_signal,binvar = `PKS-NRP_Hybrids`)
soil_isolates_PKS.NRP_Hybrids#0.42;0;0

soil_isolates_PKSother<-phylo.d(soil_isolates_phy_signal,binvar = PKSother)
soil_isolates_PKSother#0.31;0;0.0


soil_isolates_PKSI<-phylo.d(soil_isolates_phy_signal,binvar = PKSI)
soil_isolates_PKSI#0.18; 0; 0.008

soil_isolates_Others<-phylo.d(soil_isolates_phy_signal,binvar = Others)
soil_isolates_Others#0.38;0;0

####soil MAGs
library(readxl)
library(dplyr)
library(caper)
soil_MAGs<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",
                          col_names = T,skip = 0,sheet = "soil_MAGs")
head(soil_MAGs)



soil_MAGs_iqtree<-read.tree("soil_MAGs.tre")
length(soil_MAGs_iqtree$tip.label)
#make the genomes as rownames for soil isolate data
soil_MAGs<-as.data.frame(soil_MAGs)
rownames(soil_MAGs)<-soil_MAGs$genome_id
soil_MAGs<-soil_MAGs[,-1]
head(soil_MAGs)
#Make it a binary data
soil_MAGs_bin<-soil_MAGs%>%mutate_if(is.numeric,~1*(.>0))
head(soil_MAGs_bin)
soil_MAGs_bin$genomes<-rownames(soil_MAGs_bin)
class(soil_MAGs_bin)
soil_MAGs_phy_signal<-comparative.data(soil_MAGs_iqtree,soil_MAGs_bin,genomes)

#arlypolyene
soil_MAGs_aryl<-phylo.d(soil_MAGs_phy_signal,binvar = arylpolyene)
soil_MAGs_aryl#0.55;0;0
#betalactone
soil_MAGs_betalactone<-phylo.d(soil_MAGs_phy_signal,binvar = betalactone)
soil_MAGs_betalactone#1.00;0.505;0

#hserlactone
soil_MAGs_hserlactone<-phylo.d(soil_MAGs_phy_signal,binvar = hserlactone)
soil_MAGs_hserlactone#0.69;0;0
#NRPS
soil_MAGs_NRPS<-phylo.d(soil_MAGs_phy_signal,binvar = NRPS)
soil_MAGs_NRPS#0.709;0;0
#RiPPs
soil_MAGs_RiPPs<-phylo.d(soil_MAGs_phy_signal,binvar = RiPPs)
soil_MAGs_RiPPs#0.85;0;0
#siderophore
soil_MAGs_siderophore<-phylo.d(soil_MAGs_phy_signal,binvar = siderophore)
soil_MAGs_siderophore#0.81;0;0

#Terpene
soil_MAGs_Terpene<-phylo.d(soil_MAGs_phy_signal,binvar = Terpene)
soil_MAGs_Terpene#0.84;0;0

#PKS.NRP_Hybrids
soil_MAGs_PKS.NRP_Hybrids<-phylo.d(soil_MAGs_phy_signal,binvar = `PKS-NRP_Hybrids`)
soil_MAGs_PKS.NRP_Hybrids#0.76;0;0

#PKSother
soil_MAGs_PKSother<-phylo.d(soil_MAGs_phy_signal,binvar = PKSother)
soil_MAGs_PKSother#0.79;0;0
#PKSI
soil_MAGs_PKSI<-phylo.d(soil_MAGs_phy_signal,binvar = PKSI)
soil_MAGs_PKSI#0.75;0;0

#Others
soil_MAGs_Others<-phylo.d(soil_MAGs_phy_signal,binvar = Others)
soil_MAGs_Others#0.0.79;0;0



####plant MAGs
library(readxl)
library(dplyr)
library(caper)
plant_MAGs<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",
                      col_names = T,skip = 0,sheet = "plant_MAGs")
head(plant_MAGs)



plant_MAGs_iqtree<-read.tree("plant_MAGs.tre")
length(plant_MAGs_iqtree$tip.label)
#make the genomes as rownames for soil isolate data
plant_MAGs<-as.data.frame(plant_MAGs)
rownames(plant_MAGs)<-plant_MAGs$genome_id
plant_MAGs<-plant_MAGs[,-1]
head(plant_MAGs)
#Make it a binary data
plant_MAGs_bin<-plant_MAGs%>%mutate_if(is.numeric,~1*(.>0))
head(plant_MAGs_bin)
plant_MAGs_bin$genomes<-rownames(plant_MAGs_bin)
class(plant_MAGs_bin)
plant_MAGs_phy_signal<-comparative.data(plant_MAGs_iqtree,plant_MAGs_bin,genomes)

#arlypolyene
plant_MAGs_aryl<-phylo.d(plant_MAGs_phy_signal,binvar = arylpolyene)
plant_MAGs_aryl#0.55;0;0
#betalactone
plant_MAGs_betalactone<-phylo.d(plant_MAGs_phy_signal,binvar = betalactone)
plant_MAGs_betalactone#1;0.5;0

#hserlactone
plant_MAGs_hserlactone<-phylo.d(plant_MAGs_phy_signal,binvar = hserlactone)
plant_MAGs_hserlactone#0.7;0;0
#NRPS
plant_MAGs_NRPS<-phylo.d(plant_MAGs_phy_signal,binvar = NRPS)
plant_MAGs_NRPS#0.78;0;0
#RiPPs
plant_MAGs_RiPPs<-phylo.d(plant_MAGs_phy_signal,binvar = RiPPs)
plant_MAGs_RiPPs#0.77;0;0
#siderophore
plant_MAGs_siderophore<-phylo.d(plant_MAGs_phy_signal,binvar = siderophore)
plant_MAGs_siderophore#0.88;0;0

#Terpene
plant_MAGs_Terpene<-phylo.d(plant_MAGs_phy_signal,binvar = Terpene)
plant_MAGs_Terpene#0.66;0;0

#PKS.NRP_Hybrids
plant_MAGs_PKS.NRP_Hybrids<-phylo.d(plant_MAGs_phy_signal,binvar = PKS.NRP_Hybrids)
plant_MAGs_PKS.NRP_Hybrids#0.88;0;0

#PKSother
plant_MAGs_PKSother<-phylo.d(plant_MAGs_phy_signal,binvar = PKSother)
plant_MAGs_PKSother#0.76;0;0
#PKSI
plant_MAGs_PKSI<-phylo.d(plant_MAGs_phy_signal,binvar = PKSI)
plant_MAGs_PKSI#.77;0;0

#Others
plant_MAGs_Others<-phylo.d(plant_MAGs_phy_signal,binvar = Others)
plant_MAGs_Others#0.8;0;0

















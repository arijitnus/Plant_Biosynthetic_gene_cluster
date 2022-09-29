#This is the script for calculating D-values 


library(readxl)
library(dplyr)
library(caper)
dat<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/new_report (1).xlsx",
                col_names = T,skip = 0,sheet = "isolate_BGC_table")
dat_98<-dat%>%filter(ANI_98=="Y")

dim(dat_98)#1395

#Read trees
tree_98<-read.tree("plant_isolates.tre")
#subset the dataframes based on the MAGs and only the BGC columns
head(dat_98)
dat_98_sub<-dat_98[,6:17]

####for the 98 ANI data
dat_98_sub_binary<-dat_98_sub%>% mutate_if(is.numeric, ~1 * (. > 0))
head(dat_98_sub_binary)
dat_98_sub_binary$genomes<-dat_98$Genome
dat_98_sub_binary<-as.data.frame(dat_98_sub_binary)
length(tree_98$tip.label)
tree_98$node.label<-NULL
dat_98_phy_signal<-comparative.data(tree_98,dat_98_sub_binary,genomes)

#calculate D-values
dat_98_aryl<-phylo.d(dat_98_phy_signal,binvar = arylpolyene)
dat_98_aryl#0.09

dat_98_betalactone<-phylo.d(dat_98_phy_signal,binvar = betalactone)
dat_98_betalactone#0.16

dat_98_hserlactone<-phylo.d(dat_98_phy_signal,binvar = hserlactone)
dat_98_hserlactone#0.12

dat_98_NRPS<-phylo.d(dat_98_phy_signal,binvar = NRPS)
dat_98_NRPS#0.35


dat_98_RiPPs<-phylo.d(dat_98_phy_signal,binvar = RiPPs)
dat_98_RiPPs#0.34

dat_98_siderophore<-phylo.d(dat_98_phy_signal,binvar = siderophore)
dat_98_siderophore#0.25


dat_98_terpene<-phylo.d(dat_98_phy_signal,binvar = Terpene)
dat_98_terpene#0.05

dat_98_PKS_NRP_hybrids<-phylo.d(dat_98_phy_signal,binvar = PKS.NRP_Hybrids)
dat_98_PKS_NRP_hybrids#0.45

dat_98_PKS_other<-phylo.d(dat_98_phy_signal,binvar = PKSother)
dat_98_PKS_other#0.125

dat_98_PKSI<-phylo.d(dat_98_phy_signal,binvar = PKSI)
dat_98_PKSI#0.126

dat_98_ectoine<-phylo.d(dat_98_phy_signal,binvar = ectoine)
dat_98_ectoine#0.088

dat_98_other<-phylo.d(dat_98_phy_signal,binvar = Others)
dat_98_other#0.33


#######===================Plant MAGs all three ANIs; Calculate D values=======


library(readxl)
library(dplyr)
library(caper)
dat_plant_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/new_report (1).xlsx",
                           col_names = T,skip = 0,sheet = "Plant_MAGs_curated")
head(dat_plant_mags)
plant_mags_98<-dat_plant_mags%>%filter(ANI_98=="Y")


dim(plant_mags_98)#573


plant_tree_98<-read.tree("plantmags.tre")


names(plant_mags_98)
plant_mags_98_sub<-plant_mags_98[,21:31]


plant_mags_98_sub_binary<-plant_mags_98_sub%>% mutate_if(is.numeric, ~1 * (. > 0))
head(plant_mags_98_sub_binary)
plant_mags_98_sub_binary$genomes<-plant_mags_98$genome_id
plant_mags_98_sub_binary<-as.data.frame(plant_mags_98_sub_binary)
length(plant_tree_98$tip.label)
plant_tree_98$node.label<-NULL
plant_98_phy_signal<-comparative.data(plant_tree_98,plant_mags_98_sub_binary,genomes)

#calculate D-values
plant_98_aryl<-phylo.d(plant_98_phy_signal,binvar = arylpolyene)
plant_98_aryl#0.58

plant_98_betalactone<-phylo.d(plant_98_phy_signal,binvar = betalactone)
plant_98_betalactone#0.93

plant_98_hserlactone<-phylo.d(plant_98_phy_signal,binvar = hserlactone)
plant_98_hserlactone#0.7

plant_98_NRPS<-phylo.d(plant_98_phy_signal,binvar = NRPS)
plant_98_NRPS#0.77

plant_98_RiPPs<-phylo.d(plant_98_phy_signal,binvar = RiPPs)
plant_98_RiPPs#0.77

plant_98_siderophore<-phylo.d(plant_98_phy_signal,binvar = siderophore)
plant_98_siderophore#0.89

plant_98_terpene<-phylo.d(plant_98_phy_signal,binvar = Terpene)
plant_98_terpene#0.68

plant_98_PKS_NRP_Hybrids<-phylo.d(plant_98_phy_signal,binvar =`PKS-NRP_Hybrids` )
plant_98_PKS_NRP_Hybrids#0.89

plant_98_PKSother<-phylo.d(plant_98_phy_signal,binvar = PKSother )
plant_98_PKSother#0.77

plant_98_PKSI<-phylo.d(plant_98_phy_signal,binvar = PKSI )
plant_98_PKSI#0.82

plant_98_other<-phylo.d(plant_98_phy_signal,binvar = Others)
plant_98_other#0.8



#######===================Soil MAGs all three ANIs; Calculate D values=======


library(readxl)
library(dplyr)
library(caper)
dat_soil_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/new_report (1).xlsx",
                          col_names = T,skip = 0,sheet = "Soil_MAGs_curated")
head(dat_soil_mags)
soil_mags_98<-dat_soil_mags%>%filter(ANI_98=="Y")


dim(soil_mags_98)#783


soil_tree_98<-read.tree("gtdbtk_rooted_98_soilmags.tre")


names(soil_mags_98)
soil_mags_98_sub<-soil_mags_98[,20:30]


soil_mags_98_sub_binary<-soil_mags_98_sub%>% mutate_if(is.numeric, ~1 * (. > 0))
head(soil_mags_98_sub_binary)
soil_mags_98_sub_binary$genomes<-soil_mags_98$genome_id
soil_mags_98_sub_binary<-as.data.frame(soil_mags_98_sub_binary)
length(soil_tree_98$tip.label)
soil_tree_98$node.label<-NULL
soil_98_phy_signal<-comparative.data(soil_tree_98,soil_mags_98_sub_binary,genomes)


#calculate D-values at 98 ANI
soil_98_aryl<-phylo.d(soil_98_phy_signal,binvar = arylpolyene)
soil_98_aryl#0.74

soil_98_betalactone<-phylo.d(soil_98_phy_signal,binvar = betalactone)
soil_98_betalactone#0.86

soil_98_hserlactone<-phylo.d(soil_98_phy_signal,binvar = hserlactone)
soil_98_hserlactone#0.65

soil_98_NRPS<-phylo.d(soil_98_phy_signal,binvar = NRPS)
soil_98_NRPS#0.7

soil_98_RiPPs<-phylo.d(soil_98_phy_signal,binvar = RiPPs)
soil_98_RiPPs#0.84

soil_98_siderophore<-phylo.d(soil_98_phy_signal,binvar = siderophore)
soil_98_siderophore#0.81

soil_98_terpene<-phylo.d(soil_98_phy_signal,binvar = Terpene)
soil_98_terpene#0.85

soil_98_PKS_NRP_Hybrids<-phylo.d(soil_98_phy_signal,binvar =PKS.NRP_Hybrids )
soil_98_PKS_NRP_Hybrids#0.76


soil_98_PKSother<-phylo.d(soil_98_phy_signal,binvar =PKSother )
soil_98_PKSother#0.77

soil_98_PKSI<-phylo.d(soil_98_phy_signal,binvar =PKSI )
soil_98_PKSI#0.75

soil_98_others<-phylo.d(soil_98_phy_signal,binvar =Others )
soil_98_others#0.78













#This is the script for calculating Phylogenetic D values for each of the BGC category
#MAGS data
library(dplyr)
library(caper)
library(readxl)
MAGs<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/MAGS_allinfo.xlsx",
                              sheet = "final",col_names = T,skip = 0)
names(MAGs)
#Read MAGs tre file
MAGtree<-read.tree("zerobranchlength.tre")
MAGs_sub<-MAGs[,10:20]
MAGs_sub$Assembly<-MAGs$MAGs
names(MAGs_sub)
MAGs_sub_binary<-MAGs_sub%>% mutate_if(is.numeric, ~1 * (. > 0))
dim(MAGs_sub_binary)
MAGs_sub_binary<-as.data.frame(MAGs_sub_binary)
length(MAGtree$tip.label)
MAGtree$node.label<-NULL
MAGtree_assmebly<-MAGtree$tip.label
MAGs_sub_binary<-MAGs_sub_binary%>%filter(Assembly%in%MAGtree_assmebly)
dim(MAGs_sub_binary)#there is a problem of number of 
MAG_phy_signal<-comparative.data(MAGtree,MAGs_sub_binary,Assembly)

MAG_aryl<-phylo.d(MAG_phy_signal,binvar = arylpolyene)
MAG_betalactone<-phylo.d(MAG_phy_signal,binvar = betalactone)
MAG_hserlactone<-phylo.d(MAG_phy_signal,binvar = hserlactone)
MAG_NRPS<-phylo.d(MAG_phy_signal,binvar = NRPS)
MAG_others<-phylo.d(MAG_phy_signal,binvar = Others)
MAGs_sub_binary$`PKS-NRP_Hybrids`
MAG_PKSI<-phylo.d(MAG_phy_signal,binvar = PKSI)
MAG_PKSother<-phylo.d(MAG_phy_signal,binvar = PKSother)
MAG_RiPPs<-phylo.d(MAG_phy_signal,binvar = RiPPs)
MAG_siderophore<-phylo.d(MAG_phy_signal,binvar = siderophore)
MAG_terpene<-phylo.d(MAG_phy_signal,binvar = Terpene)

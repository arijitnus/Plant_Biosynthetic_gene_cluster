# The following is the script for consenTRAIT analysis of soil isolates data##
library(readxl)
library(dplyr)
library(caper)
library(castor)
tree_98<-read.tree("soil_98.tree")


#Read in the 98 file for plant MAGs
soil_98<-read_excel("/Users/arijitmukherjee/Documents/soil_D_stat/soil_isolates.xlsx",
                    col_names = T,skip = 0,sheet = "table")
aryl_98<-get_trait_depth(tree_98,
                         soil_98$arylpolyene,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
aryl_98$mean_depth#0.030;P=0.11
aryl_98$P
aryl_98

terpene_98<-get_trait_depth(tree_98,
                            soil_98$Terpene,
                            min_fraction = 0.9,
                            count_singletons = TRUE,
                            singleton_resolution= 0,
                            weighted = FALSE,
                            Npermutations = 1000)
terpene_98$mean_depth#0.050; P=0.001
terpene_98$P
terpene_98



betalactone_98<-get_trait_depth(tree_98,
                                soil_98$betalactone,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
betalactone_98$mean_depth#0.033
betalactone_98$P#0.006
betalactone_98
hserlactone_98<-get_trait_depth(tree_98,
                                soil_98$hserlactone,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
hserlactone_98$mean_depth#0.032; P=0.006
hserlactone_98$P
hserlactone_98

NRPS_98<-get_trait_depth(tree_98,
                         soil_98$NRPS,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
NRPS_98$mean_depth#0.042; P=0
NRPS_98$P
NRPS_98
PKS_NRPhybrids_98<-get_trait_depth(tree_98,
                                   soil_98$`PKS-NRP_Hybrids`,
                                   min_fraction = 0.9,
                                   count_singletons = TRUE,
                                   singleton_resolution= 0,
                                   weighted = FALSE,
                                   Npermutations = 1000)
PKS_NRPhybrids_98$mean_depth#0.032; P=0.019
PKS_NRPhybrids_98$P
PKS_NRPhybrids_98

PKSI_98<-get_trait_depth(tree_98,
                         soil_98$PKSI,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
PKSI_98$mean_depth#0.035; P=0
PKSI_98$P
PKSI_98

PKSother_98<-get_trait_depth(tree_98,
                             soil_98$PKSother,
                             min_fraction = 0.9,
                             count_singletons = TRUE,
                             singleton_resolution= 0,
                             weighted = FALSE,
                             Npermutations = 1000)
PKSother_98$mean_depth#0.040; P=0
PKSother_98$P
PKSother_98
RiPPs_98<-get_trait_depth(tree_98,
                          soil_98$RiPPs,
                          min_fraction = 0.9,
                          count_singletons = TRUE,
                          singleton_resolution= 0,
                          weighted = FALSE,
                          Npermutations = 1000)
RiPPs_98$mean_depth#0.044; P=0
RiPPs_98$P
RiPPs_98

siderophore_98<-get_trait_depth(tree_98,
                                soil_98$siderophore,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
siderophore_98$mean_depth#0.035; P=0
siderophore_98$P
siderophore_98
others_98<-get_trait_depth(tree_98,
                           soil_98$Others,
                           min_fraction = 0.9,
                           count_singletons = TRUE,
                           singleton_resolution= 0,
                           weighted = FALSE,
                           Npermutations = 1000)
others_98$mean_depth#0.040; P=0
others_98$P
others_98$max_depth


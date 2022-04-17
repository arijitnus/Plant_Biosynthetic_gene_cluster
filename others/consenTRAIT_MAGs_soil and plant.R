# The following is the script for consenTRAIT analysis of MAGs data##
##the input file is Supplementary table  for plant MAGs and soil MAGs
library(readxl)
library(dplyr)
library(caper)
library(castor)


plant_tree_98<-read.tree("gtdbtk_rooted_98_plantmags.tre")


#Read in the 98 file for plant MAGs
plant_mags_98_matched<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/Fig.2/new_report_consentrait.xlsx",col_names = T,skip = 0,sheet = "Plant_MAGs_98")

aryl_98<-get_trait_depth(plant_tree_98,
                         plant_mags_98_matched$arylpolyene,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
aryl_98$mean_depth#0.034;P=0.997
aryl_98$P
aryl_98

terpene_98<-get_trait_depth(plant_tree_98,
                            plant_mags_98_matched$Terpene,
                            min_fraction = 0.9,
                            count_singletons = TRUE,
                            singleton_resolution= 0,
                            weighted = FALSE,
                            Npermutations = 1000)
terpene_98$mean_depth#0.053; P=0.089
terpene_98$P
terpene_98
betalactone_98<-get_trait_depth(plant_tree_98,
                                plant_mags_98_matched$betalactone,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
betalactone_98$mean_depth#NaN
betalactone_98$P
betalactone_98
hserlactone_98<-get_trait_depth(plant_tree_98,
                                plant_mags_98_matched$hserlactone,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
hserlactone_98$mean_depth#0.031; P=0.26
hserlactone_98$P
hserlactone_98

NRPS_98<-get_trait_depth(plant_tree_98,
                         plant_mags_98_matched$NRPS,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
NRPS_98$mean_depth#0.091; P=0.141
NRPS_98$P
NRPS_98
PKS_NRPhybrids_98<-get_trait_depth(plant_tree_98,
                                   plant_mags_98_matched$`PKS-NRP_Hybrids`,
                                   min_fraction = 0.9,
                                   count_singletons = TRUE,
                                   singleton_resolution= 0,
                                   weighted = FALSE,
                                   Npermutations = 1000)
PKS_NRPhybrids_98$mean_depth#0.045; P=0.26
PKS_NRPhybrids_98$P
PKS_NRPhybrids_98

PKSI_98<-get_trait_depth(plant_tree_98,
                         plant_mags_98_matched$PKSI,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
PKSI_98$mean_depth#0.06; P=0.178
PKSI_98$P
PKSI_98

PKSother_98<-get_trait_depth(plant_tree_98,
                             plant_mags_98_matched$PKSother,
                             min_fraction = 0.9,
                             count_singletons = TRUE,
                             singleton_resolution= 0,
                             weighted = FALSE,
                             Npermutations = 1000)
PKSother_98$mean_depth#0.049; P=0.6
PKSother_98$P
PKSother_98
RiPPs_98<-get_trait_depth(plant_tree_98,
                          plant_mags_98_matched$RiPPs,
                          min_fraction = 0.9,
                          count_singletons = TRUE,
                          singleton_resolution= 0,
                          weighted = FALSE,
                          Npermutations = 1000)
RiPPs_98$mean_depth#0.07; P=0.313
RiPPs_98$P
RiPPs_98

siderophore_98<-get_trait_depth(plant_tree_98,
                                plant_mags_98_matched$siderophore,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
siderophore_98$mean_depth
siderophore_98$P
siderophore_98
others_98<-get_trait_depth(plant_tree_98,
                           plant_mags_98_matched$Others,
                           min_fraction = 0.9,
                           count_singletons = TRUE,
                           singleton_resolution= 0,
                           weighted = FALSE,
                           Npermutations = 1000)
others_98$mean_depth
others_98$P
others_98$max_depth

library(readxl)
soil_mags_98_matched<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/Fig.2/new_report_consentrait.xlsx",col_names = T,skip = 0,sheet = "soil_MAGs_98")
soil_tree_98<-read.tree("gtdbtk_rooted_98_soilmags.tre")
#(soil_tree_98$tip.label,"soil_tips98.tsv",sep = "\t")
#Read in the ordered data


aryl_98_soil<-get_trait_depth(soil_tree_98,
                              soil_mags_98_matched$arylpolyene,
                              min_fraction = 0.9,
                              count_singletons = TRUE,
                              singleton_resolution= 0,
                              weighted = FALSE,
                              Npermutations = 1000)
aryl_98_soil$mean_depth#0.06;P=0.267
aryl_98_soil$P
aryl_98_soil$min_depth
aryl_98_soil$max_depth

betalactone_98_soil<-get_trait_depth(soil_tree_98,
                                     soil_mags_98_matched$betalactone,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
betalactone_98_soil$mean_depth#0.052;P=0.267
betalactone_98_soil$P
betalactone_98_soil$min_depth
betalactone_98_soil$max_depth#0.24

hserlactone_98_soil<-get_trait_depth(soil_tree_98,
                                     soil_mags_98_matched$hserlactone,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
hserlactone_98_soil$mean_depth#0.052;P=0.267
hserlactone_98_soil$P
hserlactone_98_soil$min_depth
hserlactone_98_soil$max_depth

NRPS_98_soil<-get_trait_depth(soil_tree_98,
                              soil_mags_98_matched$NRPS,
                              min_fraction = 0.9,
                              count_singletons = TRUE,
                              singleton_resolution= 0,
                              weighted = FALSE,
                              Npermutations = 1000)
NRPS_98_soil$mean_depth#0.076;P=0
NRPS_98_soil$P
NRPS_98_soil$min_depth
NRPS_98_soil$max_depth


RiPPs_98_soil<-get_trait_depth(soil_tree_98,
                               soil_mags_98_matched$RiPPs,
                               min_fraction = 0.9,
                               count_singletons = TRUE,
                               singleton_resolution= 0,
                               weighted = FALSE,
                               Npermutations = 1000)
RiPPs_98_soil$mean_depth#0.068;P=0.098
RiPPs_98_soil$P
RiPPs_98_soil$min_depth
RiPPs_98_soil$max_depth

siderophore_98_soil<-get_trait_depth(soil_tree_98,
                                     soil_mags_98_matched$siderophore,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
siderophore_98_soil$mean_depth#0.068;P=0.098
siderophore_98_soil$P
siderophore_98_soil$min_depth
siderophore_98_soil$max_depth


terpene_98_soil<-get_trait_depth(soil_tree_98,
                                 soil_mags_98_matched$Terpene,
                                 min_fraction = 0.9,
                                 count_singletons = TRUE,
                                 singleton_resolution= 0,
                                 weighted = FALSE,
                                 Npermutations = 1000)
terpene_98_soil$mean_depth#0.068;P=0.098
terpene_98_soil$P
terpene_98_soil$min_depth
terpene_98_soil$max_depth

PKS_NRP_hybrids_98_soil<-get_trait_depth(soil_tree_98,
                                         soil_mags_98_matched$PKS.NRP_Hybrids,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
PKS_NRP_hybrids_98_soil$mean_depth#0.068;P=0.098
PKS_NRP_hybrids_98_soil$P
PKS_NRP_hybrids_98_soil$min_depth
PKS_NRP_hybrids_98_soil$max_depth

PKSother_98_soil<-get_trait_depth(soil_tree_98,
                                  soil_mags_98_matched$PKSother,
                                  min_fraction = 0.9,
                                  count_singletons = TRUE,
                                  singleton_resolution= 0,
                                  weighted = FALSE,
                                  Npermutations = 1000)
PKSother_98_soil$mean_depth#0.068;P=0.098
PKSother_98_soil$P
PKSother_98_soil$min_depth
PKSother_98_soil$max_depth


PKSI_98_soil<-get_trait_depth(soil_tree_98,
                              soil_mags_98_matched$PKSI,
                              min_fraction = 0.9,
                              count_singletons = TRUE,
                              singleton_resolution= 0,
                              weighted = FALSE,
                              Npermutations = 1000)
PKSI_98_soil$mean_depth#0.068;P=0.098
PKSI_98_soil$P
PKSI_98_soil$min_depth
PKSI_98_soil$max_depth

Others_98_soil<-get_trait_depth(soil_tree_98,
                                soil_mags_98_matched$Others,
                                min_fraction = 0.9,
                                count_singletons = TRUE,
                                singleton_resolution= 0,
                                weighted = FALSE,
                                Npermutations = 1000)
Others_98_soil$mean_depth#0.068;P=0.098
Others_98_soil$P
Others_98_soil$min_depth
Others_98_soil$max_depth

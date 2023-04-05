library(readxl)
library(dplyr)
library(caper)
library(castor)


plant_isolate_tree<-read.tree("plantisolate_iq.tre")

#Plant isolates
plant_isolate<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet = "plant_isolates",col_names = T,skip = 0)
plant_isolate<-as.data.frame(plant_isolate)
rownames(plant_isolate)<-plant_isolate$Genome_ID
plant_isolate<-plant_isolate[,-1]
plant_isolate_bin<-plant_isolate%>%mutate_if(is.numeric,~1*(.>0))

aryl_plant_isolate<-consentrait_depth(plant_isolate_tree,
                         plant_isolate$arylpolyene,
                         min_fraction = 0.9,
                         count_singletons = TRUE,
                         singleton_resolution= 0,
                         weighted = FALSE,
                         Npermutations = 1000)
aryl_plant_isolate$mean_depth#0.020;P=0
aryl_plant_isolate$P
aryl_plant_isolate


betalactone_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                      plant_isolate$betalactone,
                                      min_fraction = 0.9,
                                      count_singletons = TRUE,
                                      singleton_resolution= 0,
                                      weighted = FALSE,
                                      Npermutations = 1000)
betalactone_plant_isolate$mean_depth#0.018;P=0
betalactone_plant_isolate$P
betalactone_plant_isolate


hserlactone_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                             plant_isolate$hserlactone,
                                             min_fraction = 0.9,
                                             count_singletons = TRUE,
                                             singleton_resolution= 0,
                                             weighted = FALSE,
                                             Npermutations = 1000)
hserlactone_plant_isolate$mean_depth#0.023;P=0
hserlactone_plant_isolate$P
hserlactone_plant_isolate

NRPS_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                             plant_isolate$NRPS,
                                             min_fraction = 0.9,
                                             count_singletons = TRUE,
                                             singleton_resolution= 0,
                                             weighted = FALSE,
                                             Npermutations = 1000)
NRPS_plant_isolate$mean_depth#0.026;P=0
NRPS_plant_isolate$P
NRPS_plant_isolate


RiPPs_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                      plant_isolate$RiPPs,
                                      min_fraction = 0.9,
                                      count_singletons = TRUE,
                                      singleton_resolution= 0,
                                      weighted = FALSE,
                                      Npermutations = 1000)
RiPPs_plant_isolate$mean_depth#0.032;P=0
RiPPs_plant_isolate$P
RiPPs_plant_isolate

siderophore_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                       plant_isolate$siderophore,
                                       min_fraction = 0.9,
                                       count_singletons = TRUE,
                                       singleton_resolution= 0,
                                       weighted = FALSE,
                                       Npermutations = 1000)
siderophore_plant_isolate$mean_depth#0.018;P=0
siderophore_plant_isolate$P
siderophore_plant_isolate

Terpene_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                             plant_isolate$Terpene,
                                             min_fraction = 0.9,
                                             count_singletons = TRUE,
                                             singleton_resolution= 0,
                                             weighted = FALSE,
                                             Npermutations = 1000)
Terpene_plant_isolate$mean_depth#0.036;P=0
Terpene_plant_isolate$P


Terpene_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                         plant_isolate$Terpene,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
Terpene_plant_isolate$mean_depth#0.036;P=0
Terpene_plant_isolate$P

PKS_NRP_hybrids_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                         plant_isolate$PKS.NRP_Hybrids,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
PKS_NRP_hybrids_plant_isolate$mean_depth#0.022;P=0
PKS_NRP_hybrids_plant_isolate$P

PKS_other_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                                 plant_isolate$PKSother,
                                                 min_fraction = 0.9,
                                                 count_singletons = TRUE,
                                                 singleton_resolution= 0,
                                                 weighted = FALSE,
                                                 Npermutations = 1000)
PKS_other_plant_isolate$mean_depth#0.022;P=0
PKS_other_plant_isolate$P


PKSI_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                           plant_isolate$PKSI,
                                           min_fraction = 0.9,
                                           count_singletons = TRUE,
                                           singleton_resolution= 0,
                                           weighted = FALSE,
                                           Npermutations = 1000)
PKSI_plant_isolate$mean_depth#0.016;P=0.001
PKSI_plant_isolate$P

Others_plant_isolate<-consentrait_depth(plant_isolate_tree,
                                      plant_isolate$Others,
                                      min_fraction = 0.9,
                                      count_singletons = TRUE,
                                      singleton_resolution= 0,
                                      weighted = FALSE,
                                      Npermutations = 1000)
Others_plant_isolate$mean_depth#0.019;P=0
Others_plant_isolate$P

######====================###==============
#Soil isolates:
#########==================================
soil_isolate<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet = "soil_isolates",col_names = T,skip = 0)
soil_isolate<-as.data.frame(soil_isolate)
rownames(soil_isolate)<-soil_isolate$genomeID
soil_isolate<-soil_isolate[,-1]
soil_isolate_bin<-soil_isolate%>%mutate_if(is.numeric,~1*(.>0))
soil_isolate_tree<-read.tree("soil_isolate_iq.tre")


aryl_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                     soil_isolate$arylpolyene,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
aryl_soil_isolate$mean_depth#0.021;P=0.996
aryl_soil_isolate$P
aryl_soil_isolate

betalactone_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                            soil_isolate$betalactone,
                                            min_fraction = 0.9,
                                            count_singletons = TRUE,
                                            singleton_resolution= 0,
                                            weighted = FALSE,
                                            Npermutations = 1000)
betalactone_soil_isolate$mean_depth#0.031;P=0.04
betalactone_soil_isolate$P
betalactone_soil_isolate

hserlactone_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                            soil_isolate$hserlactone,
                                            min_fraction = 0.9,
                                            count_singletons = TRUE,
                                            singleton_resolution= 0,
                                            weighted = FALSE,
                                            Npermutations = 1000)
hserlactone_soil_isolate$mean_depth#0.027;P=0.296
hserlactone_soil_isolate$P
hserlactone_soil_isolate

NRPS_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                     soil_isolate$NRPS,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
NRPS_soil_isolate$mean_depth#0.036;P=0.111
NRPS_soil_isolate$P
NRPS_soil_isolate

RiPPs_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                      soil_isolate$RiPPs,
                                      min_fraction = 0.9,
                                      count_singletons = TRUE,
                                      singleton_resolution= 0,
                                      weighted = FALSE,
                                      Npermutations = 1000)
RiPPs_soil_isolate$mean_depth#0.041;P=0.001
RiPPs_soil_isolate$P
RiPPs_soil_isolate


siderophore_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                            soil_isolate$siderophore,
                                            min_fraction = 0.9,
                                            count_singletons = TRUE,
                                            singleton_resolution= 0,
                                            weighted = FALSE,
                                            Npermutations = 1000)
siderophore_soil_isolate$mean_depth#0.028;P=0.639
siderophore_soil_isolate$P
siderophore_soil_isolate


Terpene_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                        soil_isolate$Terpene,
                                        min_fraction = 0.9,
                                        count_singletons = TRUE,
                                        singleton_resolution= 0,
                                        weighted = FALSE,
                                        Npermutations = 1000)
Terpene_soil_isolate$mean_depth#0.050;P=0
Terpene_soil_isolate$P

PKS_NRP_hybrids_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                                soil_isolate$`PKS-NRP_Hybrids`,
                                                min_fraction = 0.9,
                                                count_singletons = TRUE,
                                                singleton_resolution= 0,
                                                weighted = FALSE,
                                                Npermutations = 1000)
PKS_NRP_hybrids_soil_isolate$mean_depth#0.027;P=0.829
PKS_NRP_hybrids_soil_isolate$P

PKS_other_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                          soil_isolate$PKSother,
                                          min_fraction = 0.9,
                                          count_singletons = TRUE,
                                          singleton_resolution= 0,
                                          weighted = FALSE,
                                          Npermutations = 1000)
PKS_other_soil_isolate$mean_depth#0.034;P=0.023
PKS_other_soil_isolate$P

PKSI_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                     soil_isolate$PKSI,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
PKSI_soil_isolate$mean_depth#0.029;P=0.143
PKSI_soil_isolate$P

Others_soil_isolate<-consentrait_depth(soil_isolate_tree,
                                       soil_isolate$Others,
                                       min_fraction = 0.9,
                                       count_singletons = TRUE,
                                       singleton_resolution= 0,
                                       weighted = FALSE,
                                       Npermutations = 1000)
Others_soil_isolate$mean_depth#0.033;P=0.803
Others_soil_isolate$P


#########++++++++++++++++==============
#########+#=================================
#########+########

set.seed(1231)
#Plant MAGs
plant_MAG<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet = "plant_MAGs",col_names = T,skip = 0)
plant_MAG<-as.data.frame(plant_MAG)
rownames(plant_MAG)<-plant_MAG$genome_id
plant_MAG<-plant_MAG[,-1]
plant_MAG_bin<-plant_MAG %>% mutate_if(is.numeric,~1*(.>0))

plant_MAG_tree<-read.tree("plant_MAGsiq.tre")

aryl_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                  plant_MAG$arylpolyene,
                                  min_fraction = 0.9,
                                  count_singletons = TRUE,
                                  singleton_resolution= 0,
                                  weighted = FALSE,
                                  Npermutations = 1000)
aryl_plant_MAG$mean_depth#0.040;P=0.877
aryl_plant_MAG$P

betalactone_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                         plant_MAG$betalactone,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
betalactone_plant_MAG$mean_depth#0.052;P=0.24
betalactone_plant_MAG$P


hserlactone_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                         plant_MAG$hserlactone,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
hserlactone_plant_MAG$mean_depth#0.038;P=0.811
hserlactone_plant_MAG$P


NRPS_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                  plant_MAG$NRPS,
                                  min_fraction = 0.9,
                                  count_singletons = TRUE,
                                  singleton_resolution= 0,
                                  weighted = FALSE,
                                  Npermutations = 1000)
NRPS_plant_MAG$mean_depth#0.051;P=0.729
NRPS_plant_MAG$P


RiPPs_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                   plant_MAG$RiPPs,
                                   min_fraction = 0.9,
                                   count_singletons = TRUE,
                                   singleton_resolution= 0,
                                   weighted = FALSE,
                                   Npermutations = 1000)
RiPPs_plant_MAG$mean_depth#0.061;P=0.005
RiPPs_plant_MAG$P


siderophore_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                         plant_MAG$siderophore,
                                         min_fraction = 0.9,
                                         count_singletons = TRUE,
                                         singleton_resolution= 0,
                                         weighted = FALSE,
                                         Npermutations = 1000)
siderophore_plant_MAG$mean_depth#0.043;P=0.939
siderophore_plant_MAG$P

Terpene_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                     plant_MAG$Terpene,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
Terpene_plant_MAG$mean_depth#0.069;P=0
Terpene_plant_MAG$P

PKS_NRP_hybrids_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                     plant_MAG$`PKS-NRP_Hybrids`,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
PKS_NRP_hybrids_plant_MAG$mean_depth#0.054;P=0.09
PKS_NRP_hybrids_plant_MAG$P



PKSI_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                             plant_MAG$PKSI,
                                             min_fraction = 0.9,
                                             count_singletons = TRUE,
                                             singleton_resolution= 0,
                                             weighted = FALSE,
                                             Npermutations = 1000)
PKSI_plant_MAG$mean_depth#0.067;P=0.018
PKSI_plant_MAG$P


PKSother_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                  plant_MAG$PKSother,
                                  min_fraction = 0.9,
                                  count_singletons = TRUE,
                                  singleton_resolution= 0,
                                  weighted = FALSE,
                                  Npermutations = 1000)
PKSother_plant_MAG$mean_depth#0.039;P=0.9
PKSother_plant_MAG$P


Others_plant_MAG<-consentrait_depth(plant_MAG_tree,
                                      plant_MAG$Others,
                                      min_fraction = 0.9,
                                      count_singletons = TRUE,
                                      singleton_resolution= 0,
                                      weighted = FALSE,
                                      Npermutations = 1000)
Others_plant_MAG$mean_depth#0.060;P=0.044
Others_plant_MAG$P


##########===============
##====================
#soil MAGs
soil_MAG<-read_excel("/Users/arijitmukherjee/Downloads/D_test_iqtree.xlsx",sheet = "Soil_MAGs",col_names = T,skip = 0)
soil_MAG<-as.data.frame(soil_MAG)
rownames(soil_MAG)<-soil_MAG$genome_id
soil_MAG<-soil_MAG[,-1]
soil_MAG_bin<-soil_MAG %>% mutate_if(is.numeric,~1*(.>0))

soil_MAG_tree<-read.tree("soil_MAGs_final_iq.tre")

aryl_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                 soil_MAG$arylpolyene,
                                 min_fraction = 0.9,
                                 count_singletons = TRUE,
                                 singleton_resolution= 0,
                                 weighted = FALSE,
                                 Npermutations = 1000)
aryl_soil_MAG$mean_depth#0.049;P=0.847
aryl_soil_MAG$P

betalactone_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                        soil_MAG$betalactone,
                                        min_fraction = 0.9,
                                        count_singletons = TRUE,
                                        singleton_resolution= 0,
                                        weighted = FALSE,
                                        Npermutations = 1000)
betalactone_soil_MAG$mean_depth#0.051;P=0.51
betalactone_soil_MAG$P

hserlactone_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                        soil_MAG$hserlactone,
                                        min_fraction = 0.9,
                                        count_singletons = TRUE,
                                        singleton_resolution= 0,
                                        weighted = FALSE,
                                        Npermutations = 1000)
hserlactone_soil_MAG$mean_depth#0.044;P=0.684
hserlactone_soil_MAG$P

NRPS_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                 soil_MAG$NRPS,
                                 min_fraction = 0.9,
                                 count_singletons = TRUE,
                                 singleton_resolution= 0,
                                 weighted = FALSE,
                                 Npermutations = 1000)
NRPS_soil_MAG$mean_depth#0.057;P=0.536
NRPS_soil_MAG$P

RiPPs_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                  soil_MAG$RiPPs,
                                  min_fraction = 0.9,
                                  count_singletons = TRUE,
                                  singleton_resolution= 0,
                                  weighted = FALSE,
                                  Npermutations = 1000)
RiPPs_soil_MAG$mean_depth#0.061;P=0.243
RiPPs_soil_MAG$P

siderophore_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                        soil_MAG$siderophore,
                                        min_fraction = 0.9,
                                        count_singletons = TRUE,
                                        singleton_resolution= 0,
                                        weighted = FALSE,
                                        Npermutations = 1000)
siderophore_soil_MAG$mean_depth#0.037;P=0.923
siderophore_soil_MAG$P

Terpene_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                    soil_MAG$Terpene,
                                    min_fraction = 0.9,
                                    count_singletons = TRUE,
                                    singleton_resolution= 0,
                                    weighted = FALSE,
                                    Npermutations = 1000)
Terpene_soil_MAG$mean_depth#0.07;P=0.031
Terpene_soil_MAG$P

PKS_NRP_hybrids_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                            soil_MAG$PKS.NRP_Hybrids,
                                            min_fraction = 0.9,
                                            count_singletons = TRUE,
                                            singleton_resolution= 0,
                                            weighted = FALSE,
                                            Npermutations = 1000)
PKS_NRP_hybrids_soil_MAG$mean_depth#0.054;P=0.37
PKS_NRP_hybrids_soil_MAG$P

PKSI_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                 soil_MAG$PKSI,
                                 min_fraction = 0.9,
                                 count_singletons = TRUE,
                                 singleton_resolution= 0,
                                 weighted = FALSE,
                                 Npermutations = 1000)
PKSI_soil_MAG$mean_depth#0.059;P=0.152
PKSI_soil_MAG$P

PKSother_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                     soil_MAG$PKSother,
                                     min_fraction = 0.9,
                                     count_singletons = TRUE,
                                     singleton_resolution= 0,
                                     weighted = FALSE,
                                     Npermutations = 1000)
PKSother_soil_MAG$mean_depth#0.056;P=0.456
PKSother_soil_MAG$P

Others_soil_MAG<-consentrait_depth(soil_MAG_tree,
                                   soil_MAG$Others,
                                   min_fraction = 0.9,
                                   count_singletons = TRUE,
                                   singleton_resolution= 0,
                                   weighted = FALSE,
                                   Npermutations = 1000)
Others_soil_MAG$mean_depth#0.050;P=0.825
Others_soil_MAG$P

























#Codes to reproduce the Supplementary Figure 1B, 1D and Supplementary Figure 2
library(readxl)
library(dplyr)
all_soil_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/earthmicrobiome_bgc_data (1).xlsx",col_names = T,skip = 0,sheet = "soil_bgc")
all_plant_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/earthmicrobiome_bgc_data (1).xlsx",col_names = T,skip = 0,sheet = "majorplant")
dim(all_soil_mags)
dim(all_plant_mags)
plant_isolates<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/Network_Annotations_Full_all_isolates.xlsx",sheet = "Filtered",col_names = T,skip = 0)
soil_isolates<-read_excel("/Users/arijitmukherjee/Documents/soil_D_stat/soil_isolates.xlsx",sheet = "full_tab",col_names = T,skip = 0)

df<-data.frame(total=c(all_soil_mags$bgc,all_plant_mags$bgc,plant_isolates$total,soil_isolates$total),
               Data=c(rep("Soil MAGs",1300),rep("Plant MAGs",1453),rep("Plant isolates",2761),rep("Soil isolates",2533)))

df$Data<-factor(df$Data,levels=c("Plant isolates","Soil isolates","Plant MAGs","Soil MAGs"))



library(ggplot2)
library(viridis)
library(hrbrthemes)

#Fig.S2
#boxplot of total number of BGCs per genome from isolates and MAGs
box<-ggplot(df,aes(x=Data,y=total,fill=Data))+
  geom_jitter(color="black", size=0.2, alpha=0.2) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20)
  ) +
  ggtitle("") +
  ylab("Total BGC count")+
  xlab("")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12))
box
ggsave(
  "boxplot_all_data.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 400,
)
dev.off()



#======================#=============
#Fig. S1B
#Range for completeness in horizontal box plot; plot for completeness among three datasets
boxplot(all_soil_mags$completeness)
boxplot(all_plant_mags$Completeness)
completeness_df<-data.frame(Datasets=c(rep("Plant_isolates",2761),rep("Plant_MAGs",1453),rep("Soil_MAGs",1300)),
                            completeness=c(plant_isolates$Completeness,all_plant_mags$Completeness,all_soil_mags$completeness))
box_complete<-ggplot(completeness_df,aes(x=Datasets,y=completeness,fill=Datasets))+
  geom_boxplot(width=0.2)+
  coord_flip()+
  theme_classic()+
  ylab("Completeness (%)")+
  xlab("")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x=element_text(size=14))
box_complete
ggsave(
  "box_completeness.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 3,
  units = "in",
  dpi = 700,
)
dev.off()


#Fig. S1B=========##+==========
#Plotting contamination among all three datasets

contam_df<-data.frame(Datasets=c(rep("Plant_isolates",2761),rep("Plant_MAGs",1453),rep("Soil_MAGs",1300)),
                      contam=c(plant_isolates$Contamination,all_plant_mags$Contamination,all_soil_mags$contamination))

box_contamination<-ggplot(contam_df,aes(x=Datasets,y=contam,fill=Datasets))+
  geom_boxplot(width=0.2)+
  coord_flip()+
  theme_classic()+
  ylab("Contamination (%)")+
  xlab("")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x=element_text(size=14))
box_contamination
ggsave(
  "box_contam.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 3,
  units = "in",
  dpi = 700,
)
dev.off()

#Fig. S1B#==================================
#Genome size distribution among three datasets
GS_df<-data.frame(Datasets=c(rep("Plant_isolates",2761),rep("Plant_MAGs",1453),rep("Soil_MAGs",1300)),
                  Genome_size=c(plant_isolates$Genome_Size/10^6,all_plant_mags$`Genome_size (bp)`/10^6,all_soil_mags$genome_length/10^6))# Genome size has been converted into Mbp

box_GS<-ggplot(GS_df,aes(x=Datasets,y=Genome_size,fill=Datasets))+
  geom_boxplot(width=0.2)+
  coord_flip()+
  theme_classic()+
  ylab("Genome size (Mb)")+
  xlab("")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x=element_text(size=14))
box_GS
ggsave(
  "box_GS.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 3,
  units = "in",
  dpi = 700,
)
dev.off()

#Fig.S1D===========================
#Stack bar plot for three datasets
library(readxl)
library(dplyr)
plant_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/earthmicrobiome_bgc_data (1).xlsx",sheet = "majorplant",col_names = T,skip = 0)
soil_mags<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/earthmicrobiome_bgc_data (1).xlsx",sheet = "soil_bgc",col_names = T,skip = 0)
plant_mags%>%group_by(Phylum)%>%summarise(count=n())
soil_mags%>%group_by(Phylum)%>%summarise(count=n())
plant_isolates%>%group_by(Phylum)%>%summarise(count=n())
df<-data.frame(Datasets=c(rep("Plant MAGs",5),rep("Soil MAGs",5),rep("Plant isolates",5),rep("Soil isolates",5)),
               Phylum=c("Proteobacteria","Actinobacteria","Bacteroidetes","Firmicutes","Others"),
               Number=c(894,132,178,20,229,
                        446,209,269,59,317,
                        1800,325,86,548,2,
                        968,868,195,377,145))
df

library(ggplot2)
stack<-ggplot(df,aes(x=Datasets,y=Number,fill=Phylum))+
  geom_bar(stat = "identity",position = "fill")+
  theme_classic()+
  xlab("")+
  ylab("Proportion")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))

stack

ggsave(
  "stack_bar_new.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = "in",
  dpi = 400,
)
dev.off() 

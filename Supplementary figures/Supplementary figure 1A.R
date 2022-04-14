#code to reproduce Supplementary Figure 1A

library(tidyverse)
library(ggplot2)
library(readr)
library(maps)
library(viridis)
library(readxl)
data$longitude
df_MAGs_wm<-data%>%group_by(longitude,latitude)%>%summarise(count=n())%>%arrange(-count)
write.table(df_MAGs_wm,"df_MAGs_wm.tsv",sep = "\t")
data_mags_plants<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/Fig.1/new_report (2).xlsx",
                             sheet = "Plant_MAGs",col_names = T,skip = 0)

data_mags_soil<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/Fig.1/new_report (2).xlsx",
                           sheet = "Soil_MAGs",col_names = T,skip = 0)

df_MAGs_soil<-data_mags_soil%>%group_by(longitude,latitude)%>%summarise(count=n())%>%arrange(-count)
df_MAGs_soil
write.table(df_MAGs_soil,"df_MAGs_soil2.tsv",sep = "\t")

df_MAGs_plants<-data_mags_plants%>%group_by(longitude,latitude)%>%summarise(count=n())%>%arrange(-count)
df_MAGs_plants
write.table(df_MAGs_plants,"df_MAGs_plants2.tsv",sep = "\t")





plant_isolates<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/Network_Annotations_Full_all_isolates.xlsx",sheet = "Metadata",col_names = T,skip = 0)
plant_isolates<-plant_isolates%>%filter(Country!="NA")
dim(plant_isolates)
plant_df<-plant_isolates%>%group_by(Country)%>%summarise(count=n())%>%arrange(-count)
head(plant_df)
write.table(plant_df,"Plant_isolates_location.tsv",sep = "\t")

MAGs_df<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/Fig.1/new_report (2).xlsx",
                    sheet = "MAGs_df",col_names = T,skip = 0)
head(MAGs_df)
MAGs_df$Number<-as.factor(MAGs_df$Number)
MAGs_df$Dataset<-as.factor(MAGs_df$Dataset)

# Libraries
library(ggplot2)
library(dplyr)
world<-map_data("world")


library(maps)
wmags<-ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.8) +
  geom_point( data=MAGs_df, aes(x=longitude, y=latitude,size=Number,col=Dataset),alpha=0.6) +
  theme_void()+
  xlab("")+
  ylab("")
theme(axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12))
wmags
saveRDS(wmags,"World_MAGs.rds")
ggsave(
  "world_map_mags.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 4,
  units = "in",
  dpi = 700,
)
dev.off()


#World map isolates
library(readxl)
library(readr)
plant_isolates<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/Fig.1/new_report (2).xlsx",
                           sheet="PA_isolates",col_names = T,skip = 0)
head(plant_isolates)
library(dplyr)
library(ggplot2)
plant_isolates_wm_df<-plant_isolates%>%group_by(Latitude,Longitude)%>%summarise(count=n())%>%arrange(-count)
write.table(plant_isolates_wm_df,"plant_isolates_wm_df2.tsv",sep = "\t")

#Readin the formatted data
isolates_df<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/Fig.1/new_report (2).xlsx",
                        sheet = "isolates_location2",col_names = T,skip = 0)
isolates_df$Number<-as.factor(isolates_df$Number)

# Get the world polygon and extract UK
library(maps)
wisolates<-ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.8) +
  geom_point( data=isolates_df, aes(x=longitude, y=latitude,size=Number),color="steelblue",alpha=0.6) +
  theme_void()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
wisolates
saveRDS(wisolates,file = "World_map_isolates.rds")
ggsave(
  "world_map_isolates.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 4,
  units = "in",
  dpi = 700,
)
dev.off()








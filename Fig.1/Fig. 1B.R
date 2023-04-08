#Code to reproduce Fig. 1B
#import the data
library(readxl)
library(ggplot2)
library(dplyr)
library(MASS)
library(vegan)
data<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/Network_Annotations_Full_all_isolates.xlsx",
                 sheet = "drep_98_PcoA",col_names = T,skip = 0)

head(data)
genus_count<-data%>%group_by(Genus)%>%summarise(count=n())%>%arrange(-count)
genus_count

#from this data plot filter the top seven upto streptomyces genera

cols=c("#999999", "#006400", "#56B4E9","#999900","#3333FF","#000000","#FF9933")
data7genus<-data%>%filter(Genus%in%c("Pseudomonas_E","Rhizobium","Mesorhizobium",
                                     "Xanthomonas","Bacillus_A","Sphingomonas","Streptomyces"))
dim(data7genus)
data2<-data7genus[,9:19]
data2$genome<-data7genus$Genome_ID
names(data2)
data2<-as.data.frame(data2)
rownames(data2)<-data2$genome
data2_norm<-data2[,-12]/rowSums(data2[,-12])
dim(data2_norm)
abund_table<-as.matrix(data2_norm)
head(abund_table)
distance<-vegdist(abund_table,method = "bray")
sol_MDS<-cmdscale(distance,k=2,eig = T)

MDS.df<-data.frame(x=sol_MDS$points[,1],y=sol_MDS$points[,2],
                   Genus=as.factor(data7genus$Genus))
explainedvar1<-round(sol_MDS$eig[1]/sum(sol_MDS$eig),2)*100
explainedvar2<-round(sol_MDS$eig[2]/sum(sol_MDS$eig),2)*100
sum_var<-explainedvar1+explainedvar2
sum_var
explainedvar1
explainedvar2
adonis(distance~Genus,data = data7genus)

MDS.df
library(wesanderson)
pcoA_genus<-ggplot(MDS.df,aes(x=x,y=y,col=Genus))+
  geom_point(size=1,alpha=0.7)+
  theme_classic()+
  xlab("PCo1 (39%)")+
  ylab("PCo2 (19%)")+
  scale_color_manual(values=cols)+
  stat_ellipse(level = 0.8)+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x=element_text(size=14))+
  theme(legend.position = "top",legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.key.height = unit(0.8,'cm'))+
  theme(legend.key.width = unit(0.8,'cm'))+
  annotate(geom = "text",x=0.4,y=0.4,label=expression("R"^2*"=0.61; P=0.001"))
pcoA_genus
ggsave(
  "pcoa_final.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7.5,
  height = 7,
  units = "in",
  dpi = 400,
)
dev.off()

saveRDS(pcoA_genus,"PcoA_genus_final.rds")











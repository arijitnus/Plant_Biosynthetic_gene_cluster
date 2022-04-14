#Selecting appropriate model; we tested four different models based on count data
#We tested Poisson, NB, ZINB or ZIP model and tested for each of the BGC category 
#This was done for isolates dataset since the variation explained was higher in isolates dataset 
#This was done at the genus level since higher taxonomic resolution explained better variance as shown In Fig.3A
library(readxl)
library(dplyr)
library(ggplot2)
library(MASS)
library(pscl)
##
data<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/Network_Annotations_Full_all_isolates.xlsx",
                 sheet = "drep_98_PcoA",col_names = T,skip = 0)

head(data)
View(data)
dim(data)
data%>%group_by(Genus)%>%summarise(count=n())%>%arrange(-count)




top7genus<-c("Pseudomonas_E","Rhizobium","Mesorhizobium","Xanthomonas",
             "Bacillus_A","Sphingomonas","Streptomyces")
data7genus<-data%>%filter(Genus%in%top7genus)
dim(data7genus)
#fit NB models
#arly

data7genus
#arylpolyene
aryl_poisson<-glm(arylpolyene~Genus,data = data7genus)
aryl_nb<-glm.nb(arylpolyene~Genus,data = data7genus)
aryl_zinb<-zeroinfl(arylpolyene~Genus,data = data7genus,dist = "negbin")
aryl_zip<-zeroinfl(arylpolyene~Genus,data = data7genus,dist = "poisson")
AIC(aryl_poisson,aryl_nb,aryl_zip,aryl_zinb)

#betalactone
betalactone_poisson<-glm(betalactone~Genus,data = data7genus)
betalactone_nb<-glm.nb(betalactone~Genus,data = data7genus)
betalactone_zinb<-zeroinfl(betalactone~Genus,data = data7genus,dist = "negbin")
betalactone_zip<-zeroinfl(betalactone~Genus,data = data7genus,dist = "poisson")
AIC(betalactone_poisson,betalactone_nb,betalactone_zip,betalactone_zinb)

#hserlactone
hserlactone_poisson<-glm(hserlactone~Genus,data = data7genus)
hserlactone_nb<-glm.nb(hserlactone~Genus,data = data7genus)
hserlactone_zinb<-zeroinfl(hserlactone~Genus,data = data7genus,dist = "negbin")
hserlactone_zip<-zeroinfl(hserlactone~Genus,data = data7genus,dist = "poisson")
AIC(hserlactone_poisson,hserlactone_nb,hserlactone_zip,hserlactone_zinb)

#NRPS
NRPS_poisson<-glm(NRPS~Genus,data = data7genus)
NRPS_nb<-glm.nb(NRPS~Genus,data = data7genus)
NRPS_zinb<-zeroinfl(NRPS~Genus,data = data7genus,dist = "negbin")
NRPS_zip<-zeroinfl(NRPS~Genus,data = data7genus,dist = "poisson")
AIC(NRPS_poisson,NRPS_nb,NRPS_zip,NRPS_zinb)
data$PKS.NRP_Hybrids
#PKS-NRPE Hydbrids
PKS.NRP_Hybrids_poisson<-glm(PKS.NRP_Hybrids~Genus,data = data7genus)
PKS.NRP_Hybrids_nb<-glm.nb(PKS.NRP_Hybrids~Genus,data = data7genus)
PKS.NRP_Hybrids_zinb<-zeroinfl(PKS.NRP_Hybrids~Genus,data = data7genus,dist = "negbin")
PKS.NRP_Hybrids_zip<-zeroinfl(PKS.NRP_Hybrids~Genus,data = data7genus,dist = "poisson")
AIC(PKS.NRP_Hybrids_poisson,PKS.NRP_Hybrids_nb,PKS.NRP_Hybrids_zip,PKS.NRP_Hybrids_zinb)
#PKSI
PKSI_poisson<-glm(PKSI~Genus,data = data7genus)
PKSI_nb<-glm.nb(PKSI~Genus,data = data7genus)
PKSI_zinb<-zeroinfl(PKSI~Genus,data = data7genus,dist = "negbin")
PKSI_zip<-zeroinfl(PKSI~Genus,data = data7genus,dist = "poisson")
AIC(PKSI_poisson,PKSI_nb,PKSI_zip,PKSI_zinb)
#PKSother
PKSother_poisson<-glm(PKSother~Genus,data = data7genus)
PKSother_nb<-glm.nb(PKSother~Genus,data = data7genus)
PKSother_zinb<-zeroinfl(PKSother~Genus,data = data7genus,dist = "negbin")
PKSother_zip<-zeroinfl(PKSother~Genus,data = data7genus,dist = "poisson")
AIC(PKSother_poisson,PKSother_nb,PKSother_zip,PKSother_zinb)
#RiPPs
RiPPs_poisson<-glm(RiPPs~Genus,data = data7genus)
RiPPs_nb<-glm.nb(RiPPs~Genus,data = data7genus)
RiPPs_zinb<-zeroinfl(RiPPs~Genus,data = data7genus,dist = "negbin")
RiPPs_zip<-zeroinfl(RiPPs~Genus,data = data7genus,dist = "poisson")
AIC(RiPPs_poisson,RiPPs_nb,RiPPs_zip,RiPPs_zinb)
#Siderophore
siderophore_poisson<-glm(siderophore~Genus,data = data7genus)
siderophore_nb<-glm.nb(siderophore~Genus,data = data7genus)
siderophore_zinb<-zeroinfl(siderophore~Genus,data = data7genus,dist = "negbin")
siderophore_zip<-zeroinfl(siderophore~Genus,data = data7genus,dist = "poisson")
AIC(siderophore_poisson,siderophore_nb,siderophore_zip,siderophore_zinb)
data6genus$Terpene
#terpene
terpene_poisson<-glm(Terpene~Genus,data = data7genus)
terpene_nb<-glm.nb(Terpene~Genus,data = data7genus)
terpene_zinb<-zeroinfl(Terpene~Genus,data = data7genus,dist = "negbin")
terpene_zip<-zeroinfl(Terpene~Genus,data = data7genus,dist = "poisson")
AIC(terpene_poisson,terpene_nb,terpene_zip,terpene_zinb)
#other
dim(data6genus)
data6genus$Others
Others_poisson<-glm(Others~ngenera,data = data7genus)
Others_nb<-glm.nb(Others~ngenera,data = data7genus)
Others_zinb<-zeroinfl(Others~ngenera,data = data7genus,dist = "negbin")
Others_zip<-zeroinfl(Others~ngenera,data = data7genus,dist = "poisson")
AIC(Others_poisson,Others_nb,Others_zip,Others_zinb)


summary(aryl_nb)
summary(betalactone_nb)
summary(hserlactone_nb)
summary(NRPS_nb)
summary(PKS.NRP_Hybrids_nb)
summary(PKSother_nb)
summary(PKSI_nb)
summary(RiPPs_nb)
summary(siderophore_nb)
summary(terpene_nb)

models_nb_diff<-list(aryl_nb,betalactone_nb,hserlactone_nb,NRPS_nb,
                     PKS.NRP_Hybrids_nb,PKSI_nb,PKSother_nb,RiPPs_nb,
                     siderophore_nb,terpene_nb)


summary(aryl_nb)

saveRDS(models_nb_diff,"diff_enrichment_nb_models.rds")

summary(aryl_nb)$coefficients[,1:4]

df=vector("list",10)
for (i in 1:10) {
  df[[i]]=summary(models_nb_diff[[i]])$coefficients[,1:4]
}






dim(df)
#forest plot for Genus
library(readxl)
library(ggplot2)
cols=c("#999999", "#006400", "#56B4E9","#FF6347","#3333FF","#000000","#FF9933")
fp_genus_nb<-read_excel("/Users/arijitmukherjee/Documents/Phytobiome/BGC_isolates/New_final_data/MAGs/Earth_MAGs/latest_ANI_analysis/Supplementary figures/new_report (4).xlsx",
                          skip = 0,col_names = T,sheet = "input_FP")

p = ggplot(data=fp_genus_nb,
           aes(x = Genus,y = Estimate, ymin = lower, ymax = upper ))+
  geom_pointrange(aes(col=Genus))+
  scale_color_manual(values = cols)+
  #geom_hline(aes(fill=Phylum),yintercept =1, linetype=2)+
  xlab('')+ ylab(expression("Log"[10]*" mean count per genome at 95% confidence interval"))+
  geom_errorbar(aes(ymin=lower, ymax=upper,col=Genus),width=0.3,cex=0.6)+ 
  facet_wrap(~BGC,strip.position="left",nrow=4,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_text(size = 14),
        axis.ticks.y=element_text(size = 14),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))
q<-p+theme_bw()+
  coord_flip()
q
saveRDS(q,"Forest_plot_genus.rds")
ggsave(
  "fores_genus_final.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 350,
)
dev.off()

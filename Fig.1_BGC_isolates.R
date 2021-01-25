#Fig.1 BGC manuscript
#Reading metadata file
rm(list=ls())
metadata<-read.table("isolate_metadata.tsv",header = T,sep = "\t")
head(metadata)
str(metadata)
dim(metadata)#Total 4648 genomes and 14 columns 
#First chech whether there is any genome with less than 10 percent contamination
#or <50% completedness
any(metadata$Completeness<50)
any(metadata$Contamination>10)#ok proceed
# Check the distribution of genome size,contamination and completeness
range(metadata$GC,na.rm = T)#30.4 to 41.6
range(metadata$Genome_size..bp./10^6,na.rm = T)#0.159 to 4.08 Mb
range(metadata$Completeness)#54.6 to 100 percent
range(metadata$Contamination)#0 to 9.82
#plot the boxplots 
#create a vector of GC removing NA containing rows
ggplot(metadata, aes(x=GC)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  theme_bw()+
  xlab("GC content (%)")+
  ylab("Density")
#This gives the plot for distribution of GC content
#Plotting the genome size distribution
metadata$gmbp<-metadata$Genome_size..bp./10^6
ggplot(metadata,aes(x=gmbp))+
  geom_histogram(binwidth = 0.1,col="blue",fill="lightblue")+
  theme_bw()+
  xlab("Genome size (Mbp)")+
  ylab("No. of isolates")
#plot for completeness and contamination to be done later
length(metadata$Completeness[metadata$Completeness<95])# 76 genomes have completeness less than 95
length(metadata$Completeness[metadata$Completeness>=95])#4572
length(metadata$Completeness[metadata$Completeness>=96])#4550
length(metadata$Completeness[metadata$Completeness>=97])#4495
length(metadata$Completeness[metadata$Completeness>=98])#4337
length(metadata$Completeness[metadata$Completeness>=99])#3374
#76 genomes <95, 1198 genomes 95-99 and 3374 genomes>=99
#1.63% <95, 25.77% genomes 95-99 and 72.6% >=99
#Make a pie chart
# Simple Pie Chart
slices <- c(76,1198,3374)
lbls <- c("<95% completeness","95-99% completeness",">99% completeness")
pie(slices, labels = lbls)

#observing many isolates don't have information on country and plant type



#==================================
#==================================
# now for the world map
# load packages
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
# load data
world <- ne_countries(scale = "medium", returnclass = "sf")
# generic world map
ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$admin)), " countries)"))+
  geom_point(sovereignt)


#we want to add data points of countries from where we obtained the
#isolates grouped by country and plant family
length(unique(metadata$Host_Scientific_name))#42 plants
length(unique(metadata$Country))# 84 countries covered
#Group_by countries and get the no of isolates
iso_count_country<-metadata%>%group_by(Country)%>%summarise(isolates_count=n())
#plot the countries with the number of isolates
iso_count_country<-iso_count_country[2:84,]
iso_count_country
library(ggmap)
library(maptools)
library(maps)
library(countrycode)
visited <- iso_count_country$Country
geocode("Brazil")
ll.visited <- geocode(visited)
visit.x <- ll.visited$lon
visit.y <- ll.visited$lat
iso_count_country$ISO_3<-countrycode(iso_count_country$Country, origin = 'country.name', destination = 'iso3c')
iso_count_country_f<-iso_count_country[-c(24,81),]
#has No NA
# create data frame with iso3 country codes and number of visits
datatable(iso_count_country_f, rownames = FALSE, options = list(pageLength = 5, scrollX=T), filter = "none")
install.packages(c("RgoogleMaps", "ggmap", "mapproj", "sf",
                   "dplyr", "OpenStreetMap", "devtools", "DT"))

#create a vector for the condition
#USA itself has 1097 isolates
library(arules)
discretize(iso_count_country_f$Country, method = "", breaks = 3, 
           labels = NULL, include.lowest = TRUE, right = FALSE, dig.lab = 3,
           ordered_result = FALSE, infinity = FALSE, onlycuts = FALSE, 
           categories, ...)






library(RgoogleMaps)
library(ggmap)
library(mapproj)
library(sf)
library(DT)
library(devtools)
library(OpenStreetMap)
library(rworldmap)
library(readxl)
iso_country_df<-read_excel("~/Documents/Phytobiome/BGC_isolates/Final_metadata_figures/Country_data_fig1.xlsx",
                           skip = 0,col_names = T,sheet = "Sheet1")
# combine data frame with map
visitedMap <- joinCountryData2Map(iso_country_df, 
                                  joinCode = "ISO3",
                                  nameJoinColumn = "ISO_3")

# def. map parameters, e.g. def. colors
mapParams <- mapCountryData(visitedMap, 
                            nameColumnToPlot="count",
                            oceanCol = "azure2",
                            catMethod = "categorical",
                            missingCountryCol = gray(.8),
                            colourPalette = c("coral",
                                              "coral2",
                                              "coral3", "orangered" 
                                              ),
                            addLegend = F,
                            mapTitle = "",
                            aspect = 1,
                            border = NA)
do.call(addMapLegendBoxes, c(mapParams,
                             x = 'top',
                             title = "No. of isolates",
                             horiz = TRUE,
                             bg = "transparent",
                             bty = "n"))

#Create a variable where you store the number as factors
#=======================================
#==============plot the number of isolates from each plant family
library(readxl)
library(dplyr)
library(ggplot2)
plant_isolates<-read_excel("~/Documents/Phytobiome/BGC_isolates/Final_metadata_figures/Isolates_metadata.xlsx",
                           skip=0,sheet="Sheet1",col_names=T)
head(plant_isolates)
plant_isolates<-na.omit(plant_isolates)
dim(plant_isolates)
plant_isolates<-filter(plant_isolates,Completeness>=50&Contamination<10)
dim(plant_isolates)
sc_name<-plant_isolates%>%group_by(Host_Scientific_name)%>%summarise(count=n())
top_10_plants<-sc_name%>%arrange(desc(count))
top_10_plants
top_10_plants<-top_10_plants[1:10,]
top_10_plants
barplot(top_10_plants$count,col ="Steelblue")













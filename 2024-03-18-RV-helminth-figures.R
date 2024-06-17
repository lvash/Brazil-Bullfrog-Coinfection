# The creation of these figures corresponds with the article in International Journal of Parasitology: Parasites and Wildlife:
# Ranavirus and helminth parasite co-infection in invasive American bullfrogs in the Atlantic Forest, Brazil
# Authors: Lauren V. Ash, Karla Magalhães Campião, Cauê Pinheiro Teixeira, Nicholas J. Gotelli

# LVA
# 18 March 2024

# --------------------------------------

### Libraries
library(tidyverse)
library(colorspace)
library(colorblindr)
library(cowplot)
library(viridis)
library(ggeffects)
library(interactions)
library(ggiraphExtra)
library(jtools)
library(EcoSimR) #remotes::install_github("GotelliLab/EcoSimR")

# Mapping libraries
library(terra)
library(tidyterra)
library(sf)
library(raster)
library(spData)
library(tmap)    
library(leaflet)
library(ggmap)
require(maps)
library(scatterpie)

# --------------------------------------
### Figure 1: map of Brasil

## Read in data
data <- read.csv("Data/parasite-RV-final-data.csv", header=T, comment.char = "#")
data$RanavirusFac <- as.factor(data$Ranavirus)
data$LengthScaled <- LengthScaled <- scale(data$Length, center=T, scale=T)

## Read in site data
sites <- read.csv("Data/Brazil-Sites.csv", header=T, comment.char = "#")

sites <- sites %>%
  mutate(RelAbund = TotAbund/N, RV_NegN = N - RV_PosN, RV_load_log = log(RV_load + 1))

sites$InfStatus <- factor(ifelse(sites$RV_PosN > 0, 1, 0))

# API for google maps
# https://console.cloud.google.com/google/maps-apis/credentials
register_google(key="AIzaSyCzYNgRRB3sQJNDOwLuP7uFrkLumDYhur0")

brazilFar<-ggmap::get_googlemap(center=c(lon=-58, lat = -15), zoom=4, maptype="satellite") 

# Download country borders
dsn <- "/Bullfrogs_Brasil/gadm41_BRA_shp"
lay <- "gadm41_BRA_1"
brazilBorders1 <- read_sf(dsn = dsn, layer = lay)

# The borders of only states that were sampled in
sample.states <- brazilBorders1[brazilBorders1$NAME_1 == "Paraná" | brazilBorders1$NAME_1 == "Santa Catarina" | brazilBorders1$NAME_1 == "São Paulo",]

sampleStateBorders <- st_transform(sample.states, CRS("+proj=longlat +datum=WGS84"))
sampleStateBorders <- fortify(sampleStateBorders)


brazilFarRast<- rast(brazilFar) #convert to raster
bordsAll <- vect(brazilBorders1)
br.far.only <- terra::mask(brazilFarRast, bordsAll)

# RESOLVE Ecoregions and Biomes shapefile
# retrieved here: https://hub.arcgis.com/datasets/esri::resolve-ecoregions-and-biomes/about
biomes <- read_sf("Biomes/Ecoregions2017.shp")
TropFor_biome <- biomes[biomes$BIOME_NAME=="Tropical & Subtropical Moist Broadleaf Forests",]
rm(biomes)

# tropical forest biome for only brazil
template <- rast(vect(TropFor_biome), res=0.005)
biomesRast <- rasterize(vect(TropFor_biome), template) # takes a minute
biomes.br.only <- terra::mask(biomesRast, bordsAll) # also takes a minute
biomes.sf <- terra::as.polygons(biomes.br.only)

# Mapping whole country (Figure 1) code

plot_Far <- ggplot() +
  geom_spatraster_rgb(data = br.far.only) +
  geom_spatvector(data = bordsAll, color = "black", size = 8, na.rm=TRUE, show.legend = F) +
  geom_sf(data=biomes.sf, fill = "seagreen4", alpha=0.4, col=alpha("lightgray",0), linetype="11", inherit.aes = FALSE) + 
  coord_sf(xlim = c(-74.9, -30.9), ylim = c(-33.9, 7.9), expand = FALSE) +
  geom_sf(data=sampleStateBorders, fill=NA, col="gray10", linewidth=1) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="", y="") +
  theme(axis.text=element_text(size=18))
plot_Far

#ggsave(filename="BrazilFig1A.png", width = 8, height=6, units="in")

# --------------------------------------
# Figure 1: Close-up map of states where sampling occurred with helminth and Ranavirus prevalence pie graphs

brazilClose <-ggmap::get_googlemap(center=c(lon=-49.4, lat = -25.52), zoom=7, maptype="satellite") 

brazilrast<-rast(brazilClose) # convert to raster
bords <- vect(sample.states) # borders of sampling states
br.only <- terra::mask(brazilrast, bords) # clip to bounds 

### Base map
plot_ggplot <- ggplot() +
  geom_spatraster_rgb(data = br.only) +
  geom_spatvector(data = bords, color = "black", size = 8, na.rm=TRUE, show.legend = F) +
  geom_sf(data=TropFor_biome, fill = "seagreen4", alpha=0.4, col=alpha("lightgray",0), linetype="11", inherit.aes = FALSE) + 
  coord_sf(xlim = c(-52.9, -45.9), ylim = c(-28.7, -22.3), expand = FALSE) +
  geom_point(data = sites, aes(Longitude, Latitude, group=InfStatus, fill=InfStatus), color="black",shape=21, size = 4, alpha=0.75) +
  scale_fill_manual(values=c("white","darkred"),name="Ranavirus status", labels=c("Uninfected", "Infected")) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key=element_rect(fill="white")) +
  labs(x="", y="") +
  theme(axis.text=element_text(size=18))
plot_ggplot

### Helminth prevalence pie charts

cols <- c("#440154FF","#35B779FF","#31688EFF","#FDE725FF")

ggplot()+ geom_scatterpie(aes(x=LongMap, y=LatMap, group=Locality, r=RelAbund/85), data=sites, legend_name = "Helminth.Taxa", cols=c("Acantho","Nema","Penta","Platy")) + 
  xlab("") + ylab("")  + 
  scale_fill_manual(values=cols, name="Helminth taxa", labels=c ("Acanthocephala", "Nematoda", "Pentastomida", "Platyhelminthes")) +
  theme(panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"))

#ggsave("Helminth-transparent.png", bg="transparent")

### Ranavirus prevalence pie charts

cols2<-c("darkred","white")

ggplot() + 
  geom_scatterpie(aes(x=Longitude, y=Latitude, group=Locality, r=RV_load_log/30), data=sites, legend_name = "Ranavirus", cols=c("RV_PosN", "RV_NegN")) + 
  xlab("") + ylab("") +
  scale_fill_manual(values=cols2, name="Ranavirus status", 
                    labels=c("Infected", "Uninfected")) + 
  theme(legend.position="left", 
  panel.background = element_rect(fill = "transparent",  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  plot.background = element_rect(fill = "transparent",
  colour = NA_character_), 
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent"))

## save as transparent
# ggsave(filename = "RV_pie-charts.png", bg = "transparent")

# --------------------------------------
### Figure 2: EcoSimR (for full analysis see coinfection analysis script)

### Create matrix
parasiteMat1 <- data %>%
  mutate(AcanthoPres=ifelse(Acantho>0, 1, 0),
         NemaPres= ifelse(Nema>0,1,0),
         TremaPres= ifelse(Trema>0,1,0),
         CestPres= ifelse(Cest>0,1,0),
         PentaPres=ifelse(Pentastomida_larva>0, 1, 0)) %>%
  dplyr::select(Ranavirus, AcanthoPres, NemaPres, TremaPres, CestPres, PentaPres)

parasiteMat<-t(parasiteMat1)

# Cleaning: remove NA in col 5 of Length data (and matrix)
pM <- parasiteMat[,-5]
len <- data$Length[-5]

# sim10 algorithm with bullfrogs weighted by size 
myModsim10 <- cooc_null_model(pM, algo = "sim10", suppressProg=TRUE, algoOpts = list(colWeights=len))
summary(myModsim10)
plot(myModsim10, type="cooc")
plot(myModsim10, type="hist")

# --------------------------------------
### Figure 3a
### Negative Relationship between Ranavirus load and nematode abundance in infected individuals

data1 <- data[complete.cases(data),] # no NAs 
data1$Log_Nema <- log10(ifelse(data1$Nema == 0, 1, data1$Nema)) # calculate (log) of abundance (not zeros) to normalize
data1$RV <- round(data1$RV_viral.load) # Poisson needs integers
data1$RV_viral.load <- log(data1$RV) # easier to show relationship when viral load is log-transformed and given Poisson distribution

# significant relationship with poisson distribution too
summary(glm(RV ~ Log_Nema, family= "poisson", data1)) # p<2e-16 ***

### Figure 3a
fig3a <- ggplot(data1, aes(x= Log_Nema, y= round(RV_viral.load))) + geom_point() + ylab("Log(Ranavirus viral load)") + xlab("Log(Nematoda abundance)") +
  ylim(0,12) + 
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  xlim(0,2) +
  stat_smooth(method="glm",color="black", linewidth=0.75, method.args = list(family = "poisson"), 
              se = FALSE) + theme_classic(base_size = 16)
fig3a 

#ggsave(filename="Fig3a.png",fig3a, width=5.5, height=5, units= "in")

# --------------------------------------
### Figure 3b
### Effect of Ranavirus infection status on the relationship between bullfrog size (length) and total helminth abundance (Figure 3b) 

### Data subset so only frogs that had the opportunity to be infected are compared (bullfrogs from infected sites)

posTowns<- data %>%
  group_by(Town) %>%
  dplyr::summarize(Pos=sum(Ranavirus))%>%
  filter(Pos>0)

dataPlot <- data %>%
  filter(Town %in% posTowns$Town)

# Total abundance

p <- ggplot(dataPlot, aes(x=Length, y=log(Total_abund+1), group=RanavirusFac, col=RanavirusFac)) + geom_point(size=2.5, pch=21, aes(fill=RanavirusFac), color="black") + theme_classic() + 
  geom_smooth(method = "lm", se = F) +
  ylab("Log(Total macroparasite abundance + 1)") + xlab("Length (cm)") +
  scale_fill_manual("Ranavirus infection status", values=viridis(2)) +
  scale_color_manual(name = "Ranavirus infection status", values=viridis(2)) + theme(legend.position = c(0.3,0.8), axis.text=element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14), legend.title = element_text(size=16)) 
p

gg_des <- edit_colors(p, desaturate)
Fig3b <-cowplot::ggdraw(gg_des)
Fig3b
#ggsave("Fig3b.png", width=7, height=5, units="in")

#-----------------------------------------
### Richness (Supplementary Figure): Effect of Ranavirus infection status on the relationship between bullfrog size (length) and helminth richness

rich <- ggplot(dataPlot, aes(x=Length, y= Total_rich, col = RanavirusFac)) + 
  geom_point(size = 2) + 
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"))

richPlot <- rich + geom_point(size=2.5, pch=21, aes(fill=RanavirusFac), color="black") + theme_classic() + 
  ylab("Macroparasite taxa richness") + xlab("Length (cm)") +
  scale_fill_manual("Ranavirus infection status", values=viridis(2)) +
  scale_color_manual(name = "Ranavirus infection status", values=viridis(2)) + theme(legend.position = c(0.3,0.8), axis.text=element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14), legend.title = element_text(size=16)) 

rich_des <- edit_colors(richPlot, desaturate)
FigS2<-cowplot::ggdraw(rich_des)
FigS2
#ggsave("FigS2.png", width=7, height=5, units="in")


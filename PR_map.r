
library(malariaAtlas)
library(ggplot2)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(maps)
library(tibble)

### How to download empty map/districts ###
# how to download single country (example = Kenya)


prad0 <- read_sf("/Users/lsh1513859/OneDrive - London School of Hygiene and Tropical Medicine/vector_genomics2/PuertoRico/pri_adm_2019_shp/pri_admbnda_adm0_2019.shp")

ggplot() +
  geom_sf(data = world) +
  geom_sf(data = prad0, fill = "firebrick") +
  coord_sf(xlim = c(-40, -100), ylim = c(0, 50), expand = FALSE) +
  theme_classic()

plot1 <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = prad0, fill = "#40476D", alpha = 0.9) +
  coord_sf(xlim = c(-60, -80), ylim = c(10, 25), expand = FALSE) +
  theme_classic()

prad1 <- read_sf("/Users/lsh1513859/OneDrive - London School of Hygiene and Tropical Medicine/vector_genomics2/PuertoRico/pri_adm_2019_shp/pri_admbnda_adm1_2019.shp")

prad1 <-   st_crop(prad1, xmin = -65.2, xmax = -67.3, ymin = 17, ymax = 19)
plot(prad1)
plot2 <- ggplot() +
  geom_sf(data = prad1, fill = "white") +
  geom_sf(data = prad1[prad1$ADM1_ES == "San Juan",], fill = "#40476D", alpha = 0.9) +
  geom_sf(data = prad1[prad1$ADM1_ES == "Ponce",], fill = "#40476D", alpha = 0.9) +
  #geom_sf(data = prad1[prad1$ADM1_ES == "Culebra",], fill = "#40476D") +
  geom_sf(data = prad1[prad1$ADM1_ES == "Bayamón",], fill = "#40476D", alpha = 0.9) +
  geom_sf(data = prad1[prad1$ADM1_ES == "Dorado",], fill = "#40476D", alpha = 0.9) +
  geom_sf(data = prad1[prad1$ADM1_ES == "Guánica",], fill = "#40476D", alpha = 0.9) +
  geom_sf(data = prad1[prad1$ADM1_ES == "Culebra",], fill = "#40476D", alpha = 0.9) +
  #ggtitle("Puerto Rico mosquito collection regions") +
  theme_classic() 

library(patchwork)

plot_combo <- plot1+ plot2 +plot_layout(nrow = 2)
plot_combo
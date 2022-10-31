### Map Notes ###
library(tidyverse)
library(knitr)
library(sf)
library(raster)
library(dplyr)
library(spData)

library(tmap)
library(ggplot2)

tm_shape(nz) + tm_polygons() + tm_grid()

eastcoast=st_bbox(c(xmin = -85, xmax = -65,
                    ymin = 30, ymax = 45), 
                  crs = st_crs(us_states)) %>% st_as_sfc()
tm_shape(us_states, bbox = eastcoast) + tm_graticules() + tm_polygons()

#Add points
spisloc <- read_csv("/Users/hannah/gitHub/Spisula/Work/map_spisula_input.csv")
#Format of this csv:



spisloc %>%
  ggplot() +
  geom_point(aes(x=Lat, y=Lon))
spissites <- st_as_sf(spisloc, coords = c("Lon","Lat"))
#making the coordinate reference system the same
st_crs(spissites) <- st_crs(us_states)

tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_grid(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 0.9, col = "Species", alpha = 0.5,  palette = "Set2") 

#Adding continental shelf contours
#USGS FILE FROM:
#https://www.sciencebase.gov/catalog/item/4fb69a6be4b03ad19d64b4c2
bary_iso <- st_read("/Users/hannah/gitHub/Spisula/Work/shape/bathymetry_l_v2.shp")
tm_shape(bary_iso) + tm_lines()
#Make specific to USA
bary_us <- bary_iso[-3171,]
tm_shape(us_states) + tm_polygons() + tm_shape(bary_us) + tm_lines()

#Plot all together
tm_shape(us_states, bbox=eastcoast) + tm_fill(alpha = 0) +
  tm_shape(bary_us, bbox=eastcoast) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states, bbox=eastcoast) +
  tm_fill(col = "grey86") + tm_borders(col="grey30") + tm_layout(bg.color = "grey99") +
  tm_graticules(col = "white", lwd = 0.7) +
  tm_shape(spissites) +  
  tm_symbols(size = 0.2, jitter=0.05,col = "Species",  palette = "Set2", border.col = "grey30")




#Plot for Paper
eastcoast=st_bbox(c(xmin = -85, xmax = -65,
                    ymin = 30, ymax = 45)) %>% st_as_sfc()
SNE=st_bbox(c(xmin = -76, xmax = -65,
              ymin = 37, ymax = 43)) %>% st_as_sfc()
GAonly=st_bbox(c(xmin = -83, xmax = -79,
                 ymin = 29, ymax = 33)) %>% st_as_sfc()
#Set the new pallete 
pal2022 <- c("#98D494","#9B70F8")

#Add more SsLI sites
#And the 1999 samples
spisloc2 <- read_csv("/Users/hannah/gitHub/Spisula/Work/map_spisula_input_mar222.csv")
spisloc2 %>%
  ggplot() +
  geom_point(aes(x=Lat, y=Lon))
#my data
spissites2 <- st_as_sf(spisloc2, coords = c("Lon","Lat")) #LON FIRST - if confused, view whole world
#making the coordinate reference system the same
st_crs(spissites2) <- st_crs(us_states)

tm_shape(us_states, bbox=SNE) + tm_fill(alpha = 0) +
  tm_shape(bary_us, bbox=SNE) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states, bbox=SNE) +
  tm_fill(col = "grey86") + tm_borders(col="grey30") + tm_layout(bg.color = "grey99") +
  tm_graticules(col = "white", lwd = 0.7) +
  tm_shape(spissites2) +  
  tm_symbols(size = 0.4, jitter=0.01, col = "Species",  palette = pal2022, border.col = "grey30") +
  tm_layout(scale=1.3)

SNE=st_bbox(c(xmin = -76, xmax = -65,
              ymin = 36.75, ymax = 42.75)) %>% st_as_sfc()
GAonly=st_bbox(c(xmin = -83, xmax = -79,
                 ymin = 29, ymax = 33)) %>% st_as_sfc()

tm_shape(us_states, bbox=SNE) + tm_fill(alpha = 0) +
  tm_shape(bary_us, bbox=SNE) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states, bbox=SNE) +
  tm_fill(col = "grey90") + 
  tm_borders(col="grey10", lwd=0.7) + tm_layout(bg.color = "grey99") +
  tm_graticules(col = "grey99", lwd = 1, labels.size = 0.8) +
  tm_shape(bary_us, bbox=SNE) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(spissites2) +  
  tm_symbols(size = 0.4, jitter=0.01, col = "Species",  palette = pal2022, border.col = "grey30") +
  tm_layout(scale=2.5)
#2000 width

tm_shape(us_states, bbox=GAonly) + tm_fill(alpha = 0) +
  tm_shape(bary_us, bbox=GAonly) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states, bbox=GAonly) +
  tm_fill(col = "grey90") + 
  tm_borders(col="grey10", lwd=0.7) + tm_layout(bg.color = "grey99") +
  tm_graticules(col = "grey99", lwd = 1) +
  tm_shape(bary_us, bbox=GAonly) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(spissites2, legend.show = FALSE) +  
  tm_symbols(size = 1, jitter=0.01, col = "Species",  palette = pal2022, border.col = "grey30") +
  tm_layout(scale=2.5, legend.show = FALSE)
#1500 width


tm_shape(us_states, bbox=SNE) + tm_fill(alpha =0) +
  tm_shape(us_states, bbox=SNE) +
  #tm_fill(col = "grey86") + 
  tm_borders(col="grey10", lwd=0.7) + tm_layout(bg.color = "grey99") +
  tm_layout(scale=1.3) +
  tm_view(alpha=1)


library(usmap)
library(ggplot2)
p <- plot_usmap(include = c("NY", "NJ", "MA", "VA","MD","DE","RI","CT","DC","PA","NH","VT","ME"),fill=NA)
ggsave(p, filename = "testmapagain.png",  bg = "transparent")







#Note: I might be able to make these nicer by using a finerscale polymap now that I don't need the whole US at once?
us_states2163 <- st_transform(us_states, 2163)
tm_shape(us_states2163) + tm_polygons() 

bary_iso <- st_read("/Users/hannah/gitHub/Spisula/Work/shape/bathymetry_l_v2.shp")
tm_shape(bary_iso) + tm_lines()
#Make specific to USA
bary_us <- bary_iso[-3171,]
bary_us_2 <- bary_iso[-2163,]
tm_shape(bary_us_2 ) + tm_lines()

tm_shape(us_states2163) + tm_polygons() + tm_shape(bary_us_2) + tm_lines() +  tm_shape(spissites2) +  
  tm_symbols(size = 0.4, jitter=0.01, col = "Species",  palette = pal2022, border.col = "grey30") +
  tm_layout(scale=1.3)

tm_shape(us_states2163, bbox=GAonly) + tm_polygons() + tm_shape(bary_us_2, bbox=GAonly) + tm_lines()
#Cant get it to accept bounded boxes

tm_shape(us_states2163, bbox=SNE) + tm_fill(alpha = 0) +
  tm_shape(bary_us, bbox=SNE) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states2163, bbox=SNE) +
  tm_fill(col = "grey86") + tm_borders(col="grey30") + tm_layout(bg.color = "grey99") +
  tm_graticules(col = "white", lwd = 0.7) +
  tm_shape(spissites2) +  
  tm_symbols(size = 0.4, jitter=0.01, col = "Species",  palette = pal2022, border.col = "grey30") +
  tm_layout(scale=1.3)

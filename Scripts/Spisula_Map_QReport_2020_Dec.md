Map Spisula Sites
================

## Intro and Practice

The code used here for reference is from Harmony, Honggang
(<https://github.com/hzz0024/BioinfoTools>), and
<https://geocompr.robinlovelace.net/adv-map.html>.

My goal is to make a map of the spisula sample sites from 2012 to 2020,
with number of samples coded by the size of the circle and subspecies
coded by color.

``` r
tm_shape(nz) +
    tm_fill() 
map_nz1 = tm_shape(nz) + tm_polygons()
#add elevation
#map_nz1 = map_nz +
#  tm_shape(nz_elev) + tm_raster(alpha = 0.7) #need to define nz_elevation
#representing nz territorial waters
nz_water = st_union(nz) %>% st_buffer(22200) %>% 
  st_cast(to = "LINESTRING")
map_nz2 = map_nz1 +
  tm_shape(nz_water) + tm_lines()
#add peaks
map_nz3 = map_nz2 +
  tm_shape(nz_height) + tm_dots()  #but why don't I need to define nz_height
tmap_arrange(map_nz1, map_nz2, map_nz3)
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/regional-1.png)<!-- -->![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/regional-2.png)<!-- -->![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/regional-3.png)<!-- -->![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/regional-4.png)<!-- -->![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/regional-5.png)<!-- -->

``` r
library(sf)
library(raster)
library(dplyr)
library(spData)
library(tmap)
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ tibble  2.1.3     ✓ purrr   0.3.3
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x tidyr::extract() masks raster::extract()
    ## x dplyr::filter()  masks stats::filter()
    ## x dplyr::lag()     masks stats::lag()
    ## x dplyr::select()  masks raster::select()

``` r
#Show full map
tm_shape(us_states) + tm_polygons() + tm_graticules() 
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/formatit-1.png)<!-- -->

``` r
#create my bounds with st_bbox, piped into st as sfc
eastcoast=st_bbox(c(xmin = -85, xmax = -65,
                      ymin = 30, ymax = 45)) %>% st_as_sfc()
#use those bounds under bbox variable within tm_Shape
tm_shape(us_states, bbox = eastcoast) + tm_polygons() + tm_graticules() 
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/formatit-2.png)<!-- -->

``` r
#tm_polygons = tm_fill + tm_borders
#color and formatting
tm_shape(us_states, bbox = eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white")
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/formatit-3.png)<!-- -->
\#\# Importting
Data

``` r
spisloc <- read_csv("/Users/hannah/gitHub/Spisula/Work/map_spisula_input.csv")

spisloc %>%
  ggplot() +
  geom_point(aes(x=Lat, y=Lon))
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/importdata-1.png)<!-- -->

``` r
#Add on to map
#example from: https://www.r-bloggers.com/2019/10/make-a-map-of-your-study-site-with-tmap/
sitescoral <- read.csv(url("https://raw.githubusercontent.com/cbrown5/BenthicLatent/master/data-raw/JuvUVCSites_with_ReefTypes_16Jun2016.csv"))
sites <- st_as_sf(sitescoral, coords = c("coordx", "coordy"))
#making the coordinate reference system the same
#land <- st_read(dsn ="LandPoly")
#st_crs(sites) <- st_crs(land)

  
#my data
spissites <- st_as_sf(spisloc, coords = c("Lon","Lat")) #LON FIRST - if confused, view whole world
#making the coordinate reference system the same
st_crs(spissites) <- st_crs(us_states)


library(RColorBrewer)
#plotting
tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_grid(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 0.9, col = "Species", alpha = 0.5,  palette = "Set2") 
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/importdata-2.png)<!-- -->

``` r
#For the future: to reverse a pallette, easy as: palette = rev("Set2")
```

### Adjustments

Possible iterations of the graph: +Where color or shape indicate species
+Whiter background +Where size is number of samples taken at the site
+With nantuket similis offset

``` r
eastcoast=st_bbox(c(xmin = -85, xmax = -65,
                      ymin = 30, ymax = 45)) %>% st_as_sfc()
#Base
map_1 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 1, col = "Species", alpha = 0.5,  palette = "Set2", border.alpha = 0.2)
#Size = number of samples
map_2 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = "Samples", col = "Species", alpha = 0.5,  palette = "Set2", border.alpha = NA)
#Also changed border alpha to default
#All are too small!!

#By shapes
map_3 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 1, shape = "Species",  alpha = 0.7,  col= "Samples", palette = "Reds")

#Try to change size and jitter (for sites with both)
map_4 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = "Samples", col = "Species", alpha = 0.7,  palette = "Set2",
             scale = 2.5, jitter = 0.01) 
#default icon.scale=3 ,  icon.scale = 9 has no change
#default legend.max.symbol.size = 0.8 , on its own legend.max.symbol.size = 3 has no effect
#THIS ONE: default scale=1, changing this makes them all bigger

#Circles not filled in
map_5 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = "Samples", alpha = 1, col = "Species",shape=1, border.lwd = 2, palette = "Set2",
             scale = 2.5, jitter = 0.01) 
#Do so by choosing a different shape (hollow circle) rather than trying to set border color to variable

#map_1 at larger scale
map_6 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey90") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 1, jitter= 0.01 ,col = "Species",  palette = "Set2",alpha = 0.6,border.lwd = 0.7,  border.alpha = 0.5, scale = 1.5)

map_7 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey94") +
  tm_graticules(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 1, jitter= 0.1, shape=1 ,col = "Species",  palette = "Set2",border.lwd = 2, scale = 1.5)

#tmap_arrange(map_1,map_2,map_3)

#Base
map_1
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-1.png)<!-- -->

``` r
#By size too small
map_2
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-2.png)<!-- -->

``` r
#By shape with quantity as color
map_3
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-3.png)<!-- -->

``` r
#Scaled up by size and miniscule jitter
map_4
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-4.png)<!-- -->

``` r
#Circles scaled and jitter
map_5
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-5.png)<!-- -->

``` r
#Base redone with scaling and jitter
map_6
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-6.png)<!-- -->

``` r
#Base with empty circles
map_7
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-7.png)<!-- -->

``` r
#Filled in small circles
map_1 <- tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey95") +
  tm_grid(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 0.2, jitter=0.05,col = "Species",  palette = "Set2") 

tm_shape(world) + tm_polygons()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/finaladjust-8.png)<!-- -->

### Add continental Shelf Data

``` r
#Adding continental shelf and georges bank contours if avaliable
#Did not work
url <- "https://services5.arcgis.com/YCKZmP6DZSq3o9Gi/arcgis/rest/services/Continental_Shelf/FeatureServer/0"
shelf <- read.csv(url)
#read in shapefile

#USGS FILE FROM:
#https://www.sciencebase.gov/catalog/item/4fb69a6be4b03ad19d64b4c2

#Remember to include all the other similar files
bary_iso <- st_read("/Users/hannah/gitHub/Spisula/Work/shape/bathymetry_l_v2.shp")
```

    ## Reading layer `bathymetry_l_v2' from data source `/Users/hannah/gitHub/Spisula/Work/shape/bathymetry_l_v2.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 4197 features and 5 fields
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: -5835002 ymin: -3920000 xmax: 4462000 ymax: 4907000
    ## proj4string:    +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs

``` r
#basic tmap with shapefile
tm_shape(bary_iso) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-1.png)<!-- -->

``` r
#Works, now try to crop it to my map shape
tm_shape(bary_iso, bbox=eastcoast) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-2.png)<!-- -->

``` r
#came up blank, probably different coordinate systems

#Try with boston example
#From: https://www.kaggle.com/umeshnarayanappa/quickstart-guide-shapefiles-and-r-tmap

boston <- st_read(dsn ="/Users/hannah/gitHub/Spisula/Work/shape/boston_police_districts_f55.shp")
```

    ## Reading layer `boston_police_districts_f55' from data source `/Users/hannah/gitHub/Spisula/Work/shape/boston_police_districts_f55.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 12 features and 5 fields
    ## geometry type:  MULTIPOLYGON
    ## dimension:      XY
    ## bbox:           xmin: 739826.9 ymin: 2908285 xmax: 804052.5 ymax: 2970073
    ## CRS:            2249

``` r
#From: http://worldmap.harvard.edu/data/geonode:boston_police_districts_f55
tm_shape(us_states)+tm_polygons()+
  tm_shape(boston) + tm_polygons() + tm_graticules()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-3.png)<!-- -->

``` r
#Tiny little boston

tm_shape(us_states, bbox=eastcoast) +  
  tm_fill(col = "grey99") + tm_borders(col="grey30") + tm_layout(bg.color = "grey95") +
  tm_grid(col="white") +
  tm_shape(spissites) +  
  tm_symbols(size = 0.2, jitter=0.05,col = "Species",  palette = "Set2") +
  tm_shape(boston) + tm_polygons() + tm_graticules()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-4.png)<!-- -->

``` r
#Added boston to my map

#Try adding shelf to world
tm_shape(world) + tm_polygons() +
  tm_shape(bary_iso) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-5.png)<!-- -->

``` r
#Some issue with where they each draw the lines for the pacific

#Try with us_States
tm_shape(us_states) + tm_polygons() +
  tm_shape(bary_iso) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-6.png)<!-- -->

``` r
#Close except the pacific line is still gross.
#I could take wild guesses to remove those lines
#tm_shape(bary_iso[1,]) + tm_lines()
bary_map2 <- tm_shape(us_states) + tm_polygons() + tm_shape(bary_iso[c(1:3170,3172:4000),]) + tm_lines()
bary_map2
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-7.png)<!-- -->

``` r
#THE BAD LINE IS IN ROW 3171

#Change the shape file
bary_us <- bary_iso[-3171,]
tm_shape(us_states) + tm_polygons() + tm_shape(bary_us) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-8.png)<!-- -->

``` r
#Now try with east coast
tm_shape(us_states, bbox=eastcoast) +  tm_polygons() + tm_shape(bary_us) + tm_lines()
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/continental%20shelf-9.png)<!-- -->

``` r
#Woohooo!!
```

### Final Options

``` r
tm_shape(us_states, bbox=eastcoast) + tm_fill(alpha = 0) +
  tm_graticules(col = "white", lwd = 0.7) +
  tm_shape(bary_us, bbox=eastcoast) + tm_lines(col = "grey30", lwd = 1) +
  tm_shape(us_states, bbox=eastcoast) +
  tm_fill(col = "grey86") + tm_borders(col="grey30") + tm_layout(bg.color = "grey99") +
  tm_shape(spissites) +  
  tm_symbols(size = 0.2, jitter=0.05,col = "Species",  palette = "Set2", border.col = "grey30")
```

![](Spisula_Map_QReport_2020_Dec_files/figure-gfm/applycomments-1.png)<!-- -->

#Load the plotting package
library(ggplot2)
#Load the "simple feature" package that makes playing with GIS data easy
library(sf)

#Set directory
setwd("~/Documents/ACSF oyster restoration/HRF_2019/R.Map")
#Read in data
spats = read.csv("Spat_data.csv", header = T)
head(spats)

#Load shapefiles
#First the polygon file with of NY estuaries
Water = st_read(dsn=getwd(), layer = "NY_EstuaryOcean", as_tibble = F, quiet = T)
#Then the line file of NYC bridges
Bridges = st_read(dsn=getwd(), layer = "NYC_Bridges", as_tibble = F, quiet = T)


#Then build a plot that's going to be saved as "p1" (for ease of reference when it comes to exporting it as a png)
p1 = ggplot()+
  #Lay down the estuary first, color it grey (both fill and outline)
  geom_sf(data = Water, fill = 'grey80', col = 'grey80' )+
  #Then the bridges. Let's make them a darker grey. No need to worry about fill since it's a line
  geom_sf(data = Bridges, size =0.75, col = 'grey50' )+
  #Next our points, set point size to be a function of spat count
  geom_point(data = spats, aes(x = Lon, y = Lat, size = Count), fill = 'black', col = 'black')+
  #Now we'll subset the spat points which are zeroes and put a white-filled point on top of the ones we just did, to make them easy to see
  geom_point(data = spats[which(spats$Count == 0),], aes(x = Lon, y = Lat), pch = 21, fill = 'white', size = 2)+
  
  #Within points label stuff is hard to make work right, commented out for now
  #geom_text(data = spats[which(spats$Count == 0),], aes(x = Lon, y = Lat, label = paste(Sal, sep = '')), size = 1.5, col = 'black')+
  #geom_text(data = spats[which(spats$Count >0),], aes(x = Lon, y = Lat, label = paste(Count,'\n',Sal, sep = '')), size = 1.5, col = 'white')+
  
  #Set the breaks for our proportional sizes, and the range of possible sizes. Suppress label (for legend). c(2,6) controls size range for circles
  scale_size_continuous("",breaks = c(0, 10, 50, 100, 250, 500, 1000),range = c(2, 7))+
  #Make the X breaks a little easier to understand. Suppress label
  scale_x_continuous("", breaks = seq(-74.2,-73.8, by = 0.1))+
  #Suppress Y label
  scale_y_continuous("")+
  #Make the map show only the region around the spat points, plus a little bit on the top and right (otherwise would default to extent of largest layer)
  #coord_sf(xlim = c(min(spats$Lon), max(spats$Lon)+ 0.02), ylim = c(min(spats$Lat), max(spats$Lat)+ 0.02))+
  coord_sf(xlim = c(-74.2, -73.75), ylim = c(40.50, 41.08))+
  #Make our legend show the zero as a white point
  guides(size = guide_legend(override.aes = list(shape = 21, fill = c('white', 1,1,1,1,1))))+
  #Format the legend (put it in the top left, stick a box around it)
  theme(legend.position = c(0.15,0.75), legend.title = element_text(colour="black", size=3),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

#Show the saved plot
p1

#Launch a new graphics object that's a png file
png("Map_draft_1.png", height = 8, width = 5, units = 'in', res = 600)
#Write the plot to it
p1
#Finish the file and close it
dev.off()

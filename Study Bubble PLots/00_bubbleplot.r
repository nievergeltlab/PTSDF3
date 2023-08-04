
# Libraries
library(packcircles)
library(ggplot2)
 library(viridis)
 library(data.table)
 
 
 blue=rgb(153,191,254,max=255)
 red=rgb(230,185,184,max=255)
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)
 
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)

adam_blue <-  rgb(60,160,255,maxColorValue=255)
#beauty_red <- rgb(244,122,124,maxColorValue=255)


 d1 <- fread('f3_analyzed_sets_V2.csv',data.table=F)
 d1$color <- "grey"
 d1[which(d1$Race=="EUR"),]$color <- adam_blue
 d1[which(d1$Race=="AAM"),]$color <- beauty_red
 d1[which(d1$Race=="HNA"),]$color <- green
 
 d1$value <- d1$N
 
 d1$id <- d1$Abbr

 PGC <- subset(d1,Group=="PGC2")
 MVP <- subset(d1, Group=="MVP")
 EHR <- subset(d1,Group=="EHR")

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value

#PLOT MVP 

data <- MVP 
packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
#jpeg("mvp.jpeg",7,7)
jpeg("mvp.jpg",1000,1000,quality=100)

ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  

#Plot PGC


#data <- PGC[order(PGC$N,decreasing=TRUE),]
        
#put other cohorts in differnet order.. UKBB at end, others random

data <- data[c(sample(1:(nrow(PGC)-1),replace=F),nrow(PGC)),]

data <- PGC[order(PGC$N,decreasing=FALSE),]
   
data <- data[c(1,sample(2:(nrow(PGC)),replace=F)),]


packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("PGC.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  
#plot EHR


data <- EHR[order(EHR$N,decreasing=FALSE),]
packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("EHR.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  
#Plot the metas

data1 <- subset(d1,Race=="EUR")
data <- data1[order(data1$N,decreasing=TRUE),]
packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=150)
#make plots for each race
# Make the plot
jpeg("eur_meta.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  
  
data1 <- subset(d1,Race=="AAM")
data <- data1[order(data1$N,decreasing=TRUE),]

packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("aam_meta.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
 #lat meta
 
 
data1 <- subset(d1,Race=="HNA")
data <- data1[order(data1$N,decreasing=TRUE),]

packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("hna_meta.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  #trans meta
  
  data1 <- subset(d1,Race=="HNA")
data <- data1[order(data1$N,decreasing=TRUE),]

packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("hna_meta.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  
  #everyone
    data1 <- subset(d1)
data <- data1[order(data1$N,decreasing=TRUE),]

packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)
# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)
 
# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)
#make plots for each race
# Make the plot
jpeg("all_meta.jpeg",1000,1000,quality=100)
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 1) +
  scale_fill_manual(values = data$color) +
 
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = Abbr)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
  
  dev.off()
  
  
  
  
 
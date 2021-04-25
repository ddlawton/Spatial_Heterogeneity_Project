#####
# Formatting 
#   for Fractal
#    analysis
###
rm(list=ls())
library(tidyverse)
library(data.table)
library(lubridate)
library(raster)
library(sp)
library(sf)
library(svMisc) #oh this is a cool package I just discovered that gives a progress bar to be used in for loops

#Fractal function (from Cyril Piou)

fractalD <-
  function(locRASTER,threshold=NULL,output=FALSE) 
  {
    smlDim=min(c(nrow(locRASTER),ncol(locRASTER)))
    m=matrix(locRASTER,nrow=nrow(locRASTER),ncol=ncol(locRASTER))
    if(nrow(locRASTER)!=ncol(locRASTER)){
      if(nrow(locRASTER)>ncol(locRASTER)){m<-m[1:smlDim,]}else{m<-m[,1:smlDim]}
    }
    powerl<- 1
    while(2^powerl<smlDim)powerl<-powerl+1
    powers <- 1:(powerl-1)
    boxsize <- 2^powers
    imgsize <- max(boxsize) 
    if(is.null(threshold)){threshold<-mean(m)}
    if(!any(m>=threshold)){return(0)}else{
      tsb<-m
      tsb[m >= threshold]<-1
      tsb[m <  threshold]<-0
      bv <- c()
      if(output)print('Box Size   Count   Total Box')
      for (b in boxsize) {
        bcount <- 0
        bstart <- seq(1,imgsize, by=b)
        totalbox <- 0
        for (x in bstart) {
          for (y in bstart) {
            totalbox <- totalbox+1
            xend <- x + b - 1
            yend <- y + b - 1
            if (any( tsb[ x:xend, y:yend ] > 0)) {
              bcount <- bcount + 1
            }
          }
        }
        if(output)print(paste(b,bcount,totalbox))
        bv <- c(bv, bcount)
      }
      if(output)plot(log(1/boxsize), log(bv))
      dfl <- data.frame(x=log(1/boxsize), y=log(bv))
      dfl<-dfl[is.finite(dfl$y),]
      if(nrow(dfl)<=2){return(0)}
      model <- lm(y ~  x, dfl)
      if(output){abline(model)
        print(paste('fractal dimension =',model$coeff[2]))}
      #list(fdim=model$coeff[2],model=model)
      return(ifelse((is.na(coef(summary(model))[2,4]) | coef(summary(model))[2,4] > 0.05),0,model$coeff[2]))
    }
  }


## Getting prior raster dates lined up with observation dates and going back 10 images

CT <- read.csv("data/raw/CT_unixdates.csv") #This is the APL dataset
dates <- read.csv("data/raw/MODIS Image/Dates.csv") #This is the dates and band # for the associated raster 

CT$Date <- as.Date(CT$Date,format="%m/%d/%y")


dates$date <- as.Date(dates$date,format="%Y-%m-%dT%H:%M:%S")

#convert sample data to data.table
# This is just to find the nearest band date to observation

setDT(CT) 
setDT(dates)

#convert dates to 'real' dates
CT[, date := as.IDate(Date) ]
dates[, date := as.IDate(date) ]

#update CT by reference with a rolling join
CT[, BandNumber_1 := dates[ CT, x.BandNumber, on = .(date), roll = Inf ] ]

#Now, I want to go back 9 additional images from the closest image to observation date

CT <- CT %>% mutate(
  BandNumber_2 = BandNumber_1 - 1,
  BandNumber_3 = BandNumber_1 - 2,
  BandNumber_4 = BandNumber_1 - 3,
  BandNumber_5 = BandNumber_1 - 4,
  BandNumber_6 = BandNumber_1 - 5,
  BandNumber_7 = BandNumber_1 - 6,
  BandNumber_8 = BandNumber_1 - 7,
  BandNumber_9 = BandNumber_1 - 8,
  BandNumber_10 = BandNumber_1 - 9
) %>% tidyr::pivot_longer((starts_with("BandNumber_")), names_to = "Band_time_step",values_to = "BandNumber")


CT <- left_join(CT,dates,by="BandNumber") #This gives me the image collection dates


CT <- CT %>% rename(   
  Band_date = date.y
) %>% drop_na(BandNumber) #this drops any supposed band numbers from the mutate() command that do not actually show up in raster

CT$diff_days <- as.numeric(difftime(CT$Band_date, CT$Date, units = "days")) #this calculates the difference between locust observation date and image collection dates


### Now starting to compute Fractal stats

spdf <- SpatialPointsDataFrame(coords=CT[,4:5],data=CT,
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

spdf2 <- spTransform(spdf,crs("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))




MODIS_raster_location <- "data/raw/MODIS Image/MODIS_NDVI.tif"

spdf3 <- st_as_sf(spdf2)
spdf3_buffed <- st_buffer(spdf3, dist = 2500, endCapStyle = 'SQUARE') #this creates a 5km x 5km square buffer around the point. For CT this might be too big?

# we are calculating Morans I, Fractal D, and Means in one big for-loop. Is there a better way to do this? probably....
# Moran's I is coming from the raster package
#Fractal D is from your provided function.

spdf3_buffed$FractID <- 0
spdf3_buffed$meanNDVI$Moran <- 0
spdf3_buffed$meanNDVI <- 0

for (i in 1:nrow(spdf3_buffed)){
  progress(i) #calculates progress
  spSubset <- spdf3_buffed[i,]
  BandNumber <-spSubset$BandNumber
  r <- raster(MODIS_raster_location,band=BandNumber) 
  cropped_raster <- crop(r, st_bbox(spSubset))
  b <- as.numeric(fractalD(cropped_raster))
  c <- as.numeric(Moran(cropped_raster))
  layermeans <- cellStats(cropped_raster, stat='mean', na.rm=TRUE)
  u <- mean(layermeans)/100000 #divided by 100,000 to get into normal NDVI terms.
  spdf3_buffed$FractID[[i]] <- b
  spdf3_buffed$Moran[[i]] <- c
  spdf3_buffed$meanNDVI[[i]] <- u
}

write.csv(spdf3_buffed,file="test.csv")


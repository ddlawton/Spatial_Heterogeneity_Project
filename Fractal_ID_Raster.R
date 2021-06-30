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
library(rgeos)
library(rgdal)
library(svMisc) #oh this is a cool package I just discovered that gives a progress bar to be used in for loops

#Fractal function (from Cyril Piou)

fractalD <-
  function(locRASTER,threshold=NULL,output=FALSE){
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

names(CT)
target <- c("11","10")
CT2 <- CT %>% filter(Species %in% target)

### Now starting to compute Fractal statsa

spdf <- SpatialPointsDataFrame(coords=CT2[,4:5],data=CT2,
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

spdf2 <- spTransform(spdf,crs("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

plot(r)
plot(CT2, add=TRUE)

MODIS_raster_location <- "data/raw/MODIS Image/MODIS_NDVI.tif"

r <- raster(MODIS_raster_location, band=1)

spdf3 <- st_as_sf(spdf2)
spdf3_buffed <- st_buffer(spdf3, dist = 2500, endCapStyle = 'SQUARE') #this creates a 5km x 5km square buffer around the point. For CT this might be too big?




library(exactextractr)

spSubset <- spdf3_buffed %>% filter(BandNumber == 788)
BandNumber <-spSubset$BandNumber
r_b1 <- raster(MODIS_raster_location,band=BandNumber) 

exact_extract(
  x=r_b1,
  y=spSubset,
  progress=TRUE,
  summarize_df = TRUE,
  fun = function(x) Moran(x)
)




spdf3_buffed <- spdf3_buffed %>% rowid_to_column()


library(spex)



projection(spdf3_buffed) == projection(p2)
spsel <- st_intersection(spdf3_buffed, p2)

plot(r)
plot(spsel,add=TRUE)

# we are calculating Morans I, Fractal D, and Means in one big for-loop. Is there a better way to do this? probably....
# Moran's I is coming from the raster package
#Fractal D is from your provided function.

spdf3_buffed$FractID <- 0
spdf3_buffed$meanNDVI$Moran <- 0
spdf3_buffed$meanNDVI <- 0

for (i in 1:nrow(spdf3_buffed)){
  print(i) #calculates progress
  spSubset <- spdf3_buffed[i,]
  BandNumber <-spSubset$BandNumber
  r <- raster(MODIS_raster_location,band=BandNumber) 
  cropped_raster <- crop(r, (spSubset))
  b <- as.numeric(fractalD(cropped_raster))
  c <- as.numeric(Moran(cropped_raster))
  layermeans <- cellStats(cropped_raster, stat='mean', na.rm=TRUE)
  u <- mean(layermeans)/100000 #divided by 100,000 to get into normal NDVI terms.
  spdf3_buffed$FractID[[i]] <- b
  spdf3_buffed$Moran[[i]] <- c
  spdf3_buffed$meanNDVI[[i]] <- u
}



# Putting it into parallel
library(doParallel)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

registerDoParallel(cl <- makeCluster(7))
sampled <- spdf3_buffed %>% sample_n(size=100)
start <- Sys.time()
fract_list <-list()
Moran_list <-list()
NDVI_mean_list <-list()
data <- foreach(i = 1:nrow(sampled),.packages=c('raster','sf'),.combine="comb", .multicombine=TRUE,.init=list(list(), list(), list())) %dopar% {
  spSubset <- sampled[i,]
  BandNumber <-spSubset$BandNumber
  r <- raster(MODIS_raster_location,band=BandNumber) 
  cropped_raster <- crop(r, st_bbox(spSubset))
  fract_list <- fractalD(cropped_raster)
  Moran_list <- Moran(cropped_raster)
  NDVI_mean_list <- mean(cellStats(cropped_raster, stat='mean', na.rm=TRUE))
  list(fract_list,Moran_list,NDVI_mean_list)
}

stopCluster(cl)
end <- Sys.time()
end-start


oper1 <- data[[1]]
oper2 <- data[[2]]


## playing around
library(doSNOW)


sampled <- spdf3_buffed %>% sample_n(size=100) %>% dplyr::select(!c("FractID","Moran","meanNDVI"))
r_stacked <- stack("data/raw/MODIS Image/MODIS_NDVI.tif")
library(doSNOW)
cl <- makeCluster(7)
registerDoSNOW(cl)
iterations <- nrow(sampled)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result <- foreach(i = 1:nrow(sampled),.packages=c('raster','sf'),.combine="comb", .multicombine=TRUE,.init=list(list(), list(), list()),
                  .options.snow = opts) %dopar%{
    spSubset <- sampled[1,]
    BandNumber <-spSubset$BandNumber
    r <- r_stacked[[i]]
    cropped_raster <- crop(r, st_bbox(spSubset))
    cropped_raster2 <- mask(cropped_raster, (spSubset))
    cropped_raster2[cropped_raster2 < 10000] <- NA
    fract_list <- fractalD(cropped_raster)
    Moran_list <- Moran(cropped_raster)
    NDVI_mean_list <- mean(cellStats(cropped_raster, stat='mean', na.rm=TRUE))
    list(fract_list,Moran_list,NDVI_mean_list)
}
close(pb)
stopCluster(cl) 

save(result, file="results.RData")

result_df<-as.data.frame(t((do.call(rbind.data.frame, t(result)))))

sampled$Fract <- result_df$`2`
sampled$Moran <- result_df$`21`
sampled$Mean_NDVI <- result_df$`3`

write.csv(sampled,file="Sample2.csv")
st_write(spdf3_buffed, paste0(getwd(), "/", "Test_apr26.csv"))

str(spdf3_buffed)
fwrite(spdf3_buffed)


str(spdf3_buffed)
# Eh
library(future.apply)

sampled <- spdf3_buffed %>% sample_n(size=100)

r_stacked <- stack("data/raw/MODIS Image/MODIS_NDVI.tif")
r_stacked[[1]]



s <- lapply(1:nrow(sampled), function(i) {
  cat("processing", i, "of", nrow(sampled), "\n")
  Sys.sleep(0.01)
  flush.console()
  rsub <- crop(r_stacked[[sampled[i,]$BandNumber]], extent(sampled[i,]))
  fd <- as.data.frame(landscapemetrics::lsm_c_frac_mn(rsub))[,6][2]
  c(fd, as.numeric(Moran(rsub)), mean(rsub[], na.rm=TRUE)) } )

s <- future_lapply(1:nrow(sampled), function(i) {
  rsub <- crop(r_stacked[[sampled[i,]$BandNumber]], extent(sampled[i,]))
  fd <- as.data.frame(landscapemetrics::lsm_c_frac_mn(rsub))[,6][2]
  c(fd, as.numeric(Moran(rsub)), mean(rsub[], na.rm=TRUE)) } )





#This document shows our pre-processing of the New York Citi Bike data.
#Data is available at https://ride.citibikenyc.com/system-data.
#Running this R file requires Google map API. You can download NYbike.rda directly to obtain the pre-processed data.
NYbike <- read.table("/201908-citibike-tripdata.csv",header = TRUE,sep = ',')
Aug1 <- NYbike[as.Date(NYbike$starttime) == "2019-08-01",]
Aug1.del <- Aug1[Aug1$start.station.id != Aug1$end.station.id & Aug1$tripduration > 60 & Aug1$tripduration < 3*60*60,c(2,4,8)]
Aug1.del$starttime <- (as.numeric(strptime(Aug1.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-01 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

#R package ggmap (need Google map API)
library(ggmap)
library(httr)
library(ggplot2)
library(ggpubr)
register_google(key = "********", write = TRUE)

#782 stations and 85854 rides
stations <- unique(sort(as.numeric(c(Aug1.del$start.station.id,Aug1.del$end.station.id))))
location <- NULL
for(i in stations)
  location <- rbind(location,Aug1[Aug1$start.station.id == i,c(7,6)][1,])
names(location) <- c("lon", "lat")

borough <- rep(0,782)
for(i in 1:782){
  add <- revgeocode(c(location[i,1],location[i,2]),output='all')$results[[1]]$address_components
  for(j in 1:length(add)){
    if(add[[j]]$short_name == "Manhattan")
      borough[i] <- "Manhattan"
    if(add[[j]]$short_name == "Brooklyn")
      borough[i] <- "Brooklyn"
    if(add[[j]]$short_name == "Queens")
      borough[i] <- "Queens"
  }
}
borough[which(borough == 0)] <- "Manhattan"
location$borough <- factor(borough)

#plot real locations
set_config(
  use_proxy(url="127.0.0.1", port=1080)
)
mp <- ggmap(get_googlemap(center=c(lon=-73.95,lat=40.734),zoom=12,size=c(640,640),maptype='terrain'))
mp <- mp + geom_point(aes(x=lon,y=lat,col=borough,shape=borough),data = location,size=0.8) 
mp <- mp + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD'))
mp <- mp + scale_shape_manual(values = c(18,16,15))
mp <- mp + theme(legend.position = c(0.86,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size = unit(12, "pt"))
mp <- mp + labs(x='longitude',y='latitude')

#get network data
NYbike <- NULL
for(i in stations)
  for(j in stations){
    if(i == j)
      NYbike <- c(NYbike,list(numeric(0)))
    else
      NYbike <- c(NYbike,list(Aug1.del$starttime[Aug1.del$start.station.id==i & Aug1.del$end.station.id==j]))
  }

#discrete time
NYbike1 <- NULL
for(t in 0:23)
  NYbike1 <- c(NYbike1,sapply(NYbike,function(x) sum(I(x > t & x <= (t+1)))))
NYbike1 <- array(NYbike1,dim = c(782,782,24))

#save the pre-processed data in NYbike.rda
save(NYbike1,borough,file = "NYbike.rda")

NYbike <- read.table("C:/Users/ÌïÓî°º/Desktop/SBM/real data/NYbike/201908-citibike-tripdata.csv",header = TRUE,sep = ',')
Aug1 <- NYbike[as.Date(NYbike$starttime) == "2019-08-01",]
Aug1.del <- Aug1[Aug1$start.station.id != Aug1$end.station.id & Aug1$tripduration > 60 & Aug1$tripduration < 3*60*60,c(2,4,8)]
Aug1.del$starttime <- (as.numeric(strptime(Aug1.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-01 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

#R package ggmap (need google map API)
library(ggmap)
library(httr)
library(ggplot2)
library(ggpubr)
register_google(key = "AIzaSyA-sBZmhdTejXLNUuU_i2OhqJ2GLXdXm3U", write = TRUE)

#782 stations 85854 rides
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

#estimate and rotate
NYbike.est1 <- PGD.panel(NYbike1,1)
NYbike.est2 <- PGD.panel(NYbike1,2)
NYbike.est3 <- PGD.panel(NYbike1,3)
NYbike.est4 <- PGD.panel(NYbike1,4)

#choose k=3?
logLiki <- c(NYbike.est1$logL,NYbike.est2$logL,NYbike.est3$logL,NYbike.est4$logL)
k <- 1:4
plot(k,-2*logLiki + (k*781 - k*(k-1)/2)*log(782*781*24),xlab = 'k',ylab = 'BIC')
#plot(NYbike.est3$Z[,1],NYbike.est3$Z[,2])

theta <- (1-1/24) * pi
Z.rotate <- NYbike.est2$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate <- as.data.frame(Z.rotate)
names(Z.rotate) <- c("Z1", "Z2")
Z.rotate$borough <- factor(borough)
mp2 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough),data = Z.rotate,size=0.8) + xlim(-3.5,3.5)
mp2 <- mp2 + scale_color_manual(values = c('#FF0000','#6495ED','#9A32CD')) + theme(legend.position="none")

#method without diag
NYbike.est21 <- PGD.panel2(NYbike1,1)
NYbike.est22 <- PGD.panel2(NYbike1,2)
NYbike.est23 <- PGD.panel2(NYbike1,3)
NYbike.est24 <- PGD.panel2(NYbike1,4)

logLiki2 <- c(NYbike.est21$logL,NYbike.est22$logL,NYbike.est23$logL,NYbike.est24$logL)
k <- 1:4
plot(k,-2*logLiki2 + (k*781 - k*(k-1)/2)*log(782*781*24),xlab = 'k',ylab = 'BIC')

theta <- (1 - 1/24) * pi
Z.rotate2 <- NYbike.est22$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate2 <- as.data.frame(Z.rotate2)
names(Z.rotate2) <- c("Z1", "Z2")
Z.rotate2$borough <- factor(borough)
mp3 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate2,size = 0.8) + xlim(-3.1,4.9)
mp3 <- mp3 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp3 <- mp3 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp3 <- mp3 + theme(legend.position = c(0.8,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))
#NY1 <- ggarrange(mp, mp3, ncol = 2, nrow = 1, widths = c(1.5,1))
#ggsave("NYbike1.eps", NY1, device = cairo_ps)

#other est
NYbike.ests <- PGD.panel2(apply(NYbike1,c(1,2),mean),2)
Z.rotate3 <- NYbike.ests$Z %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate3 <- as.data.frame(Z.rotate3)
names(Z.rotate3) <- c("Z1", "Z2")
Z.rotate3$borough <- factor(borough)
mp4 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate3,size = 0.8) + xlim(-3.5,4.5)
mp4 <- mp4 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp4 <- mp4 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp4 <- mp4 + theme(legend.position = c(0.81,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

NYbike.time9 <- PGD.panel2(NYbike1[,,9],2)
theta <- -1/12 * pi
Z.rotate5 <- NYbike.time9$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate5 <- as.data.frame(Z.rotate5)
names(Z.rotate5) <- c("Z1", "Z2")
Z.rotate5$borough <- factor(borough)
mp5 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate5,size = 0.8) + xlim(-3.5,4.5)
mp5 <- mp5 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp5 <- mp5 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp5 <- mp5 + theme(legend.position = c(0.81,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

NYbike.time19 <- PGD.panel2(NYbike1[,,19],2)
theta <- (1-1/24)*pi
Z.rotate6 <- NYbike.time19$Z %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate6 <- as.data.frame(Z.rotate6)
names(Z.rotate6) <- c("Z1", "Z2")
Z.rotate6$borough <- factor(borough)
mp6 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate6,size = 0.8) + xlim(-3.5,4.5)
mp6 <- mp6 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp6 <- mp6 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp6 <- mp6 + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

NY2 <- ggarrange(mp5, mp6, ncol = 2, nrow = 1)
ggsave("NYbike2.eps", NY2, device = cairo_ps)

NYbike.time22 <- PGD.panel2(NYbike1[,,22],2)
theta <- (1-1/24)*pi
Z.rotate7 <- NYbike.time22$Z  %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate7 <- as.data.frame(Z.rotate7)
names(Z.rotate7) <- c("Z1", "Z2")
Z.rotate7$borough <- factor(borough)
mp7 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate7,size = 0.8) + xlim(-3.5,4.5) + ylim(-4.15,4.15)
mp7 <- mp7 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp7 <- mp7 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp7 <- mp7 + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

NYbike.time24 <- PGD.panel2(NYbike1[,,24],2)
theta <- (1/12-1) * pi
Z.rotate8 <- NYbike.time24$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate8 <- as.data.frame(Z.rotate8)
names(Z.rotate8) <- c("Z1", "Z2")
Z.rotate8$borough <- factor(borough)
mp8 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate8,size = 0.8) + xlim(-3.5,4.5) 
mp8 <- mp8 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp8 <- mp8 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp8 <- mp8 + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

NY21 <- ggarrange(mp7, mp8, ncol = 2, nrow = 1)
ggsave("NYbike21.eps", NY21, device = cairo_ps)

#plot alpha
NY_t <- apply(NYbike1,3,sum)
alpha_t <- apply(NYbike.est22$alpha,2,mean)
alpha_sd <- apply(NYbike.est22$alpha,2,sd)
ealpha_t <- apply(exp(NYbike.est22$alpha),2,mean)
ealpha_sd <- apply(exp(NYbike.est22$alpha),2,sd)
expalpha_t <- apply(NYbike.est22$alpha,2,function(x) mean(as.vector(exp(outer(x,x,'+')))))
expalpha_sd <- apply(NYbike.est22$alpha,2,function(x) sd(as.vector(exp(outer(x,x,'+')))))

NY_overt <- data.frame(data = NY_t,est = alpha_t,sd = alpha_sd,e_est = ealpha_t,e_sd = ealpha_sd,exp_est = expalpha_t,exp_sd = expalpha_sd,T = 0.5:23.5)
data_t <- ggplot(aes(x=T,y=data),data = NY_overt) + geom_point() + geom_line()  + labs(x='Time',y="Total number of rides")
data_t <- data_t + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))
alpha_est <- ggplot(aes(x=T,y=est),data = NY_overt) + geom_point() + geom_line() + labs(x='Time',y=expression(paste("Mean of ",hat(alpha[i]))))
alpha_est <- alpha_est #+ geom_errorbar(aes(ymin = est - sd,ymax = est + sd))
alpha_est <- alpha_est + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))
ealpha_est <- ggplot(aes(x=T,y=e_est),data = NY_overt) + geom_point() + geom_line() + labs(x='Time',y=expression(paste("Mean of ",exp(hat(alpha[i])))))
ealpha_est <- ealpha_est #+ geom_errorbar(aes(ymin = e_est - e_sd,ymax = e_est + e_sd))
ealpha_est <- ealpha_est + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))
expalpha_est <- ggplot(aes(x=T,y=exp_est),data = NY_overt) + geom_point() + geom_line() + labs(x='Time',y=expression(paste("Mean of ",exp(hat(alpha[i]) + hat(alpha[j])))))
expalpha_est <- expalpha_est #+ geom_errorbar(aes(ymin = exp_est - exp_sd,ymax = exp_est + exp_sd))
expalpha_est <- expalpha_est + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))

p <- ggplot(NY_overt, aes(x=T)) +
  geom_line(aes(y=NY_t/NY_t[1], color="Ride Counts")) +
  geom_line(aes(y=expalpha_t/expalpha_t[1], color="Estimated Baseline")) +
  labs(x="Time", y="Value") + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00')) + 
  scale_color_manual(values=c("red","blue")) # ÉèÖÃÕÛÏßÑÕÉ«
p <- p + guides(color=guide_legend(title=NULL)) + theme(legend.position = c(0.21,0.91)) 

pdf(file='Bike_alpha2.pdf', width=4, height=4)
p
dev.off()

#hist for all stations
NY_n <- data.frame(data = apply(NYbike1,1,sum))
data_n <- ggplot() + geom_histogram(data = NY_n,aes(x=data),stat = 'bin',bins = 100,fill='darkgreen',color='gray')
data_n <- data_n + labs(x = 'Total number of rides',y = 'Count of nodes')

#est k
est_k1 <- USVT(NYbike1)
est_k1 <- data.frame(value = eigen(G_hat1)$values[1:10],num = 1:10)
USVT_estk <- ggplot(aes(x=num,y=value),data = est_k1) + geom_point() + geom_line() + labs(x="Component Number",y="Eigenvalue")
USVT_estk <- USVT_estk + scale_x_continuous(breaks = c(2,4,6,8,10))
G_hat2 <- PGD.G(NYbike1,sqrt(782*24))$G
est_k2 <- data.frame(value = eigen(G_hat2)$values[1:10],num = 1:10)
PMLE_estk <- ggplot(aes(x=num,y=value),data = est_k2) + geom_point() + geom_line() + labs(x="Component Number",y="Eigenvalue")
PMLE_estk <- PMLE_estk + scale_x_continuous(breaks = c(2,4,6,8,10))
NY4 <- ggarrange(USVT_estk, PMLE_estk, ncol = 2, nrow = 1)
ggsave("NYbike4.eps", NY4, device = cairo_ps)

#estimate using PMLE
NYbike_PMLE <- PGD.G(NYbike1,0.1*sqrt(782*24))
NYbike_PMLE1 <- PGD.G(NYbike1,0.01*sqrt(782*24))
G.eigen <- eigen(NYbike_PMLE1$G)
Z_PMLE <- G.eigen$vectors  %*% rbind(diag(sqrt(G.eigen$values[1:2])),matrix(0,780,2))
theta <-  -1/4 * pi
Z.PMLE <- Z_PMLE %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.PMLE <- as.data.frame(Z.PMLE)
names(Z.PMLE) <- c("Z1", "Z2")
Z.PMLE$borough <- factor(borough)
mp_pmle <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.PMLE,size = 0.8) + xlim(-1.66,2.5)
mp_pmle <- mp_pmle + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mp_pmle <- mp_pmle + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mp_pmle <- mp_pmle + theme(legend.position = c(0.81,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

pdf(file='Bike2_Z.pdf', width=3, height=4)
mp_pmle
dev.off()

expalpha_t <- apply(NYbike_PMLE1$alpha,2,function(x) mean(as.vector(exp(outer(x,x,'+')))))
NY_overt <- data.frame(exp_est = expalpha_t,T = 0.5:23.5)
expalpha_est <- ggplot(aes(x=T,y=exp_est),data = NY_overt) + geom_point() + geom_line() + labs(x='Time',y=expression(paste("Mean of ",exp(hat(alpha[i]) + hat(alpha[j])))))
expalpha_est <- expalpha_est + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))

pdf(file='Bike2_alpha.pdf', width=4, height=4)
expalpha_est
dev.off()

#pdf version
par(mar=c(4.5, 5, 1.4, 0.5))

pdf(file='Bike_map.pdf', width=4, height=4)
mp
dev.off()

pdf(file='Bike_datat.pdf', width=4, height=4)
data_t
dev.off()

pdf(file='Bike_datan.pdf', width=4, height=4)
data_n
dev.off()

pdf(file='Bike_estjoint.pdf', width=3, height=4)
mp3
dev.off()

pdf(file='Bike_merge.pdf', width=3, height=4)
mp4
dev.off()

pdf(file='Bike_est9.pdf', width=3, height=4)
mp5
dev.off()

pdf(file='Bike_est19.pdf', width=3, height=4)
mp6
dev.off()

pdf(file='Bike_est22.pdf', width=3, height=4)
mp7
dev.off()

pdf(file='Bike_est24.pdf', width=3, height=4)
mp8
dev.off()

pdf(file='Bike_alpha.pdf', width=4, height=4)
alpha_est
dev.off()

pdf(file='Bike_expalpha.pdf', width=4, height=4)
ealpha_est
dev.off()

pdf(file='Bike_expalpha1.pdf', width=4, height=4)
expalpha_est
dev.off()

#additive model
NYbike.am <- PGD.am(NYbike1,2)
theta <- (1-1/24)*pi
Z.rotate7 <- NYbike.am$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate7 <- as.data.frame(Z.rotate7)
names(Z.rotate7) <- c("Z1", "Z2")
Z.rotate7$borough <- factor(borough)
mpam <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate7,size = 0.8)# + xlim(-3.5,4.5) + ylim(-4.15,4.15)
mpam <- mpam + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mpam <- mpam + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mpam <- mpam + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

pdf(file='Bike_am.pdf', width=3, height=4)
mpam
dev.off()

NYbike.cosie <- MASE(NYbike1,2)
theta <- -1/12 * pi
Z.rotate7 <- NYbike.cosie$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate7 <- as.data.frame(Z.rotate7)
names(Z.rotate7) <- c("Z1", "Z2")
Z.rotate7$borough <- factor(borough)
mpco <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate7,size = 0.8)# + xlim(-3.5,4.5) + ylim(-4.15,4.15)
mpco <- mpco + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
mpco <- mpco + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
mpco <- mpco + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

pdf(file='Bike_cosie.pdf', width=3, height=4)
mpco
dev.off()

#data for longer time
Aug2 <- NYbike[as.Date(NYbike$starttime) == "2019-08-02",]
Aug2.del <- Aug2[Aug2$start.station.id != Aug2$end.station.id & Aug2$tripduration > 60 & Aug2$tripduration < 3*60*60,c(2,4,8)]
Aug2.del$starttime <- (as.numeric(strptime(Aug2.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-02 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug3 <- NYbike[as.Date(NYbike$starttime) == "2019-08-03",]
Aug3.del <- Aug3[Aug3$start.station.id != Aug3$end.station.id & Aug3$tripduration > 60 & Aug3$tripduration < 3*60*60,c(2,4,8)]
Aug3.del$starttime <- (as.numeric(strptime(Aug3.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-03 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug4 <- NYbike[as.Date(NYbike$starttime) == "2019-08-04",]
Aug4.del <- Aug4[Aug4$start.station.id != Aug4$end.station.id & Aug4$tripduration > 60 & Aug4$tripduration < 3*60*60,c(2,4,8)]
Aug4.del$starttime <- (as.numeric(strptime(Aug4.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-04 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug5 <- NYbike[as.Date(NYbike$starttime) == "2019-08-05",]
Aug5.del <- Aug5[Aug5$start.station.id != Aug5$end.station.id & Aug5$tripduration > 60 & Aug5$tripduration < 3*60*60,c(2,4,8)]
Aug5.del$starttime <- (as.numeric(strptime(Aug5.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-05 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug6 <- NYbike[as.Date(NYbike$starttime) == "2019-08-06",]
Aug6.del <- Aug6[Aug6$start.station.id != Aug6$end.station.id & Aug6$tripduration > 60 & Aug6$tripduration < 3*60*60,c(2,4,8)]
Aug6.del$starttime <- (as.numeric(strptime(Aug6.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-06 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug7 <- NYbike[as.Date(NYbike$starttime) == "2019-08-07",]
Aug7.del <- Aug7[Aug7$start.station.id != Aug7$end.station.id & Aug7$tripduration > 60 & Aug7$tripduration < 3*60*60,c(2,4,8)]
Aug7.del$starttime <- (as.numeric(strptime(Aug7.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-07 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

ride_num1 <- sapply(0:23,function(x) sum(floor(Aug1.del$starttime)==x))
ride_num2 <- sapply(0:23,function(x) sum(floor(Aug2.del$starttime)==x))
ride_num3 <- sapply(0:23,function(x) sum(floor(Aug3.del$starttime)==x))
ride_num4 <- sapply(0:23,function(x) sum(floor(Aug4.del$starttime)==x))
ride_num5 <- sapply(0:23,function(x) sum(floor(Aug5.del$starttime)==x))
ride_num6 <- sapply(0:23,function(x) sum(floor(Aug6.del$starttime)==x))
ride_num7 <- sapply(0:23,function(x) sum(floor(Aug7.del$starttime)==x))

NY_t <- c(ride_num1,ride_num2,ride_num3,ride_num4,ride_num5,ride_num6,ride_num7)
NY_overt <- data.frame(data = NY_t,T = 0.5:167.5)
data_t <- ggplot(aes(x=T,y=data),data = NY_overt) + geom_point() + geom_line()  + labs(x='Time',y="Total number of rides")
data_t <- data_t + scale_x_continuous(breaks = c(12,36,60,84,108,132,156),labels=c('August 1st','August 2nd','August 3rd','August 4th','August 5th','August 6th','August 7th'))

pdf(file='Bike_longt.pdf', width=12, height=4)
data_t
dev.off()

Aug8 <- NYbike[as.Date(NYbike$starttime) == "2019-08-08",]
Aug8.del <- Aug8[Aug8$start.station.id != Aug8$end.station.id & Aug8$tripduration > 60 & Aug8$tripduration < 3*60*60,c(2,4,8)]
Aug8.del$starttime <- (as.numeric(strptime(Aug8.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-08 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug9 <- NYbike[as.Date(NYbike$starttime) == "2019-08-09",]
Aug9.del <- Aug9[Aug9$start.station.id != Aug9$end.station.id & Aug9$tripduration > 60 & Aug9$tripduration < 3*60*60,c(2,4,8)]
Aug9.del$starttime <- (as.numeric(strptime(Aug9.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-09 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug10 <- NYbike[as.Date(NYbike$starttime) == "2019-08-10",]
Aug10.del <- Aug10[Aug10$start.station.id != Aug10$end.station.id & Aug10$tripduration > 60 & Aug10$tripduration < 3*60*60,c(2,4,8)]
Aug10.del$starttime <- (as.numeric(strptime(Aug10.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-10 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug11 <- NYbike[as.Date(NYbike$starttime) == "2019-08-11",]
Aug11.del <- Aug11[Aug11$start.station.id != Aug11$end.station.id & Aug11$tripduration > 60 & Aug11$tripduration < 3*60*60,c(2,4,8)]
Aug11.del$starttime <- (as.numeric(strptime(Aug11.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-11 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug12 <- NYbike[as.Date(NYbike$starttime) == "2019-08-12",]
Aug12.del <- Aug12[Aug12$start.station.id != Aug12$end.station.id & Aug12$tripduration > 60 & Aug12$tripduration < 3*60*60,c(2,4,8)]
Aug12.del$starttime <- (as.numeric(strptime(Aug12.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-12 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug13 <- NYbike[as.Date(NYbike$starttime) == "2019-08-13",]
Aug13.del <- Aug13[Aug13$start.station.id != Aug13$end.station.id & Aug13$tripduration > 60 & Aug13$tripduration < 3*60*60,c(2,4,8)]
Aug13.del$starttime <- (as.numeric(strptime(Aug13.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-13 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug14 <- NYbike[as.Date(NYbike$starttime) == "2019-08-14",]
Aug14.del <- Aug14[Aug14$start.station.id != Aug14$end.station.id & Aug14$tripduration > 60 & Aug14$tripduration < 3*60*60,c(2,4,8)]
Aug14.del$starttime <- (as.numeric(strptime(Aug14.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-14 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

ride_num8 <- sapply(0:23,function(x) sum(floor(Aug8.del$starttime)==x))
ride_num9 <- sapply(0:23,function(x) sum(floor(Aug9.del$starttime)==x))
ride_num10 <- sapply(0:23,function(x) sum(floor(Aug10.del$starttime)==x))
ride_num11 <- sapply(0:23,function(x) sum(floor(Aug11.del$starttime)==x))
ride_num12 <- sapply(0:23,function(x) sum(floor(Aug12.del$starttime)==x))
ride_num13 <- sapply(0:23,function(x) sum(floor(Aug13.del$starttime)==x))
ride_num14 <- sapply(0:23,function(x) sum(floor(Aug14.del$starttime)==x))

NY_t <- c(ride_num8,ride_num9,ride_num10,ride_num11,ride_num12,ride_num13,ride_num14)
NY_overt <- data.frame(data = NY_t,T = 0.5:167.5)
data_t <- ggplot(aes(x=T,y=data),data = NY_overt) + geom_point() + geom_line()  + labs(x='Time',y="Total number of rides")
data_t <- data_t + scale_x_continuous(breaks = c(12,36,60,84,108,132,156),labels=c('August 8th','August 9th','August 10th','August 11th','August 12th','August 13th','August 14th'))

pdf(file='Bike_longt2.pdf', width=12, height=4)
data_t
dev.off()

Aug15 <- NYbike[as.Date(NYbike$starttime) == "2019-08-15",]
Aug15.del <- Aug15[Aug15$start.station.id != Aug15$end.station.id & Aug15$tripduration > 60 & Aug15$tripduration < 3*60*60,c(2,4,8)]
Aug15.del$starttime <- (as.numeric(strptime(Aug15.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-15 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug16 <- NYbike[as.Date(NYbike$starttime) == "2019-08-16",]
Aug16.del <- Aug16[Aug16$start.station.id != Aug16$end.station.id & Aug16$tripduration > 60 & Aug16$tripduration < 3*60*60,c(2,4,8)]
Aug16.del$starttime <- (as.numeric(strptime(Aug16.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-16 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug17 <- NYbike[as.Date(NYbike$starttime) == "2019-08-17",]
Aug17.del <- Aug17[Aug17$start.station.id != Aug17$end.station.id & Aug17$tripduration > 60 & Aug17$tripduration < 3*60*60,c(2,4,8)]
Aug17.del$starttime <- (as.numeric(strptime(Aug17.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-17 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug18 <- NYbike[as.Date(NYbike$starttime) == "2019-08-18",]
Aug18.del <- Aug18[Aug18$start.station.id != Aug18$end.station.id & Aug18$tripduration > 60 & Aug18$tripduration < 3*60*60,c(2,4,8)]
Aug18.del$starttime <- (as.numeric(strptime(Aug18.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-18 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug19 <- NYbike[as.Date(NYbike$starttime) == "2019-08-19",]
Aug19.del <- Aug19[Aug19$start.station.id != Aug19$end.station.id & Aug19$tripduration > 60 & Aug19$tripduration < 3*60*60,c(2,4,8)]
Aug19.del$starttime <- (as.numeric(strptime(Aug19.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-19 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug20 <- NYbike[as.Date(NYbike$starttime) == "2019-08-20",]
Aug20.del <- Aug20[Aug20$start.station.id != Aug20$end.station.id & Aug20$tripduration > 60 & Aug20$tripduration < 3*60*60,c(2,4,8)]
Aug20.del$starttime <- (as.numeric(strptime(Aug20.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-20 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

Aug21 <- NYbike[as.Date(NYbike$starttime) == "2019-08-21",]
Aug21.del <- Aug21[Aug21$start.station.id != Aug21$end.station.id & Aug21$tripduration > 60 & Aug21$tripduration < 3*60*60,c(2,4,8)]
Aug21.del$starttime <- (as.numeric(strptime(Aug21.del$starttime,"%Y-%m-%d %H:%M:%S")) - as.numeric(strptime("2019-08-21 00:00:00","%Y-%m-%d %H:%M:%S")))/(60*60)

ride_num15 <- sapply(0:23,function(x) sum(floor(Aug15.del$starttime)==x))
ride_num16 <- sapply(0:23,function(x) sum(floor(Aug16.del$starttime)==x))
ride_num17 <- sapply(0:23,function(x) sum(floor(Aug17.del$starttime)==x))
ride_num18 <- sapply(0:23,function(x) sum(floor(Aug18.del$starttime)==x))
ride_num19 <- sapply(0:23,function(x) sum(floor(Aug19.del$starttime)==x))
ride_num20 <- sapply(0:23,function(x) sum(floor(Aug20.del$starttime)==x))
ride_num21 <- sapply(0:23,function(x) sum(floor(Aug21.del$starttime)==x))

NY_t <- c(ride_num15,ride_num16,ride_num17,ride_num18,ride_num19,ride_num20,ride_num21)
NY_overt <- data.frame(data = NY_t,T = 0.5:167.5)
data_t <- ggplot(aes(x=T,y=data),data = NY_overt) + geom_point() + geom_line()  + labs(x='Time',y="Total number of rides")
data_t <- data_t + scale_x_continuous(breaks = c(12,36,60,84,108,132,156),labels=c('August 15th','August 16th','August 17th','August 18th','August 19th','August 20th','August 21st'))

pdf(file='Bike_longt3.pdf', width=12, height=4)
data_t
dev.off()
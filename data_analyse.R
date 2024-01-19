#This document shows our analysis of the New York Citi Bike data.
#You need to download the pre-processed data NYbike.rda and the R file functions.R in your working path.
load("NYbike.rda")
source("functions.R")
library(ggplot2)

#obtain undirected networks
NYbike1 <- (NYbike1 + aperm(NYbike1,c(2,1,3)))/2

#We use the R function without diagonal data for analyse
NYbike.est <- PGD.panel2(NYbike1,2)

#plot the estimated Z (after a rotation)
theta <- (1 - 1/24) * pi
Z.rotate <- NYbike.est$Z %*% diag(c(1,-1)) %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate <- as.data.frame(Z.rotate)
names(Z.rotate) <- c("Z1", "Z2")
Z.rotate$borough <- factor(borough)
Z_est <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate,size = 0.8) + xlim(-3.1,4.9)
Z_est <- Z_est + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
Z_est <- Z_est + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
Z_est <- Z_est + theme(legend.position = c(0.8,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

#plot alpha
NY_t <- apply(NYbike1,3,sum)
alpha_est <- apply(NYbike.est$alpha,2,function(x) mean(as.vector(exp(outer(x,x,'+')))))
NY_overt <- data.frame(data = NY_t,exp_est = alpha_est,T = 0.5:23.5)
data_t <- ggplot(aes(x=T,y=data),data = NY_overt) + geom_point() + geom_line()  + labs(x='Time',y="Total number of rides")
data_t <- data_t + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))
alpha_est <- ggplot(aes(x=T,y=exp_est),data = NY_overt) + geom_point() + geom_line() + labs(x='Time',y=expression(paste("Mean of ",exp(hat(alpha[i]) + hat(alpha[j])))))
alpha_est <- alpha_est + scale_x_continuous(breaks = c(0,6,12,18,24),labels=c('0:00','6:00','12:00','18:00','24:00'))

#hist for all stations
NY_n <- data.frame(data = apply(NYbike1,1,sum))
data_n <- ggplot() + geom_histogram(data = NY_n,aes(x=data),stat = 'bin',bins = 100,fill='darkgreen',color='gray')
data_n <- data_n + labs(x = 'Total number of rides',y = 'Count of nodes')

#estimate using one-hour data
NYbike.time22 <- PGD.panel2(NYbike1[,,22],2)
theta <- (1-1/24)*pi
Z.rotate22 <- NYbike.time22$Z  %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
Z.rotate22 <- as.data.frame(Z.rotate22)
names(Z.rotate22) <- c("Z1", "Z2")
Z.rotate22$borough <- factor(borough)
Z_est22 <- ggplot() + geom_point(aes(x=Z1,y=Z2,col=borough,shape=borough),data = Z.rotate22,size = 0.8) + xlim(-3.5,4.5) + ylim(-4.15,4.15)
Z_est22 <- Z_est22 + scale_color_manual(values = c('#6495ED','#FF0000','#9A32CD')) 
Z_est22 <- Z_est22 + labs(x="1st component",y="2nd component") + scale_shape_manual(values = c(18,16,15))
Z_est22 <- Z_est22 + theme(legend.position = c(0.82,0.91)) + theme(legend.title=element_blank()) + theme(legend.key.size=unit(12,"pt"))

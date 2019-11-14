t0<-5
v0<-0.6
gamma<-4
sigma0<-0.1
tf<-200
alpha<-0.8




vt<-matrix(ncol=1,nrow=196,data=NA)
dt<-matrix(ncol=1,nrow=196,data=NA)
for( i in 1:1)
{
set.seed(floor(runif(1, min=0, max=1000)))
rn<-rnorm(n=200,m=0,sd=0.1)
dp<-cumsum(rn)
options(scipen=999)
temp<-c()
temp2<-c()

for( t in t0:tf )
{
	wf<-dp[t]
	wt<-dp[200]
	temp<-c(temp,exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma))+wt-((t-t0)/(tf-t0))*wf))
	#integrand <- function(t) {exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma))+wt-((t-t0)/(tf-t0))*wf)}
	#temp2<-c(temp2,integrate(integrand,t0,t)$value)
}
vt[,i]<-temp
dt[,i]<-cumsum(temp)
}
a<-cumsum(temp)
#plot()

vt_result<-apply(vt,1,mean)
#dt_result<-apply(dt,1,mean)

plot(5:200,vt_result,type="l",col="green",lwd=3,xlab="Exposure time(years)",ylab="Corrosion rate(mm/year)",xaxt="n", panel.first=grid(lwd=4),ylim=c(0,1))
axis(1, at = seq(0, 200, by = 20), las=1)
legend("topright",legend="1000 times average of simulation ",col="green",lty=1,lwd=2)
	   
plot(5:200,a[1:196],type="l",col="green",lwd=3,xlab="Exposure time(years)",ylab="Corrosion depth(mm)",xaxt="n", panel.first=grid(lwd=4))
axis(1, at = seq(0, 200, by = 20), las=1)







## 람다티
rambdaT<-matrix(ncol=1,nrow=196,data=NA)

for( i in 1:1)
{
set.seed(floor(runif(1, min=0, max=1000)))
rn<-rnorm(n=200,m=0,sd=0.1)
dp<-cumsum(rn)
options(scipen=999)
temp<-c()
temp2<-c()

for( t in t0:tf )
{
	wf<-dp[t]
	wt<-dp[200]
	temp<-c(temp,exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma))))
	#integrand <- function(t) {exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma))+wt-((t-t0)/(tf-t0))*wf)}
	#temp2<-c(temp2,integrate(integrand,t0,t)$value)
}
rambdaT[,i]<-temp
}

write.csv(rambdaT[,1],"C:/Users/user/Desktop/업무/lambdat.csv")



## 위너 프로세스
rnf<-rnorm(n=200,m=0,sd=0.1)
dpf<-cumsum(rnf)





### 람다티 적분

parameters <- read.csv("C:/Users/user/Desktop/업무/Parameter.csv")
set.seed(1234)
color <- round(runif(20,0,600))
i=1
temp<-c()
t0<-parameters[i,6]
v0<-parameters[i,2]
gamma<-parameters[i,4]	
alpha<-parameters[i,3]
for(j in t0:200)
{
	integrand <- function(t) {exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma)))}
	temp<-c(temp,integrate(integrand, lower =t0 , upper = j)$value)
}
plot(t0:200,temp,type="l", panel.first=grid(lwd=4),xlab="Exposure time(years)",ylab="Corrosion depth(mm)",main="Integral of RambdaT",ylim=c(0,35),col=colors()[color[i]])

for(i in 2:20)
{
	temp<-c()
	t0<-parameters[i,6]
	v0<-parameters[i,2]
	gamma<-parameters[i,4]	
	alpha<-parameters[i,3]
	for(j in t0:200)
	{
		integrand <- function(t) {exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma)))}
		temp<-c(temp,integrate(integrand, lower =t0 , upper = j)$value)
	}
	lines(t0:200,temp,col=colors()[color[i]])
}




## 람다티 적분 실제값 그래프 동시에

par(mfrow=c(4,5))
Data<- read.csv("C:/Users/user/Desktop/업무/시뮬레이션 데이터.csv")
absError<-matrix(ncol=3,nrow=20,data=NA)
absError<-as.data.frame(absError)
names(absError) <- c("50년","60년","70년")
rmse<-matrix(ncol=3,nrow=20,data=NA)
rmse<-as.data.frame(rmse)
names(rmse) <- c("27~50년","27~60년","27~70년")

library(hydroGOF)

for(i in 1:20)
{
	temp<-c()
	t0<-parameters[i,6]
	v0<-parameters[i,2]
	gamma<-parameters[i,4]	
	alpha<-parameters[i,3]
	for(j in t0:200)
	{
		integrand <- function(t) {exp(log(v0)+(alpha-1)*log(1+((t-t0)/gamma)))}
		temp<-c(temp,integrate(integrand, lower =t0 , upper = j)$value)
	}
	plot(t0:200,temp,type="l", panel.first=grid(lwd=4),xlab="Exposure time(years)",ylab="Corrosion depth(mm)",ylim=c(0,35),xlim=c(0,100),col="red",lwd=3)
	temp2<-c(0,Data[,i+1])
	lines(5:201,temp2,col="blue",lwd=3)
	legend("topleft",c("Integral of RambdaT", "Integral corrosion depth"),col=c("red", "blue"), lty=1, cex=0.8)
	
	absError[i,]<-abs(approx(t0:200, temp, c(50,60,70))$y-c(temp2[50-5+1],temp2[60-5+1],temp2[70-5+1]))
	rmse[i,1]<-rmse(approx(t0:200, temp, c(27:50))$y,temp2[(27-5+1):(50-5+1)])
	rmse[i,2]<-rmse(approx(t0:200, temp, c(27:60))$y,temp2[(27-5+1):(60-5+1)])
	rmse[i,3]<-rmse(approx(t0:200, temp, c(27:70))$y,temp2[(27-5+1):(70-5+1)])

}

write.csv(absError,"C:/Users/user/Desktop/업무/년도별에러.csv")
write.csv(rmse,"C:/Users/user/Desktop/업무/년도별rmse.csv")

####


## 데이터간 상관관계
Data<- read.csv("C:/Users/user/Desktop/업무/시뮬레이션 데이터.csv")
Data<-Data[,-1]

corMat<-cor(Data)


corMat2<-matrix(ncol=20,nrow=20,data=NA)
corMat2<-as.data.frame(corMat2)
names(corMat2) <- 0:19
for(i in 1:20){
corMat2[i,]<-as.numeric(gsub("V","",(names(sort(corMat[i,],decreasing=T)))))-1
}

write.csv(corMat2,"C:/Users/user/Desktop/업무/correlationmatrix.csv")


##


## contour plot
Data<- read.csv("C:/Users/user/Desktop/업무/결합개수 학습시간에 따른 rmse stacked method.csv")

Data<-Data[complete.cases(Data), ]
Data[,1]<-c(10,15,20,25,30,35)

Data2<-matrix(ncol=7,nrow=251,data=NA)
Data2<-as.data.frame(Data2)
names(Data2)<-c("기간",0:5)
Data2[,1]<-seq(10, 35, by = 0.1)


for(i in 1:6)
{
	Data2[,i+1]<-approx(Data[,1], Data[,i+1], c(seq(10, 35, by = 0.1)))$y
}



library(reshape2)

Data_melt <- melt(data = Data2, 
                           id.vars = "기간", 
                            measure.vars = c(2:7))


names(Data_melt) <- c("TrainingTime","TrainingSet","RMSE")


color <- rep(0:5,each=251)
Data_melt<-cbind(Data_melt,color)
library(plotly)

p <- plot_ly(Data_melt, x = ~TrainingTime, y = ~TrainingSet, z = ~RMSE, type = 'scatter3d', mode = 'lines',
        opacity = 1, line = list(width = 6, color = ~color,reverscale = FALSE))
screen=list(z=40, x=-60, y=0)
cloud(RMSE~결합데이터셋수+학습기간, Data_melt, panel.3d.cloud=panel.3dbars, col.facet='grey', 
      xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
      par.settings = list(axis.line = list(col = "transparent")),screen = list(z = -60, x = -30))


	  
library(ggplot2)
plot3 <- ggplot() +
         geom_tile(data = Data_melt, aes(TrainingSet, TrainingTime, RMSE, fill = equalSpace)) +
         geom_contour(color = “white”, alpha = 0.5) +
         theme_bw() +
         xlab("Weight (1,000lbs)") +
         ylab("Horsepower") +
         scale_fill_manual(values = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5", 
                                     "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a",
                                     "#543005", "#330000"),
                           name = "¼ Mi. Time (s)", breaks = breaks, labels = breaks)
						   
						   

dimnames(Data2)[[1]]<-Data2[,1]
		Data3<-Data2[,2:7]	
		par(mfrow=c(1,1))
image(as.numeric(dimnames(Data3)[[2]]),as.numeric(dimnames(Data3)[[1]]),t(as.matrix(Data3)),xlab="결합데이터셋수",ylab="학습기간(년)")

contour(as.numeric(dimnames(Data3)[[2]]),as.numeric(dimnames(Data3)[[1]]),t(as.matrix(Data3)), add = TRUE)
	
	
library(lattice)

	
contourplot(RMSE ~ TrainingSet * TrainingTime, data = Data_melt,
            cuts = 100, region = TRUE,
            xlab = "Training Set",
            ylab = "Training Time", labcex=5,
            main = "결합데이터셋수와 학습기간에 따른 RMSE stacked LSTM contour plot",at=seq(0, 15, length.out=15), col.regions = heat.colors(100)[1:length(heat.colors(100))])

	
	
## 히트맵

library("lattice")
 
## Example data
x <- seq(1,10, length.out=20)
y <- seq(1,10, length.out=20)
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)
 
## Try it out
par(mar=c(3,4,2,2))
levelplot(RMSE ~ 결합데이터셋수*학습기간, data=Data_melt  , xlab="결합데이셋수",ylab="학습기간" , col.regions = heat.colors(100)[1:length(heat.colors(100))], main="결합데이터셋수와 학습기간에 따른 RMSE",at=seq(0, 15, length.out=15))
 



##



#### LSTM 방법별 비교 그래프


Data<- read.csv("C:/Users/user/Desktop/업무/결합개수 학습시간에 따른 rmse.csv")

Data<-Data[complete.cases(Data), ]
Data[,1]<-c(10,15,20,25,30,35)
Data1<-matrix(ncol=7,nrow=251,data=NA)
Data1<-as.data.frame(Data1)
names(Data1)<-c("기간",0:5)
Data1[,1]<-seq(10, 35, by = 0.1)


for(i in 1:6)
{
	Data1[,i+1]<-approx(Data[,1], Data[,i+1], c(seq(10, 35, by = 0.1)))$y
}


Data2<- read.csv("C:/Users/user/Desktop/업무/결합개수 학습시간에 따른 rmse window method.csv")

Data2<-Data2[complete.cases(Data2), ]
Data2[,1]<-c(10,15,20,25,30,35)

Data22<-matrix(ncol=7,nrow=251,data=NA)
Data22<-as.data.frame(Data22)
names(Data22)<-c("기간",0:5)
Data22[,1]<-seq(10, 35, by = 0.1)


for(i in 1:6)
{
	Data22[,i+1]<-approx(Data2[,1], Data2[,i+1], c(seq(10, 35, by = 0.1)))$y
}


Data3<- read.csv("C:/Users/user/Desktop/업무/결합개수 학습시간에 따른 rmse stacked method.csv")

Data3<-Data3[complete.cases(Data3), ]
Data3[,1]<-c(10,15,20,25,30,35)

Data33<-matrix(ncol=7,nrow=251,data=NA)
Data33<-as.data.frame(Data33)
names(Data33)<-c("기간",0:5)
Data33[,1]<-seq(10, 35, by = 0.1)


for(i in 1:6)
{
	Data33[,i+1]<-approx(Data3[,1], Data3[,i+1], c(seq(10, 35, by = 0.1)))$y
}



dat2<-matrix(nrow=18,ncol=3,data=NA)
dat2<-as.data.frame(dat2)
names(dat2) <- c("Method","학습기간","RMSE")
dat2[,3]<-c(Data[,5],Data2[,5],Data3[,5])
dat2[,1]<-rep(c("LSTM for Regression","LSTM for Regression Using the Window Method","Stacked LSTM for Regression"),each=6)
dat2[,2]<-rep(c(5,10,15,20,25,30),3)


dat1 <- data.frame(
    Method = factor(c("LSTM for Regression","LSTM for Regression Using the Window Method","Stacked LSTM for Regression")),
    결합데이터셋수 = factor(c("0","0","0","1","1","1","2","2","2","3","3","3","4","4","4","5","5","5"), levels=c("0","1","2","3","4","5")),
    RMSE = c(7.113293, 6.274513,5.220718, 5.334092, 4.421036,4.827472,4.934456,3.411418,4.545326,4.935661,2.900692,4.161088,4.702985,3.272655,3.669254,4.444592,2.967262,3.773968)
)


library(ggplot2)
ggplot(data=dat1, aes(x=결합데이터셋수, y=RMSE, group=Method, colour=Method)) +
    geom_line(size=2.5) +
    geom_point(size=4.5)+ theme_bw()+ggtitle("5~27 years training")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(1,1),legend.justification=c(1,1),legend.box.just = c("top"),legend.background = element_rect(fill=alpha('blue', 0.1)))+ labs(x = "Training set")

number_ticks <- function(n) {function(limits) pretty(limits, n)}
ggplot(data=dat2, aes(x=학습기간, y=RMSE, group=Method, colour=Method)) +geom_point(size=1)+ggtitle("3 Training sets")+theme(plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=seq(0,30,1), limits=c(5,30))+geom_smooth(method = "loess", size = 1.5)+ theme_bw()+ggtitle("3 Training sets")+theme(plot.title = element_text(hjust = 0.5),legend.position=c(1,1),legend.justification=c(1,1),legend.box.just = c("top"),legend.background = element_rect(fill=alpha('blue', 0.1)))+ labs(x = "Training time")




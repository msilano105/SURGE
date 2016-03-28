


```{r}
tide_height = read.csv("~/SURGE/tide height data/2015NHA.txt", sep="", skip="10")
tide<-as.numeric(levels(height$f)[height$f])
surge<-as.numeric(levels(height$f2)[height$f2])
less.surge<-(tide+surge)
sorted_height<-tide_height[order(tide_height$yyyy.mm.dd),]
dates <- as.Date(sorted_height$yyyy.mm.dd)
date1 <- as.Date("2015-01-01")
date2 <- Sys.Date()-365
desired_rows <- which(dates >= date1 & dates <= date2)
height<-sorted_height[desired_rows, ]
library(chron)
tod<-chron(times=height$hh.mi.ssf)
dtod<-paste(dates[desired_rows] , tod)
x<- strptime(dtod, format="%Y-%m-%d %H:%M:%S")
y <-less.surge[desired_rows]

library(rjags)
nchain = 3

y.1=as.numeric(y)
data=list(n=length(desired_rows), z=y.1)





_______________________________________________________
## Random Walk Historical Surge Time Series

y=y

```{r}

RandomWalk = "
model{
#### Data Model
for(i in 1:n){
y[i] ~ dnorm(x[i],tau_obs)
}
#### Process Model
for(i in 2:n){
x[i]~dnorm(x[i-1],tau_add)
}
#### Priors
x[1] ~ dunif(-1,10)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
}
"
```


## MCMC
```{r}
data <- list(y=y,n=length(y),a_obs=1,r_obs=1,a_add=1,r_add=1)
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),tau_obs=5/var(y.samp))
}
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 100)

plot(jags.out)
```

```{r}
## after convergence
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000,
                            thin = 10)
time = 1:length(y)
time.rng = c(1,length(time)) ## adjust to zoom in and out
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
out <- as.matrix(jags.out)
ci <- apply(out[,3:ncol(out)],2,quantile,c(0.025,0.5,0.975))

plot(time,ci[2,],type='l',ylim=range(y,na.rm=TRUE),ylab="Surge Height",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ciEnvelope(time,ci[1,],ci[3,],col="lightBlue")
points(time,y,pch="+",cex=0.5)
```

```
layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])
plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])
```

write.csv(out, "Hist_Ran_Time_output.csv", row.names = T)
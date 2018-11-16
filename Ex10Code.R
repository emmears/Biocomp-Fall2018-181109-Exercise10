############ Tumor Cells Challenge from 11/12/18 ##################
#########1. In lecture, we used maximum likelihood and a likelihood ratio test####### 
#to complete a t-test. We can actually use a likelihood ratio test 
#to compare two models as long as one model is a subset of the other model. 
#For example, we can ask whether y is a hump-shaped vs. linear function of x 
#by comparing a quadratic (a+bx+cx2) vs. linear (a+bx) model. 
#Generate a script that evaluates which model is more appropriate 
#for the data in data.txt.

#likelihood ratio test

# load data
data=read.csv('data.txt',header=TRUE)

#create likelihood functions
linearMod<-function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  pred=B0+B1*x
  n=-sum(dnorm(x=y,mean=pred,sd=sigma,log=TRUE))
  
  return(n)
}

quadraticMod<-function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma=exp(p[2])
  
  pred=B0+B1*x+B2*x*x
  n=-sum(dnorm(x=y,mean=pred,sd=sigma,log=TRUE))
  
  return(n)
}

# estimate parameters
QuadraticGuess=c(30,35,25,1)
LinearGuess=c(80,40,1) 

fitlinear=optim(par=LinearGuess,fn=linearMod,x=data$x,y=data$y)
fitquadratic=optim(par=QuadraticGuess,fn=quadraticMod,x=data$x,y=data$y)

# run likelihood ratio test
teststat=2*(fitquadratic$value-fitlinear$value)

df=length(fitlinear$par)-length(fitquadratic$par)

1-pchisq(teststat,df)




###########2. A classic model of competition between two species was developed###### 
#by Lotka & Volterra. This model has two state variables described 
#by two di???erential equations:
  #dN1 dt = R1(1???N1??11???N2??12)N1 
  #dN2 dt = R2(1???N2??22???N1??21)N2 
#The criteria for coexistence of two species in the Lotka-Volterra competition 
#model is ??12 < ??11 and ??21 < ??22 
#Generate a script that uses three or more model simulations 
#to demonstratethe validity of these criteria for coexistence.

### Load the deSolve package and ggplot2 for plotting
library(deSolve)
library(ggplot2)

#model 
dnSim<-function(t,y,p){
  N=y[1]
  R=y[2]
  
  R1=p[1]
  N1=p[2]
  a11=p[3]
  a12=p[4]
  R2=p[5]
  N2=p[6]
  a21=p[7]
  a22=p[8]

  dN1dt = R1*(1-N1*a11-N2*a12)*N1 
  dN2dt = R2*(1-N2*a22-N1*a21)*N2
  
  return(list(c(dN1dt,dN2dt)))
}

# case 2
times=1:100
y0=c(0.1,0.1)
params2=c(0.5,10,2,0.5,10,0.5)
sim2=ode(y=y0,times=times,func=dnSim,parms=params2)
out2=data.frame(time=sim2[,1],dn1=sim2[,2],dn2=sim2[,3])
ggplot(out2,aes(x=time,y=dn1))+geom_line()+geom_line(data=out2,mapping=aes(x=time,y=dn2),col='red')+theme_classic()

# case 3
times=1:400
y0=c(0.05,0.3)
params3=c(0.5,10,0.5,0.5,10,0.5)
sim3=ode(y=y0,times=times,func=dnSim,parms=params3)
out3=data.frame(time=sim3[,1],dn1=sim3[,2],dn2=sim3[,3])
ggplot(out3,aes(x=time,y=dn1))+geom_line()+geom_line(data=out3,mapping=aes(x=time,y=dn2),col='red')+theme_classic()

# case 4
times=1:100
y0=c(0.1,0.1)
params4=c(0.5,10,0.5,0.5,10,2)
sim4=ode(y=y0,times=times,func=dnSim,parms=params4)
out4=data.frame(time=sim4[,1],dn1=sim4[,2],dn2=sim4[,3])
ggplot(out4,aes(x=time,y=dn1))+geom_line()+geom_line(data=out4,mapping=aes(x=time,y=dn2),col='red')+theme_classic()




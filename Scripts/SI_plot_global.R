#### Script belonging to the paper "Multiscale selection in spatially structured populations" by Hilje M. Doekes and Rutger Hermsen ####

### PLOT GLOBAL CHARACTERISTICS OF SI-MODEL SIMULATIONS (Fig 3c) ### 

## Data ##

global_default = read.table("Data/Default/countglobal.txt", sep="\t", header=T)  # Default
global_starthigh = read.table("Data/StartHigh/countglobal.txt", sep="\t", header=T)  # Start from higher transmissibility
global_mix = read.table("Data/Mix/countglobal.txt", sep="\t", header=T)  # Well-mixed


## Plot mean transmissibility ##

maxtime = 20000  # equivalent to 1000 generations; 1 generation is approximately 20 time steps (1/dS = 20)
plot(global_default$Time[global_default$Time<=maxtime],global_default$mean_inf[global_default$Time<=maxtime],type="l",lwd=3,col="steelblue3",xlim=c(0,maxtime+1),ylim=c(1,50),xaxs="i",yaxs="i",xlab="Time (time steps)",ylab="Mean transmissibility")
lines(global_starthigh$Time[global_starthigh$Time<=maxtime],global_starthigh$mean_inf[global_starthigh$Time<=maxtime],type="l",lwd=3,col="darkblue")
lines(global_mix$Time[global_mix$Time<=maxtime],global_mix$mean_inf[global_mix$Time<=maxtime],type="l",lwd=3,col="grey30")
legend("topleft",c("Structured","Structured","Well-Mixed"),lty=c(1,1,1),lwd=c(3,3),col=c("steelblue3","darkblue","grey30"))


## Plot number of susceptible and infected individuals ##

maxtime = 50000 # equivalent to 2500 generations

# Default (init = 5)
plot(global_default$Time[global_default$Time<=maxtime],global_default$Sus[global_default$Time<=maxtime],type="l",lwd=3,col="black",xlim=c(0,maxtime+1), ylim=c(0,500000),xaxs="i",yaxs="i",main="Default",xlab="Time (time steps)",ylab="Number of individuals")
lines(global_default$Time[global_default$Time<=maxtime], global_default$Infected[global_default$Time<=maxtime],type="l",lwd=3,col="red")
legend("topright",c("Susceptible","Infected"),lty=c(1,1),lwd=c(3,3),col=c("black","red"))

# Starting at high transmissibility (init = 12)
plot(global_starthigh$Time[global_starthigh$Time<=maxtime],global_starthigh$Sus[global_starthigh$Time<=maxtime],type="l",lwd=3,col="black",xlim=c(0,maxtime+1), ylim=c(0,500000),xaxs="i",yaxs="i",main="StartHigh",xlab="Time (time steps)",ylab="Number of individuals")
lines(global_starthigh$Time[global_starthigh$Time<=maxtime], global_starthigh$Infected[global_starthigh$Time<=maxtime],type="l",lwd=3,col="red")
legend("topright",c("Susceptible","Infected"),lty=c(1,1),lwd=c(3,3),col=c("black","red"))

# Well-mixed population
plot(global_mix$Time[global_mix$Time<=maxtime],global_mix$Sus[global_mix$Time<=maxtime],type="l",lwd=3,col="black",xlim=c(0,maxtime+1), ylim=c(0,500000),xaxs="i",yaxs="i",main="Well-mixed",xlab="Time (time steps)",ylab="Number of individuals")
lines(global_mix$Time[global_mix$Time<=maxtime], global_mix$Infected[global_mix$Time<=maxtime],type="l",lwd=3,col=rgb(1,0,0,0.7))
legend("topright",c("Susceptible","Infected"),lty=c(1,1),lwd=c(3,3),col=c("black","red"))

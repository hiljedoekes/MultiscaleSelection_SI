#### Script belonging to the paper "Multiscale selection in spatially structured populations" by Hilje M. Doekes and Rutger Hermsen ####

### PLOT SELECTION CURVES OF SI-MODEL SIMULATIONS (Fig 3d) ### 

### Data ###

tmin = 200000  # Time periode over which mean and standard error are calculated
tmax = 210000

## Default simulation (gamma = 0.05)

seldata = read.table("Data/Default/selection.txt", header=T)       # Read data
tvec = which((seldata$Time>=tmin & seldata$Time<=tmax),arr.ind=T)  # Construct vector with time points in data
ntimes = length(tvec)

# Initialise vectors
Sigma = c(seq(1.0,2.5,0.5), 3:50, seq(60,150,10))   # Values of local neighbourhood size tested (r in text)
mean_Slocal = vector("numeric")          # Mean value of S_local, will be calculated for varying values of Sigma/r
mean_Slocal_lowerr = vector("numeric")   # Confidence interval of S_local (+/- 2 SEM)
mean_Slocal_uperr = vector("numeric")
mean_Sinter = vector("numeric")          # Mean value of S_interlocal, will be calculated for varying values of Sigma/r
mean_Sinter_lowerr = vector("numeric")   # Confidence interval of S_interlocal (+/- 2 SEM)
mean_Sinter_uperr = vector("numeric")
mean_Stot = vector("numeric")            # Mean value of total selection (S)
Sat = vector("numeric")                  # Value of s(r), contribution of scale r to selection.

# For each length scale Sigma, calculate S_local, S_interlocal, s(r)
for (i in 1:length(Sigma)){
  curr_sigma = Sigma[i]
  tempname = paste("Sbelow_",format(round(curr_sigma,2),nsmall = 2),sep="")  # Name of column in dataframe seldata
  temp_mean_below = mean(seldata[tvec,tempname])          # Calculate mean S_local
  temp_sem = sqrt( var(seldata[tvec,tempname]) / ntimes ) # Calculate SEM of S_local
  mean_Slocal = c(mean_Slocal, temp_mean_below )          # Store in output vectors
  mean_Slocal_lowerr = c( mean_Slocal_lowerr, (temp_mean_below - 2*temp_sem) )
  mean_Slocal_uperr = c( mean_Slocal_uperr, (temp_mean_below + 2*temp_sem) )
  if(i>1){ # Approximate s(r) as the difference between Slocal for this and previous r values
    diffdata = seldata[tvec,tempname] - tempdata_prev
    diffSigma = Sigma[i]-Sigma[i-1]
    Sat = c( Sat, mean(diffdata)/diffSigma )
  }
  tempdata_prev = seldata[tvec,tempname]      # Store S_local values in dummy vector
  tempname = paste("Sabove_",format(round(curr_sigma,2),nsmall=2),sep="")    # Read S_interlocal values
  temp_mean_above = mean(seldata[tvec,tempname])          # Calculate mean of S_interlocal
  temp_sem = sqrt( var(seldata[tvec,tempname]) / ntimes ) # Calculate SEM of S_interlocal
  mean_Sinter = c(mean_Sinter, temp_mean_above)           # Store in output vectors
  mean_Sinter_lowerr = c( mean_Sinter_lowerr, (temp_mean_above - 2*temp_sem) )
  mean_Sinter_uperr = c( mean_Sinter_uperr, (temp_mean_above + 2*temp_sem) )
  mean_Stot = c(mean_Stot, (temp_mean_above + temp_mean_below))
}
# Construct data frames containing all output vectors
meanseldata_default = data.frame(Sigma, mean_Slocal, mean_Slocal_lowerr, mean_Slocal_uperr, mean_Sinter, mean_Sinter_lowerr, mean_Sinter_uperr, mean_Stot)
Sigma_diff = Sigma[2:length(Sigma)] - Sigma[1:(length(Sigma)-1)]
Sigma = Sigma[1:(length(Sigma)-1)] + Sigma_diff/2   # r values at which Sat is calculated
Sat_default = data.frame(Sigma,Sat)


## Simulation with larger patterns (gamma = 0.02)

seldata = read.table("Data/Gamma02/selection.txt", header=T)       # Read data
tvec = which((seldata$Time>=tmin & seldata$Time<=tmax),arr.ind=T)  # Construct vector with time points in data
ntimes = length(tvec)

# Initialise vectors
Sigma = c(seq(1.0,2.5,0.5), 3:50, seq(60,150,10))   # Values of local neighbourhood size tested (r in text)
mean_Slocal = vector("numeric")          # Mean value of S_local, will be calculated for varying values of Sigma/r
mean_Slocal_lowerr = vector("numeric")   # Confidence interval of S_local (+/- 2 SEM)
mean_Slocal_uperr = vector("numeric")
mean_Sinter = vector("numeric")          # Mean value of S_interlocal, will be calculated for varying values of Sigma/r
mean_Sinter_lowerr = vector("numeric")   # Confidence interval of S_interlocal (+/- 2 SEM)
mean_Sinter_uperr = vector("numeric")
mean_Stot = vector("numeric")            # Mean value of total selection (S)
Sat = vector("numeric")                  # Value of s(r), contribution of scale r to selection.

# For each length scale Sigma, calculate S_local, S_interlocal, s(r)
for (i in 1:length(Sigma)){
  curr_sigma = Sigma[i]
  tempname = paste("Sbelow_",format(round(curr_sigma,2),nsmall = 2),sep="")  # Name of column in dataframe seldata
  temp_mean_below = mean(seldata[tvec,tempname])          # Calculate mean S_local
  temp_sem = sqrt( var(seldata[tvec,tempname]) / ntimes ) # Calculate SEM of S_local
  mean_Slocal = c(mean_Slocal, temp_mean_below )          # Store in output vectors
  mean_Slocal_lowerr = c( mean_Slocal_lowerr, (temp_mean_below - 2*temp_sem) )
  mean_Slocal_uperr = c( mean_Slocal_uperr, (temp_mean_below + 2*temp_sem) )
  if(i>1){ # Approximate s(r) as the difference between Slocal for this and previous r values
    diffdata = seldata[tvec,tempname] - tempdata_prev
    diffSigma = Sigma[i]-Sigma[i-1]
    Sat = c( Sat, mean(diffdata)/diffSigma )
  }
  tempdata_prev = seldata[tvec,tempname]      # Store S_local values in dummy vector
  tempname = paste("Sabove_",format(round(curr_sigma,2),nsmall=2),sep="")    # Read S_interlocal values
  temp_mean_above = mean(seldata[tvec,tempname])          # Calculate mean of S_interlocal
  temp_sem = sqrt( var(seldata[tvec,tempname]) / ntimes ) # Calculate SEM of S_interlocal
  mean_Sinter = c(mean_Sinter, temp_mean_above)           # Store in output vectors
  mean_Sinter_lowerr = c( mean_Sinter_lowerr, (temp_mean_above - 2*temp_sem) )
  mean_Sinter_uperr = c( mean_Sinter_uperr, (temp_mean_above + 2*temp_sem) )
  mean_Stot = c(mean_Stot, (temp_mean_above + temp_mean_below))
}
# Construct data frame containing all output vectors
meanseldata_large = data.frame(Sigma, mean_Slocal, mean_Slocal_lowerr, mean_Slocal_uperr, mean_Sinter, mean_Sinter_lowerr, mean_Sinter_uperr, mean_Stot)
Sigma_diff = Sigma[2:length(Sigma)] - Sigma[1:(length(Sigma)-1)]
Sigma = Sigma[1:(length(Sigma)-1)] + Sigma_diff/2   # r values at which Sat is calculated
Sat_large = data.frame(Sigma,Sat)


### Plot S_local and S_interlocal ###

# Default simulation
plot(meanseldata_default$Sigma,meanseldata_default$mean_Slocal,type="l",lwd=3,col="red",bg="red", xlab="length scale (r)",ylab="Selection differential",xlim=c(0,100),ylim=c(-0.001,0.001))
polygon(c(meanseldata_default$Sigma,rev(meanseldata_default$Sigma)), c(meanseldata_default$mean_Slocal_lowerr,rev(meanseldata_default$mean_Slocal_uperr)), col=rgb(1,0,0,0.3), border=rgb(1,0,0,0.3) )
lines(meanseldata_default$Sigma,meanseldata_default$mean_Sinter,type="l",lwd=3,col="blue",bg="blue")
polygon(c(meanseldata_default$Sigma,rev(meanseldata_default$Sigma)), c(meanseldata_default$mean_Sinter_lowerr,rev(meanseldata_default$mean_Sinter_uperr)), col=rgb(0,0,1,0.3), border=rgb(0,0,1,0.3) )
lines(meanseldata_default$Sigma,meanseldata_default$mean_Stot,type="l",lwd=3,col="black",bg="black")
abline(h=0)
legend("topright",c(expression("S"["tot"]),expression("S"["local"]),expression("S"["interlocal"])),lty=c(1,1),lwd=c(2.5,2.5),col=c("black","red","blue"),cex = 1.5)

# Add lines for simulation with gamma = 0.02 (larger patterns)
lines(meanseldata_large$Sigma,meanseldata_large$mean_Slocal,type="l",lwd=2,pch=21,col="pink",bg="pink")
polygon(c(meanseldata_large$Sigma,rev(meanseldata_large$Sigma)), c(meanseldata_large$mean_Slocal_lowerr,rev(meanseldata_large$mean_Slocal_uperr)), col=rgb(1,0.75,0.8,0.3), border=rgb(1,0.75,0.8,0.3) )
lines(meanseldata_large$Sigma,meanseldata_large$mean_Sinter,type="l",lwd=3,pch=21,col="skyblue",bg="skyblue")
polygon(c(meanseldata_large$Sigma,rev(meanseldata_large$Sigma)), c(meanseldata_large$mean_Sinter_lowerr,rev(meanseldata_large$mean_Sinter_uperr)), col=rgb(0.53,0.81,0.92,0.3), border=rgb(0.53,0.81,0.92,0.3) )


### Plot s(r) ###

plot(Sat_default$Sigma, Sat_default$Sat, type="o", pch=21, col="black", bg="black", xlab="r",ylab="Contribution to selection s(r)")
abline(h=0)
lines(Sat_large$Sigma, Sat_large$Sat, type ="o", pch=21, col="grey50", bg="grey50")

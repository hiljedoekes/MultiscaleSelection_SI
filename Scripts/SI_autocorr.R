#### Script belonging to the paper "Multiscale selection in spatially structured populations" by Hilje M. Doekes and Rutger Hermsen ####

### Calculate critical scale of selection rC and autocorrelation coefficient for simulations with varying gamma ###
## Plot rC vs patch size estimated by autocorrelation (Fig 3e) - see end of script

## DATA Processing and calculations ##

# Location of data
inputfolders = c("Data/Gamma02","Data/Gamma021","Data/Gamma022","Data/Gamma023","Data/Gamma024","Data/Gamma025","Data/Gamma027","Data/Gamma03","Data/Gamma032","Data/Gamma035","Data/Gamma04","Data/Gamma045","Data/Default")

# Initialise output vectors
reprS = c(0.02,0.021,0.022,0.023,0.024,0.025,0.027,0.03,0.032,0.035,0.04,0.045,0.05)  # Values of gamma
Sat0 = vector("numeric")   # Values of critical scale of selection rC
expAC = vector("numeric")  # Values of exponent of autocorrelation (from exponential fit)

# Loop over all input folders
for (folder in inputfolders){
  print(folder)
  ##### Calculate Sat0 #####
  print("Calculating Sat0")
  seldata = read.table(paste(folder,"/selection.txt",sep=""), header=T)       # Read data
  tmin = 200000  # Time periode over which mean Slocal is calculated - for calculation of Sat0
  tmax = 210000
  tvec = which((seldata$Time>=tmin & seldata$Time<=tmax),arr.ind=T)  # Construct vector with time points
  Sigma = c(seq(1.0,2.5,0.5), 3:50, seq(60,150,10))   # Values of local neighbourhood sizes tested
  mean_Slocal_prev = 0
  for (i in 1:length(Sigma)){
    curr_sigma = Sigma[i]
    colname = paste("Sbelow_",format(round(curr_sigma,2),nsmall = 2),sep="")  # Name of column in dataframe seldata
    mean_Slocal = mean(seldata[tvec,colname]) # Calculate mean S_local
    # Test if s(r) is zero in this interval
    if (mean_Slocal < mean_Slocal_prev){
      Sat0 = c(Sat0, Sigma[i-1])        # We say rC is equal to the last r value for which s(r) was positive.
      break;
    }
    mean_Slocal_prev = mean_Slocal
  }
  ##### Fit exponential function to autocorrelation data (= calculate exponent) #####
  print("Calculating expAC")
  ACdata = read.table(paste(folder,"/autocorr_t200000.txt",sep=""),header=T,sep="\t")   # Read data
  ACdistances.df = read.table(paste(folder,"/autocorr_distances.txt",sep=""),header=F)  # Distances at which the autocorrelation was calculated
  ACdistances = ACdistances.df[,1]  # Vector instead of data frame
  mean_ac = apply(ACdata[,3:ncol(ACdata)],2,mean)  # Calculate the mean autocorr at given distances over several time steps. Do not include distance 0.
  ac.df = data.frame(ACdistances,mean_ac)
  ACmodel = nls(mean_ac ~ (1+a*exp(b*ACdistances)), data=ac.df, start=list(a=0.1,b=-0.1))  # Fit an exponential decay function to the autocorrelation data
  # plot with the fitted model for a sanity check
  plot(ac.df$ACdistances,ac.df$mean_ac,type="p",pch=21,cex=0.25,col="grey30",bg="grey30",main = folder,xlab="Distance (lattice points)",ylab="Autocorrelation")
  abline(h=1.0)
  lines(ac.df$ACdistances,predict(ACmodel),col="purple",lwd=1.5)
  # Add fitted exponent to output vector
  expAC = c(expAC, coef(ACmodel)[2])
}

# Construct data frame with all output and save
output.df = data.frame(reprS,Sat0,expAC)
#write.table(output.df,"Data/reprS_Sat0_expAC.txt")

## PLOT critical scale of selection vs patch size estimated from autocorrelation coefficient ##

data_AC = read.table("Data/reprS_Sat0_expAC.txt", header=T) # Data previously stored - avoid having to generate it all...
b = 1/(-data_AC$expAC)
plot(data_AC$Sat0,b,type="p",xlim=c(0,17),ylim=c(0,17),pch=21,col="black",bg="black",xlab="Critical scale of selection rC", ylab="Patch size as estimated from autocorrelation",xaxs="i",yaxs="i")
abline(lm(b ~ 0 + data_AC$Sat0),lty=3,lwd=2)  # Add regression line



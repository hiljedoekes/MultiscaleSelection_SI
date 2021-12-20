#### Script belonging to the paper "Multiscale selection in spatially structured populations" by Hilje M. Doekes and Rutger Hermsen ####

### PLOT STILL OF SIMULATION LATTICE (Fig 3b) ### 

## External Function ##

#This function creates a color scale for use with e.g. the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
#Depending on the orientation, x- or y-limits may be defined that
#are different from the z-limits and will reduce the range of
#colors displayed.
# Taken from: http://menugget.blogspot.nl/2011/08/adding-scale-to-image-plot.html#more

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

## Data ##
nr = 1024   # number of rows in simulation lattice
nc = 1024   # number of columns
#plane = read.csv("Data/Default/saveplane200000.txt",sep=" ", header=F, skip=4)  # Default (gamma = 0.05)
plane = read.csv("Data/Gamma02/saveplane200000.txt",sep=" ", header=F, skip=4)  # Larger patterns (gamma = 0.02)

## Color scale for transmissibility ##
colorscale_inf = rgb(c(1,seq(0,1,length=256)),c(1,seq(0,0,length=256)),c(1,seq(1,0,length=256)))

## Construct matrices of individuals, susceptibles only, and transmissibilities of infected individuals ##
Ind_mat = matrix(plane$V1,nrow = nr,ncol = nc)    # Individuals
S_mat = ifelse(Ind_mat==2,1,0)                # Susceptibles
T_mat = matrix(plane$V2,nrow = nr, ncol = nc)     # Transmissibilities

## Construct data frames with coordinates of susceptible and infected individuals ##
S = as.data.frame((which(S_mat!=0, arr.ind = T)))  # Susceptibles
I = as.data.frame((which(T_mat!=0, arr.ind = T)))  # Infected
## Add transmissibilities to I dataframe
transvals = plane$V2[plane$V1==3]
I$transval = transvals

## Select which part of the lattice is plotted ##
xmin = 75    # 150 - 300 for default, # 75 - 225 for large patterns
xmax = 225
ymin = 850
ymax = 1000

## Set color scale for transmissibility ##
min_trans = 6
max_trans = 12
colorscale_trans = rgb(c(1,seq(0,1,length=256)),c(1,seq(0,0,length=256)),c(1,seq(1,0,length=256)))

## Plot ##
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1))
par(mar=c(3,3,3,1),xaxt='n',yaxt='n')
transbreaks=c(0,seq(min_trans,max_trans,length=257))
# Susceptibles
plot(S$row[S$row>=xmin & S$row<=xmax & S$col>=ymin & S$col<=ymax],S$col[S$row>=xmin & S$row<=xmax & S$col>=ymin & S$col<=ymax],xaxs="i",yaxs="i",xlab="",ylab="",col="gray50",pch=19,cex=0.35)
# Infected, coloured by transmissibility
cutvec = cut(I$transval[I$row>=xmin & I$row<=xmax &I$col>=ymin & I$col<=ymax],transbreaks,labels=F)
colvec = colorscale_trans[2:257]
points(I$row[I$row>=xmin & I$row<=xmax &I$col>=ymin & I$col<=ymax], I$col[I$row>=xmin & I$row<=xmax &I$col>=ymin & I$col<=ymax],col=colvec[cutvec],pch=19,cex=0.35)
image.scale(T_mat[T_mat!=0], col=colorscale_inf[2:257], breaks=transbreaks[2:258], horiz=F)

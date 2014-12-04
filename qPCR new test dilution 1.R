# New Test dilution qPCR 1

DataJPMt1 <- read.table("qPCR new test dilution 1.txt") 
# DataJPMt1[,1] <- paste(DataJPMt1[,1],DataJPMt1[,2],sep = "_")
# DataJPMt1[,2] <- NULL

DIL1raw <- DataJPMt1[,2:41]  # 2 to 41 is te data

plot(1:40,DIL1raw[1,])
for (i in 6:10){
    points(1:40,DIL1raw[i,], col=i)
}

# one way to optimize is to find the region where the data is important
# or filtrate the background

DIL1max <- apply(X = DIL1raw,MARGIN = 2,FUN = max)


# plot(1:40,DIL1max)

dydxmax <- vector()
for(i in 1:(length(DIL1max)-1)){
    dydxmax[i] <- (DIL1max[(i + 1)]-DIL1max[i])
}

# plot(1:39,dydxmax)

maxinfpoint <- which.max(x = dydxmax)

dydxmax <- dydxmax[1:maxinfpoint]
# log of dydx until inflection point


minlogdydx <- log10(min(dydxmax[dydxmax > 0]))
log_dydxmax <- vector()

log_dydxmax <- rapply(as.list(dydxmax),function(n){  # in raw data for more curves
    if(n<=0){
        n <- minlogdydx  # if the value is negative, put the log of the min
    } else {
        n <- log10(n)
    }
})

# plot(1:maxinfpoint,log_dydxmax[1:maxinfpoint])

dydxmin <- maxinfpoint - 3

while(summary(lm(dydxmin:maxinfpoint ~ log_dydxmax[dydxmin:maxinfpoint]))$r.squared >= 0.9 ){  # 0.8 to be more carful
    dydxmin <- dydxmin - 1
}

# this was to find the background level from all the data. The user can define it, but is better to calculate it
# if the user define it the algorithm should tell him what is the calculated one.


####################################################

source("qPCRV4FUNCTIONS.R")

#### FUNCTION 1  # parameters for the analysis
Parameters(rowdata = DIL1raw)

Par$Start <- dydxmin
# this was to find the background level from all the data. The user can define it, but is better to calculate it
# if the user define it the algorithm should tell him what is the calculated one.

#### FUNCTION 2  # calculate the derivate of the curve to find the inflection point
Inflection_point(DIL1raw)

## plot(2:40,dydx[3,])

#### Object 1  # We record the names of the well and the inflection points in this Object 1
inf_curves[1:Par$curves,1] <- DataJPMt1[1:Par$curves,1]        # Names of the first column


# add the max of the derivate to each curve
for (curv in 1:Par$curves){
    if(which.max(dydx[curv,]) == 1) {
        inf_curves[curv,2] <- 5
    } else {
        inf_curves[curv,2] <- which.max(dydx[curv,])
    }
}

#### FUNCTION 3

log_Data(rawData = DIL1raw,Name = "DIL1")

DIL1_log_Data[1,]

# TEST
# lm(21:25 ~ DIL1raw_log_Data[1,21:25])



# the negatives value take the min, but there are some values in the begining that are background

# Now from this points, I need to find the best line from cycle 1 to 28 (of course from 20 to 28 but I am thinking in general)
# So all the possible combination of data from 1 to 40 is given by it = ((End-dif+Start)*(End-dif))/2 = 666
# The interactions that will be used here will be less, from 1 to 28 so only 300

# this also tell you the combination
# points <- list(1:40,40:1)
# p1_p2 <- expand.grid(points)


#### FUNCTION 4

Init_Value(DIL1_log_Data)

    for (cur in 1:Par$curves){ # first loop to do the same for each curve.
        # matrix to record the R^2 of combinations
        bestline <- matrix(data = NA, nrow = Values$it, ncol = 5) #should be out of the loop and the nrow in this case is the max iteraction
        count <- 1 # From here could be possible to add the interaction values in relation to each "cur"
        
        #points <- list(1:inf_curves[cur,2],inf_curves[cur,2]:1)
        #p1_p2 <- expand.grid(points)
        #p1_p2[1,1]
        #p1_p2[1,2]    
        #we can use this for any apply but for the moment I don't know and this works
        #cur<-1
        #log_Data <- Neilraw_log_Data
        #i<-21
        #j<-25
        
        
        for (i in Par$Start:inf_curves[cur,2]) {        # in this case is inflection point
            for (j in inf_curves[cur,2]:Par$Start) {
                if ( (i < j) & (j-i >= Par$dif) ) {
                    Temp_Data_lm <- lm((i:j) ~ (DIL1_log_Data[cur,(i:j)])) # known Y ~ known X, to get the right slope
                    bestline[count,1:2] <- (c(i,j))# 3 a, 4 b
                    bestline[count,3] <- summary(Temp_Data_lm)$r.squared # 5 R^2
                    bestline[count,4] <- Temp_Data_lm$coefficients[[2]] # 6 slope(x~y)
                    bestline[count,5] <- Temp_Data_lm$coefficients[[1]] # 7 CtY0 cycle where R is 0
                    count <- count + 1
                } 
            }
        }
        bestline <- bestline[(bestline[,4] < Par$slope_min) & (bestline[,4] > Par$slope_max),] 
        bestline <- bestline[complete.cases(bestline),]
        bestline <- as.matrix(bestline)
        if(ncol(bestline) == 1) bestline <- t(bestline)
        
        if(nrow(bestline)!=0){       # To skip the error
            inf_curves[cur,3:7]  <- bestline[which.max(bestline[,3]),] #copy the best R^2 line to info curves
            inf_curves[cur,8]    <- (10^(1/inf_curves[cur,6])-1)+1 # 8 convert the slope in efficiency
            inf_curves[cur,9]    <- (inf_curves[cur,8]-1)*100 # 9 slope in %, probably unecesary
            inf_curves[cur,10]   <- 10^(-inf_curves[cur,7]/inf_curves[cur,6]) 
            # 10 Ri, is the anti log of C when Y is 0
            # or the level of emmission at the cycle 0
        } else {
            inf_curves[cur,3:10] <<- NA
        }
    }



# system.time(Init_Value(Neilraw_log_Data))
#   user  system elapsed 
# 269.52    0.51  285.98 
#  37.72    0.00   37.92 

# With this we can get the best best line accordig to the R^2 but that also is in the slope range
# this says that the data from 21 to 25 is linear or exponential
## inf_curves

#very cool staf can be measured here
# hist(inf_curves[,6], breaks = 50)
# rug(inf_curves[,6])


# and at the end we just isolate the exponential part of the curve:
# plot(21:25,Neilraw[1,21:25])
# or the linear
# plot(21:25,Neilraw_log_Data[1,21:25])

# With this we have a good stimation of the initial value Ri, then we just need to play
# in how to combine them and for that I have a few ideas but this is very simple.

# Next step, optimize this as a function and put this in a web site.

inf_curves

plot(1:50,Rod1raw[1,])
test1 <- vector()
for(i in 1:50){
    test1[i] <- inf_curves[1,10]*(inf_curves[1,8]^i)
}
points(1:50,test1,col=2)

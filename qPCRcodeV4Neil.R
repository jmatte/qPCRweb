# This is a real data. with 40 cycles
DataNeil <- read.table("qPCRDataNeil.txt") 
DataNeil[,1] <- paste(DataNeil[,1],DataNeil[,2],sep = "_")
DataNeil[,2] <- NULL

Neilraw <- DataNeil[,2:41]  # 2 to 41 is te data
## plot(1:40,Neilraw[3,])


# one way to optimize is to find the region where the data is important
# or filtrate the background

Neilmax <- apply(X = Neilraw,MARGIN = 2,FUN = max)


plot(1:40,Neilmax)

dydxmax <- vector()
for(i in 1:(length(Neilmax)-1)){
    dydxmax[i] <- (Neilmax[(i + 1)]-Neilmax[i])
}

# plot(1:39,dydxmax)

maxinfpoint <- which.max(x = dydxmax)

dydxmax <- dydxmax[1:maxinfpoint]
# log of dydx until inflection point


minlogdydx <- log10(min(dydxmax[dydxmax > 0]))
log_dydxmax <- vector()

log_dydxmax <- rapply(as.list(dydxmax),function(n){  # in raw data for more curves
    if(n<0){
        n <- minlogdydx  # if the value is negative, put the log of the min
    } else {
        n <- log10(n)
    }
})

# plot(1:maxinfpoint,log_dydxmax[1:maxinfpoint])

dydxmin <- maxinfpoint - 3

while(summary(lm(dydxmin:maxinfpoint ~ log_dydxmax[dydxmin:maxinfpoint]))$r.squared >= 0.9 ){
    dydxmin <- dydxmin - 1
}

# this was to find the background level from all the data. The user can define it, but is better to calculate it
# if the user define it the algorithm should tell him what is the calculated one.


####################################################

source("qPCRV4FUNCTIONS.R")

#### FUNCTION 1  # parameters for the analysis
Parameters(rowdata = Neilraw)

Par$Start <- dydxmin
# this was to find the background level from all the data. The user can define it, but is better to calculate it
# if the user define it the algorithm should tell him what is the calculated one.

#### FUNCTION 2  # calculate the derivate of the curve to find the inflection point
Inflection_point(Neilraw)

## plot(2:40,dydx[3,])

#### Object 1  # We record the names of the well and the inflection points in this Object 1
inf_curves[1:Par$curves,1] <- DataNeil[1:Par$curves,1]        # Names of the first column


# add the max of the derivate to each curve
for (curv in 1:Par$curves){
    if(which.max(dydx[curv,]) == 1) {
        inf_curves[curv,2] <- 5
    } else {
        inf_curves[curv,2] <- which.max(dydx[curv,])
    }
}

## now we know that the PCR is exponential from cycle 1 to 28 or inflection point, but the
## machine only detect data over certain minimum. So it could be 20-28 or 21-28 or 10-28 or any combination...


## now we have a bunch of parameters and we need to test which combination is the best for the curve
## To do that I would still use the log, because it is easy to find the best line in a linear data set...

#### FUNCTION 3

log_Data(rawData = Neilraw,Name = "Neilraw")

## plot(1:40,Neilraw_log_Data[3,])

# the negatives value take the min, but there are some values in the begining that are background

# Now from this points, I need to find the best line from cycle 1 to 28 (of course from 20 to 28 but I am thinking in general)
# So all the possible combination of data from 1 to 40 is given by it = ((End-dif+Start)*(End-dif))/2 = 666
# The interactions that will be used here will be less, from 1 to 28 so only 300

# this also tell you the combination
# points <- list(1:40,40:1)
# p1_p2 <- expand.grid(points)


#### FUNCTION 4

Init_Value(Neilraw_log_Data)

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

plot(1:40,Neilraw[144,])
test144 <- vector()
for(i in 1:40){
    test144[i] <- inf_curves[144,10]*(inf_curves[144,8]^i)
}
points(1:40,test144,col=2)





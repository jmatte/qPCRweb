# Test Rodrigo 1

DataRod1 <- read.table("RodrigoClean_1.txt") 
DataRod1[,1] <- paste(DataRod1[,1],DataRod1[,2],sep = "_")
DataRod1[,2] <- NULL

Rod1raw <- DataRod1[,2:41]  # 2 to 41 is te data

plot(1:40,Rod1raw[1,])
for (i in 1:24){
    points(1:40,Rod1raw[i,], col=i)
}

# one way to optimize is to find the region where the data is important
# or filtrate the background

Rod1max <- apply(X = Rod1raw,MARGIN = 2,FUN = max)


# plot(1:40,Rod1max)

dydxmax <- vector()
for(i in 1:(length(Rod1max)-1)){
    dydxmax[i] <- (Rod1max[(i + 1)]-Rod1max[i])
}

# plot(1:20,dydxmax)

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
Parameters(rowdata = Rod1raw)

Par$Start <- dydxmin
# this was to find the background level from all the data. The user can define it, but is better to calculate it
# if the user define it the algorithm should tell him what is the calculated one.

#### FUNCTION 2  # calculate the derivate of the curve to find the inflection point
Inflection_point(Rod1raw)

## plot(2:40,dydx[3,])

#### Object 1  # We record the names of the well and the inflection points in this Object 1
inf_curves[1:Par$curves,1] <- DataRod1[1:Par$curves,1]        # Names of the first column


# add the max of the derivate to each curve
for (curv in 1:Par$curves){
    if(which.max(dydx[curv,]) == 1) {
        inf_curves[curv,2] <- 5
    } else {
        inf_curves[curv,2] <- which.max(dydx[curv,])
    }
}

#### FUNCTION 3

log_Data(rawData = Rod1raw,Name = "Rod1")

Rod1_log_Data[1,]

# TEST
# lm(21:25 ~ Neilraw_log_Data[1,21:25])



# the negatives value take the min, but there are some values in the begining that are background

# Now from this points, I need to find the best line from cycle 1 to 28 (of course from 20 to 28 but I am thinking in general)
# So all the possible combination of data from 1 to 40 is given by it = ((End-dif+Start)*(End-dif))/2 = 666
# The interactions that will be used here will be less, from 1 to 28 so only 300

# this also tell you the combination
# points <- list(1:40,40:1)
# p1_p2 <- expand.grid(points)


#### FUNCTION 4

Init_Value(Rod1_log_Data)

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

plot(1:40,Rod1raw[1,])
test1 <- vector()
for(i in 1:40){
    test1[i] <- inf_curves[1,10]*(inf_curves[1,8]^i)
}
points(1:40,test1,col=2)

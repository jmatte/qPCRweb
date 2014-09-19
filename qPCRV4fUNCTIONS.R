############ FUNCTION 1 ############################################
# Parameters function, determine all the parameters for the analysis

Par <- list(Start = NULL, End = NULL, dif = NULL, curves = NULL, Eff_min = NULL, Eff_max = NULL, slope_min = NULL, slope_max = NULL)
Values <- list(it = NULL, Cycles = NULL)

Parameters <- function(Start = 1, End = length(rowdata), dif = 4, rowdata=NULL, Eff_min = 0.6, Eff_max = 1.2){
    Start  <- Start
    End <- End
    dif    <- dif # min number of points used in the curve
    curves <- nrow(rowdata)
    Eff_min <- Eff_min
    Eff_max <- Eff_max
    slope_min <- 1/log10(Eff_min+1)
    slope_max <- 1/log10(Eff_max+1)
    Par <<- list(Start = Start, End = End, dif = dif, curves = curves, 
                 Eff_min = Eff_min, Eff_max = Eff_max, slope_min = slope_min, slope_max = slope_max)
    Values <<- list(it = ((End-dif+Start)*(End-dif))/2, Cycles = Start:End)
}

#####################################################################

############ FUNCTION 2 ############################################
## We calculate the derivative of all the curves, so we need a matrix length x-1 with same rows as curves


## may be there is an apply to do that, but... and also should be as a function
Inflection_point <- function(rawData = NULL){
    dydx <<- matrix(nrow = Par$curves, ncol = Par$End - 1)
    for(curv in 1:Par$curves){
        for(i in 1:(Par$End - 1)){
            dydx[curv,i] <<- (rawData[curv,(i + 1)]-rawData[curv,i])
        }
    }
}

#####################################################################

############ New Object 1 ###########################################
inf_curves <- data.frame(Well = NA, inf_point = NA, a = NA, b = NA, "R^2" = NA, 
                         "slope(x~y)" = NA, CtY0 = NA, Eff = NA, "Eff'%'" = NA, Ri = NA)         # the second value is the inlfection point

#####################################################################


############ FUNCTION 3 ############################################
## very optimizable function, but this in theory will put in log any data that I have.

log_Data <- function(rawData = NULL,Name = "logData"){
    namelogdata <- paste(Name,"log_Data",sep = "_")
    log_Data <- matrix(nrow = Par$curves, ncol = Par$End)
    for (cur in 1:Par$curves){
        linetested <- rawData[cur,] # to put in the loop the curve that is tested
        minval <- min(linetested[linetested > 0])   # to check the min value on each curve 
        log_Data[cur,] <- rapply(as.list(rawData[cur,]),function(n){  # in raw data for more curves
            if(n<0){
                n <- log10(minval)  # if the value is negative, put the log of the min
            } else {
                n <- log10(n)
            }
        })
    }
    assign(namelogdata,log_Data,envir = .GlobalEnv)
}

#####################################################################


############ FUNCTION 4 ############################################
# TEST
# log_Data = Rod1_log_Data
# cur = 2
# rm(log_Data)
# rm(cur)
# rm(bestline)
# rm(count)


Init_Value <- function (log_Data = NULL){
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
                    Temp_Data_lm <- lm((i:j) ~ (log_Data[cur,(i:j)])) # known Y ~ known X, to get the right slope
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
            inf_curves[cur,3:7]  <<- bestline[which.max(bestline[,3]),] #copy the best R^2 line to info curves
            inf_curves[cur,8]    <<- (10^(1/inf_curves[cur,6])-1)+1 # 8 convert the slope in efficiency
            inf_curves[cur,9]    <<- (inf_curves[cur,8]-1)*100 # 9 slope in %, probably unecesary
            inf_curves[cur,10]   <<- 10^(-inf_curves[cur,7]/inf_curves[cur,6]) 
                                                        # 10 Ri, is the anti log of C when Y is 0
                                                        # or the level of emmission at the cycle 0
        } else {
            inf_curves[cur,3:10] <<- NA
        }
    }
}

#####################################################################



# source("qPCRcodeV4Neil.R")

# data from Neil

head(inf_curves)

Samples <- read.table("Neilsamples.txt", fill = T)
Genes <- read.table("Neilgenes.txt", fill = T, header = T)


p1s <- read.table("Neilplate1samples.txt", fill = T, header = T, nrows = 8) 
p1s[,1] <- NULL
p1s[,12] <- NA
colnames(p1s) <- 1:12
rownames(p1s) <- LETTERS[1:8]

p1g <- read.table("Neilplate1genes.txt", fill = T, header = T, nrows = 8) 
p1g[,1] <- NULL
p1g[,12] <- NA
colnames(p1g) <- 1:12
rownames(p1g) <- LETTERS[1:8]

p2s <- read.table("Neilplate2samples.txt", fill = T, header = T, nrows = 8) 
p2s[,1] <- NULL
p2s[,12] <- NA
colnames(p2s) <- 1:12
rownames(p2s) <- LETTERS[1:8]

p2g <- read.table("Neilplate2genes.txt", fill = T, header = T, nrows = 8) 
p2g[,1] <- NULL
p2g[,12] <- NA
colnames(p2g) <- 1:12
rownames(p2g) <- LETTERS[1:8]

# check technical replicates
# Put the gene and sample in a new column to check as a factor and search for tech rep
# put together the tech rep as list to plot and check and average.
plates <- inf_curves[,1]

tmp <- strsplit(inf_curves$Well,'_')
tmp <- do.call(rbind,tmp)

paste0('p',tmp[,1],'g')
tmp[,2]

find_gene <- function(plate,well){
    plate_ <- get(plate,envir=.GlobalEnv)
    x <- regmatches(well,regexpr('\\d+',well))
    y <- regmatches(well,regexpr('[a-zA-Z]',well))
    as.character(plate_[which(rownames(plate_) == y),as.numeric(x)])
}


inf_curves$gene <- mapply(find_gene,plate=paste0('p',tmp[,1],'g'),well=tmp[,2])
inf_curves$sample <- mapply(find_gene,plate=paste0('p',tmp[,1],'s'),well=tmp[,2])
inf_curves$code <- paste(inf_curves$gene,inf_curves$sample,sep = "_") # we culd use a code.. match(inf_curves$gene[1],Genes[,1])
    
head(inf_curves)

#technical replicates

tRepAv <- tapply(inf_curves$Ri,inf_curves$code,mean)
tRepSD <- tapply(inf_curves$Ri,inf_curves$code,sd)

tRepAv[2]
tRepSD[1]

levels(inf_curves$code)

#Samples 3G, 3N, 3B


# Plate
# Sample 1 (3G) Control
# Sample 2,3 (3N,3B) Treatment
# RefGene ParA
# TarGene SodA, OtsA, Rpon2, aceA, LtpD, pspA

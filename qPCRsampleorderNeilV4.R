# source("qPCRcodeV4Neil.R")

# data from Neil

head(inf_curves)

Samples <- read.table("Neilsamples.txt", fill = T)
# samples could be control or treatment, and will have biological replicates
# also the initial concentration of some samples could be known, but by default
# control initial concentration will be 1 and treatment X. T/C

Genes <- read.table("Neilgenes.txt", fill = T, header = T)
# genes could be reference or target. the reference gene had a value 1 and the target is x
# And the calculation of the target gene is in base of the reference t/r


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
tRepall <- tapply(inf_curves$Ri,inf_curves$code,function(x) x)
tRepDIF <- tapply(inf_curves$Ri,inf_curves$code,function(x) max(x)-min(x))

ref <- as.character(Genes[Genes$Condition == 'Reference',1])

comb <- unique(regmatches(names(tRepAv),regexpr('\\d+[_]\\d+',names(tRepAv))))
normalised <- sapply(comb,function(x){
    sapply(USE.NAMES=F,names(tRepAv)[grepl(x,names(tRepAv)) & !grepl(ref,names(tRepAv))],function(y){
        tRepAv[y]/tRepAv[paste0(ref,'_',x)]
    })  
})

normalised <- unlist(normalised)



gen <- "SodA"
biorep_c = '1_1'
biorep_t = '2_1'
treatment <- sapply(simplify='array',as.character(Genes[Genes$Condition == 'Target',1]),function(gen){
    selected_gen <- normalised[grepl(gen,names(normalised))]
    sapply(comb[grepl(paste0(as.character(which(Samples$Group == 'Control')),'[_]'),comb)],function(biorep_c){
        sapply(comb[!grepl(paste0(as.character(which(Samples$Group == 'Control')),'[_]'),comb)],function(biorep_t){
            normalised[grepl(paste0(gen,'_',biorep_t),names(normalised))]/normalised[grepl(paste0(gen,'_',biorep_c),names(normalised))]
        })
    })
})

results_mean <- apply(treatment,3,function(gen){
    sapply(unique(gsub('(?<=[_]).+','',rownames(gen),perl = T)),function(t){
        mean(unlist(gen[grepl(t,rownames(gen)),]),na.rm = T)
    })
})

results_mean <- reshape2::melt(results_mean)

results_sd <- apply(treatment,3,function(gen){
    sapply(unique(gsub('(?<=[_]).+','',rownames(gen),perl = T)),function(t){
        sd(unlist(gen[grepl(t,rownames(gen)),]),na.rm = T)
    })
})

results_sd <- reshape2::melt(results_sd)
to_plot <- cbind(results_mean,results_sd$value)
colnames(to_plot) <- c('Sample','Gene','mean','sd')

library(ggplot2)
ggplot(to_plot,aes(x=Gene,y=log2(mean))) + geom_bar(aes(fill = Sample), position = "dodge", stat="identity") +
    geom_errorbar(aes(ymin=log2(mean-sd), ymax=log2(mean+sd),group=Sample),width=.1,position=position_dodge(.9))


summary(as.factor(inf_curves$code))

# cy for tech rep can be CtY0


#Samples 3G, 3N, 3B


# Plate
# Sample 1 (3G) Control
# Sample 2,3 (3N,3B) Treatment
# RefGene ParA
# TarGene SodA, OtsA, Rpon2, aceA, LtpD, pspA

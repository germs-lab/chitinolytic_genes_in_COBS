get_count <- function(gename, dat, final,meta){
#make sum
su = c()
for (i in 2:dim(dat)[2]){
	su[i] <- sum(dat[,i])
}

#cal average for each crop
foo <- vector(mode="list", length=3)
names(foo) <- c("CC", "P", "FP")
foo["CC"] <- 0
foo["P"] <- 0
foo["FP"] <- 0

for (i in 2:dim(dat)[2]){
	crp = meta[meta$rast_file == names(dat)[i],][1] 
	if (foo[crp[1,1]] == 0){
		foo[crp[1,1]] <- su[i]
	}else{
		temp <- as.numeric(foo[crp[1,1]])
		temp <- temp + su[i]
		foo[crp[1,1]]  <- temp
	}
}
temp <- as.numeric(foo["CC"])
temp <- temp/sum(meta$Crop == "CC")
foo["CC"] <- temp

temp <- as.numeric(foo["P"])
temp <- temp/sum(meta$Crop == "P")
foo["P"] <- temp

temp <- as.numeric(foo["FP"])
temp <- temp/sum(meta$Crop == "PF")
foo["FP"] <- temp

gene <- c(gename,gename,gename)
va <- c("CC","P","FP")
value <- c(as.numeric(foo["CC"]),as.numeric(foo["P"]),as.numeric(foo["FP"]))
cfinal = data.frame(gene,va,value)
final = rbind(final,cfinal)

#cal average for each aggregate
agfoo <- vector(mode="list", length=5)
names(agfoo) <- c("LM", "MM", "SM", "Micro", "WS")
agfoo["LM"] <- 0
agfoo["MM"] <- 0
agfoo["SM"] <- 0
agfoo["Micro"] <- 0
agfoo["WS"] <- 0

for (i in 2:dim(dat)[2]){
	crp = meta[meta$rast_file == names(dat)[i],][2] 
	if (agfoo[crp[1,1]] == 0){
		agfoo[crp[1,1]] <- su[i]
	}else{
		temp <- as.numeric(agfoo[crp[1,1]])
		temp <- temp + su[i]
		agfoo[crp[1,1]]  <- temp
	}
}
temp <- as.numeric(agfoo["LM"])
temp <- temp/sum(meta$SoilFrac == "LM")
agfoo["LM"] <- temp

temp <- as.numeric(agfoo["MM"])
temp <- temp/sum(meta$SoilFrac == "MM")
agfoo["MM"] <- temp

temp <- as.numeric(agfoo["SM"])
temp <- temp/sum(meta$SoilFrac == "SM")
agfoo["SM"] <- temp

temp <- as.numeric(agfoo["Micro"])
temp <- temp/sum(meta$SoilFrac == "Micro")
agfoo["Micro"] <- temp

temp <- as.numeric(agfoo["WS"])
temp <- temp/sum(meta$SoilFrac == "WS")
agfoo["WS"] <- temp

gene <- c(gename,gename,gename,gename,gename)
va <- c("LM","MM","SM","Micro","WS")
value <-c(as.numeric(agfoo["LM"]),as.numeric(agfoo["MM"]),as.numeric(agfoo["SM"]),as.numeric(agfoo["Micro"]),as.numeric(agfoo["WS"]))
agfinal = data.frame(gene,va,value)
final = rbind(final,agfinal)

return(final)
}

#####################
# Here start main code  #
#####################

#read file
chia <- read.table("~/Box Sync/2016/5May/chia/chiA.count.txt",header=T)
nag <- read.table("~/Box Sync/2016/5May/chia/nag.count.txt",header=T)
meta <- read.csv("~/Box Sync/2016/5May/chia/SampleLink4.csv")
meta$rast_file <- sub('_001','',meta$rast_file)
to_remove <- subset(meta, meta$Note == "Poor Quality")
meta<- subset(meta, meta$Note != "Poor Quality")
#Change sample name
names(chia) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('.sam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('Hofmockel','H',names(chia))
names(chia) <- sub('Hofmocke','H',names(chia))
names(nag) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(nag))
names(nag) <- sub('.sam.sorted.bam.idxstats.txt','',names(nag))
names(nag) <- sub('Hofmockel','H',names(nag))
names(nag) <- sub('Hofmocke','H',names(nag))

#remove poor quality
index=names(chia) %in% to_remove$rast_file
badrows = c()
for (i in 1:length(index)){
	if(index[i]){
		badrows = union(badrows,c(i))
	}
}
chia = chia[-badrows]
nag = nag[-badrows]



#initialize final table
gene <- c()
va <- c()
value <- c()
final = data.frame(gene,va,value)

#chia
final <- get_count("ChiA", chia,final,meta)

#nag
final <- get_count("NAG",nag,final,meta)

#plot
library(ggplot2)
pdf("~/Box Sync/2016/5May/chia/cobs_count.pdf",width=6,height=8)
ggplot(final,aes(x=va, y=value, fill=gene))+geom_bar(stat="identity")
dev.off()
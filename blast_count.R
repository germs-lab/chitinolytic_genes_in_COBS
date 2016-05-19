get_count <- function(gename, dat, final,meta){
#make sum
su = c()
for (i in 2:dim(dat)[2]){
	su[i] <- sum(dat[,i])
}

cc <- c()
p <-c()
fp <- c()

for (i in 2:dim(dat)[2]){
	crp = meta[meta$rast_file == names(dat)[i],][1] 
	if (crp == "CC"){
		cc <- cbind(cc,su[i])
	}else if(crp == "P"){
		p <- cbind(p,su[i])
	}else if(crp == "PF"){
		fp <- cbind(fp,su[i])
	}
}


gene <- c(gename,gename,gename)
va <- c("CC","P","FP")
value <- c(mean(cc),mean(p),mean(fp))
stdev <- c(sd(cc),sd(p),sd(fp))
cfinal = data.frame(gene,va,value,stdev)
final = rbind(final,cfinal)

#cal average for each aggregate


lm <- c()
mm <-c()
sm <- c()
micro <-c()
ws <- c()
con <- c()
let <- c()
for (i in 2:dim(dat)[2]){
	crp = meta[meta$rast_file == names(dat)[i],][2] 
	if (crp == "LM"){
		lm <- cbind(lm,su[i])
	}else if(crp == "MM"){
		mm <- cbind(mm,su[i])
	}else if(crp == "SM"){
		sm <- cbind(sm,su[i])
	}else if(crp == "Micro"){
		micro <- cbind(micro,su[i])
	}else if(crp == "WS"){
		ws <- cbind(ws,su[i])
	}
}


gene <- c(gename,gename,gename,gename,gename)
va <- c("LM","MM","SM","Micro","WS")
value <- c(mean(lm),mean(mm),mean(sm),mean(micro),mean(ws))
stdev <- c(sd(lm),sd(mm),sd(sm),sd(micro),sd(ws))
agfinal = data.frame(gene,va,value,stdev)
final = rbind(final,agfinal)

return(final)
}

#####################
# Here start main code  #
#####################

#read file
chia <- read.table("~/Box Sync/2016/5May/chia/summary-count.tsv",header=T)
#nag <- read.table("~/Box Sync/2016/5May/chia/nag.count.txt",header=T)
meta <- read.csv("~/Box Sync/2016/5May/chia/SampleLink4.csv")
meta$rast_file <- sub('_001','',meta$rast_file)
to_remove <- subset(meta, meta$Note == "Poor Quality")
meta<- subset(meta, meta$Note != "Poor Quality")
#Change sample name
names(chia) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('.sam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('Hofmockel','H',names(chia))
names(chia) <- sub('Hofmocke','H',names(chia))
names(chia) <- sub('_R1_001','',names(chia))
names(chia) <- sub('_R2_001','',names(chia))
#names(nag) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(nag))
#names(nag) <- sub('.sam.sorted.bam.idxstats.txt','',names(nag))
#names(nag) <- sub('Hofmockel','H',names(nag))
#names(nag) <- sub('Hofmocke','H',names(nag))

#remove poor quality
index=names(chia) %in% to_remove$rast_file
badrows = c()
for (i in 1:length(index)){
	if(index[i]){
		badrows = union(badrows,c(i))
	}
}
chia = chia[-badrows]
#nag = nag[-badrows]

x = seq(1,dim(chia)[2],by=2)
for (i in x){
	temp <- chia[,x] + chia[,x+1]
}

temp


#initialize final table
gene <- c()
va <- c()
value <- c()
stdev <- c()
final = data.frame(gene,va,value,stdev)

#chia
#final <- get_count("ChiA", chia,final,meta)
final <- get_count("ChiA", temp,final,meta)


#nag
#final <- get_count("NAG",nag,final,meta)

#plot
library(ggplot2)
pdf("~/Box Sync/2016/5May/chia/cobs_blast_count.pdf",width=6,height=8)
ggplot(final,aes(x=va, y=value, fill=gene))+geom_bar(stat="identity")+geom_errorbar(aes(ymax = value + stdev, ymin = value - stdev),width=0.25)
dev.off()
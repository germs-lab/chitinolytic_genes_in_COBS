#heat map
#annotation
#chia <- read.table("~/Box Sync/2016/5May/chia/chiA.count.txt",header=T)
#blast
chia <- read.table("~/Box Sync/2016/5May/chia/summary-count.tsv",header=T)
#add 
x = seq(1,dim(chia)[2],by=2)
for (i in x){
	temp <- chia[,x] + chia[,x+1]
}
chia <- temp
meta <- read.csv("~/Box Sync/2016/5May/chia/SampleLink4.csv")

meta$rast_file <- sub('_001','',meta$rast_file)
to_remove <- subset(meta, meta$Note == "Poor Quality")
meta<- subset(meta, meta$Note != "Poor Quality")

names(chia) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('.sam.sorted.bam.idxstats.txt','',names(chia))
names(chia) <- sub('Hofmockel','H',names(chia))
names(chia) <- sub('Hofmocke','H',names(chia))
names(chia) <- sub('_R1_001','',names(chia))
names(chia) <- sub('_R2_001','',names(chia))

#remove poor quality
index=names(chia) %in% to_remove$rast_file
badrows = c()
for (i in 1:length(index)){
	if(index[i]){
		badrows = union(badrows,c(i))
	}
}
chia = chia[-badrows]


###########################
##                   Crop                   ##
###########################
cc <- c()
p <-c()
fp <- c()
final = data.frame(cc,p,fp)
dat <- chia

for (j in 1:dim(dat)[1]){
	cc <- c()
	p <-c()
	fp <- c()
	for (i in 2:dim(dat)[2]){
		crp = meta[meta$rast_file == names(dat)[i],][1] 
		if (crp == "CC"){
			cc <- cbind(cc,dat[j,i])
		}else if(crp == "P"){
			p <- cbind(p,dat[j,i])
		}else if(crp == "PF"){
			fp <- cbind(fp,dat[j,i])
		}
	}
temp = data.frame(mean(cc),mean(p),mean(fp))
final = rbind(final,temp)
}

machia <- as.matrix(final)
heatmap(machia)

library(reshape2)
melted <- melt(machia)
colnames(melted) <- c("chia_gene","crop","value")
anova(lm(value ~  chia_gene * crop, melted))


###########################
##                Aggregates            ##
###########################
lm <- c()
mm <-c()
sm <- c()
micro <-c()
ws <- c()
final = data.frame(lm,mm,sm,micro,ws)
dat <- chia

for (j in 1:dim(dat)[1]){
	lm <- c()
	mm <-c()
	sm <- c()
	micro <-c()
	ws <- c()
	for (i in 2:dim(dat)[2]){
		crp = meta[meta$rast_file == names(dat)[i],][2] 
		if (crp == "LM"){
			lm <- cbind(lm,dat[j,i])
		}else if(crp == "MM"){
			mm <- cbind(mm,dat[j,i])
		}else if(crp == "SM"){
			sm <- cbind(sm,dat[j,i])
		}else if(crp == "Micro"){
			micro <- cbind(micro,dat[j,i])
		}else if(crp == "WS"){
			ws <- cbind(ws,dat[j,i])
		}
	}
temp = data.frame(mean(lm),mean(mm),mean(sm),mean(micro),mean(ws))
final = rbind(final,temp)
}
machia <- as.matrix(final)

heatmap(machia)
library(reshape2)
melted <- melt(machia)
colnames(melted) <- c("chia_gene","aggregates","value")
anova(lm(value ~  chia_gene * aggregates, melted))



#pcoa
# loading the package used for most multivariate analyses

library(vegan)

sampleREL.dist=vegdist(t(dat),method="bray")
sampleREL.pcoa=cmdscale(sampleREL.dist)
plot(sampleREL.pcoa[,1],sampleREL.pcoa[,2],cex=0,main="standardized sample PCOA")

text(sampleREL.pcoa[,1],sampleREL.pcoa[,2],seq(1,44),cex=1.5)


## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

###########################
##              Main starts here         ##
###########################
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

#change chia gene id into name
map <- read.table("~/Box Sync/Github/germs-lab/chitinolytic_genes_in_COBS/chia.map")
rmap <- as.matrix(map)
temp <- c()
for (i in 1:length(row.names(chia))){
	for (j in 1:dim(rmap)[1]){
	if (row.names(chia)[i] == rmap[j,2]){
		#print(map[j,1])
		#row.names(chia)[i] <- map[j,1]
		temp = cbind(temp,rmap[j,1])
		}
	}	
}
row.names(chia) <- temp

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
row.names(final) <- row.names(chia)
machia <- as.matrix(final)
heatmap(machia)

library(reshape2)
melted <- melt(chia)
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
row.names(final) <- row.names(chia)
machia <- as.matrix(final)

heatmap(machia)
library(reshape2)
melted <- melt(machia)
colnames(melted) <- c("chia_gene","aggregates","value")
anova(lm(value ~  chia_gene * aggregates, melted))



###########################
##      Aggregates & Crop          ##
###########################

clm <- c()
#cmm <-c()
#csm <- c()
cmicro <-c()
cws <- c()
plm <- c()
#pmm <-c()
#psm <- c()
pmicro <-c()
pws <- c()
fplm <- c()
fpmm <-c()
fpsm <- c()
fpmicro <-c()
fpws <- c()

#final = data.frame(clm,cmm,csm,cmicro,cws,plm,pmm,psm,pmicro,pws,fplm,fpmm,fpsm,fpmicro,fpws)
final = data.frame(clm,cmicro,cws,plm,pmicro,pws,fplm,fpmm,fpsm,fpmicro,fpws)
dat <- chia

for (j in 1:dim(dat)[1]){
	clm <- c()
#cmm <-c()
#csm <- c()
cmicro <-c()
cws <- c()
plm <- c()
#pmm <-c()
#psm <- c()
pmicro <-c()
pws <- c()
fplm <- c()
fpmm <-c()
fpsm <- c()
fpmicro <-c()
fpws <- c()

	for (i in 2:dim(dat)[2]){
		crp = meta[meta$rast_file == names(dat)[i],][1] 
		agg = meta[meta$rast_file == names(dat)[i],][2] 
		if (agg == "LM"){
			if (crp == "CC"){
				clm <- cbind(clm,dat[j,i])
			}else if(crp == "P"){
				plm <- cbind(plm,dat[j,i])
			}else if (crp == "PF"){
				fplm <- cbind(fplm,dat[j,i])
			}
			
		}else if(agg == "MM"){
			if (crp == "CC"){
				#cmm <- cbind(cmm,dat[j,i])
			}else if(crp == "P"){
				#pmm <- cbind(pmm,dat[j,i])
			}else if (crp == "PF"){
				fpmm <- cbind(fpmm,dat[j,i])
			}
			
		}else if(agg == "SM"){
			if (crp == "CC"){
				#csm <- cbind(csm,dat[j,i])
			}else if(crp == "P"){
				#psm <- cbind(psm,dat[j,i])
			}else if (crp == "PF"){
				fpsm <- cbind(fpsm,dat[j,i])
			}
			
		}else if(agg == "Micro"){
			if (crp == "CC"){
				cmicro <- cbind(cmicro,dat[j,i])
			}else if(crp == "P"){
				pmicro <- cbind(pmicro,dat[j,i])
			}else if (crp == "PF"){
				fpmicro <- cbind(fpmicro,dat[j,i])
			}
			
		}else if(agg == "WS"){
			if (crp == "CC"){
				cws <- cbind(cws,dat[j,i])
			}else if(crp == "P"){
				pws <- cbind(pws,dat[j,i])
			}else if (crp == "PF"){
				fpws <- cbind(fpws,dat[j,i])
			}
			
		}
	}
#temp = data.frame(mean(clm),mean(cmm),mean(csm),mean(cmicro),mean(cws),mean(plm),mean(pmm),mean(psm),mean(pmicro),mean(pws),mean(fplm),mean(fpmm),mean(fpsm),mean(fpmicro),mean(fpws))

temp = data.frame(mean(clm),mean(cmicro),mean(cws),mean(plm),mean(pmicro),mean(pws),mean(fplm),mean(fpmm),mean(fpsm),mean(fpmicro),mean(fpws))
final = rbind(final,temp)
}

row.names(final) <- row.names(chia)
machia <- as.matrix(final)
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/heatmap1.pdf",width=17,height=10)
heatmap(machia)
dev.off()

#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")  # annHeatmap or annHeatmap2
#library(Heatplus)
#library(vegan)
#library(RColorBrewer)

#install.packages("gplots")
library(gplots) 
heatmap.2(machia)



library(reshape2)
melted <- melt(machia)
colnames(melted) <- c("chia_gene","aggregates","value")
anova(lm(value ~  chia_gene * aggregates, melted))

sum(machia[,1])
chiasum=c()
for (i in 1:dim(machia)[2]){
	chiasum = cbind(chiasum,sum(machia[,i]))
}
row.names(machia)
chiasum[,]
forgraph = data.frame(colnames(machia),chiasum[,])
library(ggplot2)
colnames(forgraph) <- c("crop","value")
ggplot(forgraph,aes(x=crop, y=value))+geom_bar(stat="identity")




###########################
##            ANOVA               ##
###########################
library(reshape2)
library(ggplot2)
dat <- chia
melted <- melt(dat)
ccrp <- c()
cagg <- c()
for (i in 1:dim(melted)[1]){
	crp <- meta[meta$rast_file == melted[i,1],][1] 
	agg = meta[meta$rast_file == melted[i,1],][2] 
	ccrp = cbind(ccrp,toString(crp$Crop))
	cagg = cbind(cagg,toString(agg$SoilFrac))
}
temp = data.frame(t(ccrp),t(cagg))
melted <- cbind(melted,temp)
colnames(melted) <- c("sample","value","crop","agg")
anova(lm(value ~ crop * agg, melted))

##plot for all
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/all_agg.pdf",width=6,height=8)
ggplot(melted,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/all_crop.pdf",width=6,height=8)
ggplot(melted,aes(x=crop,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(melted, measurevar = "value", groupvars=c("agg","sample"))
summarySE(melted, measurevar = "value", groupvars=c("agg"))

#***
subsetPF <- subset(melted, crop == "PF")
#model <- lm(value ~ agg, subsetPF)
#summary(model)
anova(lm(value ~ agg, subsetPF))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/PF_agg.pdf",width=6,height=8)
ggplot(subsetPF,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()

#summarySE(subsetPF, measurevar = "value", groupvars=c("agg"))
summarySE(subsetPF, measurevar = "value", groupvars=c("agg","sample"))
#test <-  subset(subsetPF, agg == "WS")
#mean(test$value)
#sd(test$value)


#***
subsetCC <- subset(melted, crop == "CC")
anova(lm(value ~ agg, subsetCC))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/CC_agg.pdf",width=6,height=8)
ggplot(subsetCC,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(subsetCC, measurevar = "value", groupvars=c("agg"))
summarySE(subsetCC, measurevar = "value", groupvars=c("agg","sample"))

subsetP <- subset(melted, crop == "P")
anova(lm(value ~ agg, subsetP))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/P_agg.pdf",width=6,height=8)
ggplot(subsetP,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(subsetP, measurevar = "value", groupvars=c("agg","sample"))

subsetLM <- subset(melted, agg == "LM")
anova(lm(value ~ crop,subsetLM))

#**
subsetmicro <- subset(melted, agg == "Micro")
anova(lm(value ~ crop,subsetmicro))

subsetWS <- subset(melted, agg == "WS")
anova(lm(value ~ crop,subsetWS))
#ggplot(melted,aes(x=agg,y=value))+geom_bar(stat="identity")
#pcoa
# loading the package used for most multivariate analyses


###########################
##      ANOVA   - annotation      ##
###########################
anochia <- read.table("~/Box Sync/2016/5May/chia/chiA.count.txt",header=T)

meta <- read.csv("~/Box Sync/2016/5May/chia/SampleLink4.csv")

meta$rast_file <- sub('_001','',meta$rast_file)
to_remove <- subset(meta, meta$Note == "Poor Quality")
meta<- subset(meta, meta$Note != "Poor Quality")

names(anochia) <- sub('.sam.bam.sorted.bam.idxstats.txt','',names(anochia))
names(anochia) <- sub('.sam.sorted.bam.idxstats.txt','',names(anochia))
names(anochia) <- sub('Hofmockel','H',names(anochia))
names(anochia) <- sub('Hofmocke','H',names(anochia))
names(anochia) <- sub('_R1_001','',names(anochia))
names(anochia) <- sub('_R2_001','',names(anochia))

#remove poor quality
index=names(anochia) %in% to_remove$rast_file
badrows = c()
for (i in 1:length(index)){
	if(index[i]){
		badrows = union(badrows,c(i))
	}
}
anochia = anochia[-badrows]

dat <- anochia
melted <- melt(dat)
ccrp <- c()
cagg <- c()
for (i in 1:dim(melted)[1]){
	crp <- meta[meta$rast_file == melted[i,2],][1] 
	agg = meta[meta$rast_file == melted[i,2],][2] 
	ccrp = cbind(ccrp,toString(crp$Crop))
	cagg = cbind(cagg,toString(agg$SoilFrac))
}
temp = data.frame(t(ccrp),t(cagg))
melted <- cbind(melted,temp)
colnames(melted) <- c("genome","sample","value","crop","agg")
anova(lm(value ~ crop * agg, melted))

pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/anno_all_agg.pdf",width=6,height=8)
ggplot(melted,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/anno_all_crop.pdf",width=6,height=8)
ggplot(melted,aes(x=crop,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(melted, measurevar = "value", groupvars=c("agg","sample"))
summarySE(melted, measurevar = "value", groupvars=c("agg"))

#***
subsetPF <- subset(melted, crop == "PF")
#model <- lm(value ~ agg, subsetPF)
#summary(model)
anova(lm(value ~ agg, subsetPF))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/anno_PF_agg.pdf",width=6,height=8)
ggplot(subsetPF,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
#summarySE(subsetPF, measurevar = "value", groupvars=c("agg"))
summarySE(subsetPF, measurevar = "value", groupvars=c("agg","sample"))
#test <-  subset(subsetPF, agg == "WS")
#mean(test$value)
#sd(test$value)

#***
subsetCC <- subset(melted, crop == "CC")
anova(lm(value ~ agg, subsetCC))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/anno_CC_agg.pdf",width=6,height=8)
ggplot(subsetCC,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(subsetCC, measurevar = "value", groupvars=c("agg"))
summarySE(subsetCC, measurevar = "value", groupvars=c("agg","sample"))

subsetP <- subset(melted, crop == "P")
anova(lm(value ~ agg, subsetP))
pdf("/Users/jinchoi/Box\ Sync/Github/germs-lab/chitinolytic_genes_in_COBS/anno_P_agg.pdf",width=6,height=8)
ggplot(subsetP,aes(x=agg,y=value))+stat_summary(fun.y="mean",geom="bar") + stat_summary(fun.data = mean_se , geom = "errorbar", width = .1)
dev.off()
summarySE(subsetP, measurevar = "value", groupvars=c("agg","sample"))




library(vegan)

sampleREL.dist=vegdist(t(dat),method="bray")
sampleREL.pcoa=cmdscale(sampleREL.dist)
plot(sampleREL.pcoa[,1],sampleREL.pcoa[,2],cex=0,main="standardized sample PCOA")

text(sampleREL.pcoa[,1],sampleREL.pcoa[,2],seq(1,44),cex=1.5)


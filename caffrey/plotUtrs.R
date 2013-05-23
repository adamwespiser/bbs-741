library(ggplot2)
library(plyr)
library(reshape2)


df <- read.csv("comb.csv",stringsAsFactors=FALSE)
df <- transform(df,cdLength=transLength - threeUtrLength - fiveUtrLength)
melt.df <- melt(df,id.vars=c("id","species"),measure.vars=c("threeUtrLength","fiveUtrLength","cdLength","transLength"))
ggplot(melt.df,aes(x=log(value),color=species))+geom_density()+facet_wrap(~variable)+theme_bw() + ggtitle("Density plot of transcript regions for different species")
ggsave("densityPlotBioRegion.png",height=5,width=7)


df$sumTransCount <- 1

ddply(df,.(species),numcolwise(sum)) -> specTrans.df
merge(ddply(df,.(species),numcolwise(mean))[c("species","transLength","threeUtrLength","fiveUtrLength","cdLength")],specTrans.df[c("species","sumTransCount")],by="species") -> specAverages.df

meltSpecAve.df <- melt(specAverages.df,id.var="species")

trans.df <- meltSpecAve.df[meltSpecAve.df$variable != "sumTransCount" & meltSpecAve.df$variable != "transLength",]
trans.df$variable <- ordered(trans.df$variable ,levels = c("threeUtrLength","cdLength","fiveUtrLength"))


ggplot(trans.df,aes(x=species,y=value,fill=variable,order=factor(variable,levels=c("threeUtrLength","cdLength","fiveUtrLength"))))+geom_histogram(stat="identity",position="stack")+theme_bw()+ylab("avearge transcript feature length(bp)")+ggtitle("average length of transcript regions over species")
ggsave("transRegionLength.png",height=7,width=7)





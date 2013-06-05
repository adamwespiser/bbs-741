# Adam Wespiser
# GPL v2.0



library(Rmixmod)
library("caTools")
library(DNAcopy)
genomdat <- read.table("gbm.chr7and9.01.dat",header=T, row.names=1, check.names = FALSE, sep="\t")
data = genomdat[,3]

#load data
for (i in 2:9){
  fname=paste("gbm.chr7and9.0",i,".dat",sep="") # filename
  genomdat.new = read.table(fname,header=T, row.names=1, check.names = FALSE, sep="\t")
  data = cbind(data,genomdat.new[,3])  # col. for new data.
}
for (i in 10:15){fname=paste("gbm.chr7and9.",i,".dat",sep="") # filename
                 genomdat.new = read.table(fname,header=T, row.names=1, check.names = FALSE, sep="\t")
                 data = cbind(data,genomdat.new[,3])  # col. for new data.
}


CNA.object <- CNA(as.matrix(data), genomdat[,1], as.numeric(genomdat[,2])) 
smoothed.CNA.object <- smooth.CNA(CNA.object)


# Vary over alpha values


alphaRun = 0 # set to one to run...
if (alphaRun == 1){  # lets not run this is we don't have to...
#segment with different values of the parameter alpha
alpha.vals = c(0.01,0.02,0.03,0.4,0.5) # List parameter choices for alpha
for (i in 1:5){
  segment.smoothed.CNA.object=segment(smoothed.CNA.object, min.width=5, verbose=2,alpha = alpha.vals[i]) # set each value for alpha in turn
  # Different output file for each value of alpha:
  output.file = paste("gbm.alpha",alpha.vals[i],".seg.txt",sep="")
  write.table(segment.smoothed.CNA.object$output, output.file,sep="\t",row.names=F,col.names=F, quote=F)
}
# segment.smoother.CNA.object is the segmented object for further analysis!
}


# vary over min.width params 


widthRun= 0
if (widthRun == 1){  # lets not run this is we don't have to...
#segment with different values of the parameter minimum width parameters
min.width.vals = c(2,5) # List parameter choices for min.width 
for (i in seq_along(min.width.vals)){
  segment.smoothed.CNA.object=segment(smoothed.CNA.object, min.width=min.width.vals[i], verbose=2,alpha = 0.01) # set each value for alpha in turn
  # Different output file for each value of min.value:
  output.file = paste("gbm.min.width=",min.width.vals[i],".seg.txt",sep="")
  write.table(segment.smoothed.CNA.object$output, output.file,sep="\t",row.names=F,col.names=F, quote=F)
}
}#end of width run



### Run the moving averages...

genomdat <- read.table("gbm.chr7and9.01.dat",header=T, row.names=1, check.names = FALSE, sep="\t") # cols are : "Chromosome", "PhysicalPosition", "COTES_p_TCGAaffxB8_9a_S_GenomeWideSNP_6_A03_292900"
data = genomdat[,3]
# Select k and thresholds
window.width = 100

moving.ave= runmean(x=data,k=window.width)
segmented.data = moving.ave

# Compute Tukey's five statistics using R's "fivenum" function

tukeyFive <- fivenum(segmented.data)
threshold1 <- tukeyFive[2]
threshold2 <- tukeyFive[4]

#segment data w/ Tukey's lower and upper shoulder

segmented.data[moving.ave < threshold1] = mean(data[moving.ave<threshold1])
segmented.data[moving.ave >= threshold2] = mean(data[moving.ave>threshold2])
midlevel = (moving.ave >= threshold1) & (moving.ave<threshold2)
segmented.data[midlevel] = mean(data[midlevel])

#library(ggplot2)
# assign the data to data.frames for ggplot to accept
genomeData.df= data.frame(signal=data,genomeCoord=seq_along(moving.ave))
segmentedData.df = data.frame(signal=segmented.data,genomeCoord=seq_along(segmented.data))
#Fig 6
ggplot(genomeData.df,aes(x=genomeCoord,y=signal))+geom_point(color="grey",alpha=I(0.02))+geom_line(data=segmentedData.df,aes(x=genomeCoord,y=signal,color="red"))+ggtitle("Attempt setting thresholds using \nTukey's 5 number statistic(2,4 numbers taken)")+theme_bw()

#okay, that was not a very good idea, let's try what we learned in class...

# visually determine cutoff (more objectively, we could train a gausian mixture model..)

# Fig 4
ggplot(data.frame(x=moving.ave),aes(x))+geom_density(adjust=1/10,fill="grey")+geom_vline(xintercept=0.4,color="green")+geom_vline(xintercept=-0.2,color="red")+theme_bw()+xlab("signal")+ylab("density")+ggtitle("PDFof genome signal\nwindow=100,mean\nred = lower.threshold\ngreen = upper.threshold")+xlim(-1,1)

#set thresholds:
threshold1 <- -0.2
threshold2 <- 0.4
segmented.data = moving.ave
segmented.data[moving.ave < threshold1] = mean(data[moving.ave<threshold1])
segmented.data[moving.ave >= threshold2] = mean(data[moving.ave>threshold2])
midlevel = (moving.ave >= threshold1) & (moving.ave<threshold2)
segmented.data[midlevel] = mean(data[midlevel])

#assign data to df's for ggplot
genomeData.df= data.frame(signal=data,genomeCoord=seq_along(moving.ave))
segmentedData.df = data.frame(signal=segmented.data,genomeCoord=seq_along(segmented.data))

#Fig 5
ggplot(genomeData.df,aes(x=genomeCoord,y=signal))+geom_point(alpha=I(0.02),color="grey")+geom_line(data=segmentedData.df,aes(x=genomeCoord,y=signal,color="red"))+theme_bw()+ggtitle("running mean w/ k=100,thresholds= -0.2,0.4\nGrey points = signal from genome\nRed Line = segmentation estimate")


# create a data frame with multiple runmean trials
segTest.df <- data.frame(geneCoord = rep(seq_along(data),6),
                   k = rep(c(10,25,50,100,200,500),each=length(data)),
		mean = c(runmean(x=data,10),
		         runmean(x=data,25),
		         runmean(x=data,50),
		         runmean(x=data,100),
		         runmean(x=data,200),
			 runmean(x=data,500)))
# plot the runmean trials and seperate by k(windowsize)
# Fig 7
ggplot(segTest.df,aes(x=geneCoord,y=mean))+geom_line()+facet_wrap(~k,ncol=2)+ggtitle("effect of different window sizes\non running average calculation")

# create a data frame with multiple runmed trials
# In some signals processing median is more robust to changes and less bias.
# ...I'm not sure that's the case here...
segTestmed.df <- data.frame(geneCoord = rep(seq_along(data),6),
                   k = rep(c(11,25,51,101,201,501),each=length(data)),
		mean = c(runmed(x=data,11),
		         runmed(x=data,25),
		         runmed(x=data,51),
		         runmed(x=data,101),
		         runmed(x=data,201),
			 runmed(x=data,501)))

# Fig 8, take a look at median runs
ggplot(segTestmed.df,aes(x=geneCoord,y=mean))+geom_line()+facet_wrap(~k,ncol=2)+ggtitle("effect of different window sizes\non running median calculation")


#Let's try a guassian mixture model to assign threshold parameters in a more objective way...
#library(Rmixmod)

#get the mean'd data from the previous df.
averaged.data = segTest.df[segTest.df$k==100,"mean"]

#create a mixed model with 2 to 6 clusters from the averaged data
dataMixModel=mixmodCluster(averaged.data,nbCluster=2:6)

#extract the parameters of the model, (only mean is used)
# TODO: merge means by variance if ( pr(x|model1)/pr(x|model2) < discrimRatio) all x.
variances <- unlist(dataMixModel@results[[1]]@parameters@variance)
means <-  as.vector(dataMixModel@results[[1]]@parameters@mean)

#Fig 10, visualize the data and guassian models
hist(dataMixModel)
#function to assign data points to neartest mean

# assign groups to the guassian it's most likely to belong to.
assignGroupZ <- function(x,vec.means=means,vec.sd=sd(variances)){
 z = abs(x-vec.means)/vec.sd
 means[which(z == min(z))]
}
# create a function to assign averaged data to the nearest mixture model
assignGroup <- function(x,vmeans=means,vec.sd=sd(variances)){
 diff = (means -x)^2
 means[which(diff == min(diff))]
}

gaussSegment.df <- data.frame(
segGroup = sapply(averaged.data,function(x)assignGroupZ(x)),
windowAve = averaged.data,
genomeSignal = data,
genomeCoord = seq_along(data))

#Fig 9
ggplot(gaussSegment.df, aes(x=genomeCoord,y=genomeSignal))+geom_point(alpha=I(0.02),color="grey")+geom_line(aes(x=genomeCoord,y=segGroup,color="red"))+theme_bw()+ggtitle("Mixture Model assignment of segmentation groups")

CNA.object <- CNA(as.matrix(data), genomdat[,1], as.numeric(genomdat[,2])) # The samples will automatically be named


#Fig 10
CNA.object <- CNA(data) # from geom_data[,3]
smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object=segment(smoothed.CNA.object, min.width=5, verbose=2,alpha = 0.01) # set each value for alpha in turn
  # Different output file for each value of alpha:
  output.file = paste("gbm.segTest",1,".seg.txt",sep="")
  write.table(segment.smoothed.CNA.object$output, output.file,sep="\t",row.names=F,col.names=F, quote=F)


# Unused plots I might need later...
#plots:
#plot(segment.smoothed.CNA.object,plot.type="w") #one plot per sample
#plot(segment.smoothed.CNA.object, plot.type="s",sbyc.layout=c(5,5)) #by chr
#plotSample(segment.smoothed.CNA.object,sampleid=1,chromlist=7,xmaploc=TRUE)
#zoomIntoRegion(segment.smoothed.CNA.object,sampleid=1,chrom=7,
#maploc.start=53000000,maploc.end=58000000)#by region


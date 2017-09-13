
library('scatterplot3d')
library('foreach')
library('doMC') 
library('fields')
library('MASS')

##############################
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"
#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
#SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation"/

#SimPath<-"/media/zoe/PUBLIC/DataCopy/"


#ScriptLocation<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/"
ScriptLocation<-"/home/csdori/ksCSD_2014/trunk/"
SumPlotPlace<-paste0(ScriptLocation,"plots/summary/")
source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
dir.create(paste0(ScriptLocation,"plots/"))
dir.create(SumPlotPlace)

#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"

locationsKernel<-numeric()
locationsData<-numeric()

CellType<-"Y"
SNR<-0
#lambda<-0



if(CellType=="Y"){
  #SimPath<-paste0(ScriptLocation,"simulation/")
  #SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/Sim_Yshaped/"
  SimPath<-paste0(ScriptLocation,"Sim_Yshaped/")
  
  #SimPath<-paste0("/home/zoe/Dropbox/skCSD/data/")
  
  Versions<-c(0:5)
  Yelectrodes<-rep(c(2,4,8,16),each=length(Versions))
  #Yelectrodes<-16 #rep(c(16),each=length(Versions))
  ElNumb<-4*Yelectrodes
  #locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Rot_symm_d50_ver",Versions) 
  locationsKernel<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
  #locationsKernel<-"Y_el4x16_Elrand_symm_d50_ver0_kiserlet" #"Mainen3D_ElGrid2"
  #locationsData<-paste0("Y_el4x",Yelectrodes,"_Elrand_symm_d50_ver",Versions)
  locationsData<-locationsKernel
  outname1<-"/kernelOut_4"
}


####Cross validation function
#ksCSD cross-validation
#cross.valid function gives back the value of the crossvalidation error and also the crossvalidated potential
#source("/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/alprogik/errorCalcFun.R")
source(paste0(ScriptLocation,"alprogik/Colors_BlueRed.R"))
source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))






library('scatterplot3d')
library('foreach')
library('doMC') 
#library('fields')
library('MASS')

##############################
#SimPath<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/"
#SimPath<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/"

#SumPlotPlace<-"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/plots/summary/"
#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"


#legend(max(ElNumb)+4,mean(ErrorFixElectrodes),c("Grid ",paste("Rand",1:5)),pch=20,col=1:6,title="El Geom")





#Imageplot at different time segments to compare the results
skCSDTime<-numeric()

for(Dlength in 1:length(locationsKernel)){
  # Dlength<-2
  inname<-paste0(SimPath, locationsKernel[Dlength],"/")
  outname<-paste0(SimPath, locationsData[Dlength])
  
  SimTime<-as.matrix(read.table(paste0(outname,'/time')))
  TimeInstant<-which(SimTime==5) #45
  TimeInstant<-TimeInstant+2
  #BestEMatr<-BestEMatrFix
  
  if(file.exists(paste0(outname,outname1,"/skCSDall_M",512,"_R",BestEMatr[1,Dlength], "lambda",BestEMatr[2,Dlength]))==FALSE) next
  skCSD.all<-array(0,c(86, 561))
  skCSD.all.part<-as.matrix(read.table(paste0(outname,outname1,"/skCSDall_M",512,"_R",BestEMatr[1,Dlength], "lambda",BestEMatr[2,Dlength])))
  #Reading in tge currents
  cat(Dlength)
  
  
  SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
  seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
  
  for(i in 1: seg.nb){
    sameplace<-numeric() 
    sameplace<-as.integer(unlist(strsplit(SameplaceAll," ")[i]))
    #Plot the current related ti the different segments
    #  matplot(cbind(t(C.calc[c(sameplace),])),t='l', xlab="Time steps", ylab="ksCSD", main="Values of current at the same location")
    if(length(sameplace)==1) skCSD.all[i,]<- skCSD.all.part[sameplace,]
    if(length(sameplace)>1) skCSD.all[i,]<-as.matrix(colSums( skCSD.all.part[c(sameplace),]),nrow=1)
  }
  
  
  skCSDTime<-c(skCSDTime,c(skCSD.all[,TimeInstant])) #/max(c(skCSD[,TimeInstant]
  coordsEnd<-as.matrix(read.table(paste0(outname,"/coordsend_x_y_z")))
  
}

seglength<-as.matrix(read.table(paste0(outname,'/seglength')))
funaramvonal<-function(x) x/seglength

membCurr<-as.matrix(read.table(paste0(outname,'/membcurr')))
membCurr<-apply(membCurr,2,funaramvonal)[,TimeInstant]
skCSDTimeAll<-matrix(skCSDTime,ncol=length(locationsKernel))


png(paste0(SumPlotPlace,"Y_NewCool_time5.png"), width=800, height=800, pointsize=32)
# clear objects  
#graphics.off()    
par(oma=c(2,2,2,2))  
#par(mfrow=c(1,2))
layout(rbind(c(1,2,3,4,5)),widths=c(0.5,0.2,1.5,3,0.2),respect=FALSE) 

#par(mfrow=c(1,1))
col1<-ColoursDori(membCurr)[[1]]
ExtVal<-ColoursDori(membCurr)[[2]]
par(mar=c(0.5,1,1,0.5)) 
image(0:1,1:86,t(as.matrix(membCurr)),col=col1, zlim=c(-ExtVal,ExtVal),main="CSD",xaxt="n",frame.plot=TRUE)
abline(h=0.5)
abline(v=1)
lines(c(membCurr/ExtVal/2)+0.5,1:86,col="BLACK")
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(as.matrix(membCurr)),col=col1,axis.pos=4,zlim=c(-ExtVal,ExtVal))
#layout.show(2)
col2<-ColoursDori(skCSDTimeAll)[[1]]
ExtVal<-ColoursDori(skCSDTimeAll)[[2]]
par(mar=c(0.5,3,1,0.5)) 
image(1:4,1:dim(skCSDTimeAll)[1],t(skCSDTimeAll[,seq(1,24,by=6)]),col=col2, yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD",frame.plot=TRUE)
abline(h=0.5)
abline(v=4.5)
axis(1, at=c(1:4),labels=c(8,16,32,64),las=2)
for(i in 1:4){
  lines(t(skCSDTimeAll[,seq(1,24,by=6)])[i,]/(ExtVal*2)+i,1:86)
}

col2<-ColoursDori(skCSDTimeAll)[[1]]
par(mar=c(0.5,0.5,1,1.5)) 
image(1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[1],1:dim(t(skCSDTimeAll[,-seq(1,24,by=6)]))[2],t(skCSDTimeAll[,-seq(1,24,by=6)]),col=col2,yaxt='n',zlim=c(-ExtVal,ExtVal),xaxt="n",main="ksCSD - Random ",frame.plot=TRUE)
abline(h=0.5)
axis(1, at=c(1:4)*5-2,labels=c(8,16,32,64),las=2)
#,legend.shrink=0.3,legend.width=0.8,legend.mar=3)
abline(v=c(seq(0,dim(skCSDTimeAll)[2],by=5)+0.5),col="BLACK", lty=4)

for(i in 1:4){
  lines(t(rowMeans(skCSDTimeAll[,((i-1)*6+1):(i*6)]))/(ExtVal*2)+5*(i-1)+3,1:86,col="BLACK")
}
par(mar=c(0.5,0.1,1,0.5)) 
image.scale(t(skCSDTimeAll),col=col2,axis.pos=4,zlim=c(-ExtVal,ExtVal))


dev.off()








library('scatterplot3d')
library('foreach')
library('doMC')
#require(gWidgets)
#library('fields')
library('MASS')
#options("guiToolkit"="RGtk2")
registerDoMC(cores=24)

#Directory with all the data... output of LFPy

where2save<-DirData #"/home/dcserpan/Documents/skCSD/skCSDnew/trunk/simulation/Y_el16_m200_600_d100"
outnameSub<-"kernelOut_4/"
wherearewenow<-getwd()
dir.create(where2save)

#The method calculates the solutions in case of different parameters (width and number of Gaussian functions) it is possible to set the min$
basis.width.min<-3
basis.width.max<-7
basis.width.step<-1

basis.number.min<--5#4
basis.number.max<--1#1#9
basis.number.step<-1
#M=2^k.... This is how the basis number is calculated k goes from basis.number.mn till max

kCSDversion<- "wobg" #"wbg" #"kCSDversion<- "wbg"" wobg: without background

setwd(DirData)

#reading in electrode coordinates

elec.kord <- as.matrix(read.table("elcoord_x_y_z"))
elec.kord<-matrix(elec.kord,ncol=3)

outname<-paste(where2save,"/",outnameSub,sep="")
#if(file.exists(outname)) file.remove(list.files("outname"))
dir.create(outnameSub)

#number of electrodes
#LFP
#lfp<-read.table('myLFP',colClasses='numeric' )
LFP<-as.matrix(read.table("myLFP" ))
lfpmatr<-dim(LFP)
LFP<-matrix(as.numeric(LFP),nrow=lfpmatr[1])

#Don't forget to write out also these parameters
#m(lfpmatr)
#morphology
morpho <-read.table("segcoordinates.txt")
#connection information
connections <- read.table("connections.txt")
#cat(connections)
#membrane currents
memb.currents <- as.matrix(read.table("membcurr"))
cmatr<-dim(memb.currents )
memb.currents <-matrix(as.numeric(memb.currents),nrow=cmatr[1])
#rm(lfpmatr)
#length of segments
seg.length <- as.numeric(as.matrix(read.table("seglength")))

rangeX<-range(elec.kord[,1],morpho[,1])
rangeY<-range(elec.kord[,2],morpho[,2])
rangeZ<-range(elec.kord[,3],morpho[,3])
#################################
setwd(wherearewenow)

if( kCSDversion== "wobg") {
  source('utils/kernel_basisfunctions_Regularizal.R',local=TRUE,verbose=TRUE)
  sCSD_currents<-ksCSD_all( basis.width.min, basis.width.max, basis.width.step, basis.number.min, 
                            basis.number.max, basis.number.step, LFP , elec.kord ,memb.currents,seg.length  ,where2save) #,R.V ,source.V.t)
}
if( kCSDversion== "wbg"){
  source('utils/kernel_basisfunctions_new_with_ksCSD.R',local=TRUE,verbose=TRUE)
  
  R.V<-300 #basis.width.max #width of background basis funcitons:
  loc.X<-seq(min(rangeX),max(rangeX),by=R.V)
  loc.Y<-seq(min(rangeY),max(rangeY),by=R.V)
  loc.Z<-seq(min(rangeZ),max(rangeZ),by=R.V)
  
  source.V.t<- as.matrix(expand.grid(loc.X, loc.Y, loc.Z))
  
  cat(paste0("Number of spherical basis functions:",dim(source.V.t)[1]))
  #Rprof('profilekernel_sim.out',memory.profiling=TRUE)
  ksCSD_currents<-ksCSD_all( basis.width.min, basis.width.max, basis.width.step, basis.number.min, 
                             basis.number.max, basis.number.step, LFP , elec.kord ,memb.currents,seg.length  ,where2save,R.V ,source.V.t)
  
}

#Rprof(NULL)
#image(t(ksCSD_currents))

setwd(wherearewenow)

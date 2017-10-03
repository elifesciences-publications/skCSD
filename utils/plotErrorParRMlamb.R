library(data.table)
library(MASS)
#ScriptLocation<-"/media/zoe/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/"
ScriptLocation<-"/home/zoe/agy/ksCSD_2014/trunk/"
SumPlotPlace<-paste0(ScriptLocation,"plots/summary/")
source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
dir.create(paste0(ScriptLocation,"plots/"))
dir.create(SumPlotPlace)

#outname<-"/media/BA0ED4600ED416EB/agy/ksCSD_2014/trunk/simulation/Y_dense/"

locationsKernel<-numeric()
locationsData<-numeric()

CellType<-"Domi" #"gangNew"#"GangBN"#"YVer0"#"Y"# "gangNew" #Yreg"#"BS" #"Yreg" #"BSOLD"
SNR<-0
lambda<-0


SimPath<-paste0("/home/zoe/agy/skCSD_public/simulation/") #cell_Y_symm/

#SimPath<-paste0("/home/zoe/Dropbox/skCSD/data/")

ElNumb<-32

locationsKernel<-"cell_Y_symm/"

locationsData<-locationsKernel
outname1<-paste0("skCSDreconstruct/")


source(paste0(ScriptLocation,"alprogik/Colors_BlueRed.R"))
source(paste0(ScriptLocation,"alprogik/errorCalcFun.R"))
#source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
source(paste0(ScriptLocation,"alprogik/JustErrorCalcFun.R"))




lambda.all<-c("1e-05","1e-04","0.001","0.01","0.1" ) #c(0.1,0.01,0.001)#
R.all<-2^(3:7) #c(8,16,128)#


sigma<-0.3
JustL1Reg<-function(SimPath,locationsKernel,locationsData,SNR,M){
  Reg<-TRUE
  ErrorResults <- vector("list", length(locationsKernel))
  
  for(Dlength in 1:length(locationsKernel)){
    
    El1ErrorRaw<-array(0,c(length(R.all),length(lambda.all)))
    
    inname<-paste0(SimPath, locationsKernel[Dlength],"/")
    outname<-paste0(SimPath, locationsData[Dlength],"/")
    
    if(file.exists(paste0(outname,outname1))==FALSE) dir.create(paste0(outname,outname1))
    
    SameplaceAll<-readLines(paste0(inname,outname1,"/sameplace.txt"))
    seg.nb<-max(as.matrix(read.table(paste(inname,outname1,'/connection_touched_once',sep=''))))
    LFP<-as.matrix(fread(paste(outname,"/myLFP",sep="")))
    
    if (SNR!=0){
      LFPvariance<-var(c(LFP))
      Noisevariance<-LFPvariance/SNR
      NoiseLFP<-rnorm(length(c(LFP)), mean = 0, sd = sqrt(Noisevariance))
      LFP<-LFP+NoiseLFP
    }
    LFPOriginal<-LFP
    
    ElCoords<-as.matrix(read.table(paste(inname,"elcoord_x_y_z",sep="")))
    ElCoords<-matrix(ElCoords,ncol=3)
    Xel<-unique(ElCoords[,2])
    XelNb<-length(Xel)
    Yel<-unique(ElCoords[,3])
    
    
    for(lamb in 1:length(lambda.all)){
      for(r in 1:length(R.all)){
        
        LFP<-LFPOriginal
        lambda<-lambda.all[lamb]
        R<-R.all[r]
        
        source(paste0(ScriptLocation,"alprogik/KernSmoothDistance.R"))
        
        CurrNew<-KernSmoothDistance(30,inname, outname, outname1,"30")
        
        el2ignorename<-paste(inname,outname1,"/El2Ignore_M",M,"_R",R,sep="")
        if(file.exists(el2ignorename)){
          El2Ignore<-c(as.matrix(read.table(el2ignorename)))
          LFP<-LFP[-El2Ignore,] }
        
        
        
        Kmatr<-as.matrix(read.table(paste(inname,outname1,"/K_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        Ktildematr<-as.matrix(read.table(paste(inname,outname1,"/Ktilda_M",M,"_R",R,"lambda",lambda.all[lamb],sep="")))
        
        Tmatr<-(4*pi*sigma)*t(Ktildematr)%*%ginv(Kmatr)
        
        
        
        skCSD.all.part<-Tmatr%*%LFP
        write.table(skCSD.all.part, paste(outname,outname1,"/ksCSD_Matr",M,"_R",R,"lambda",
                                          lambda.all[lamb], "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        
        
        
        skCSD.all<-array(0,c(seg.nb, dim(LFP)[2]))
        
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
        
        
        if(CellType=="BSOLD") {skCSD.all<-skCSD.all[c(2,1,3:52),]}
        El1ErrorRaw[r,lamb]<-L1Error(skCSD.all,CurrNew)[[1]]
        
        
        #      cvErrorOut<-cross.valid(LFP,Tmatr,1:dim(Tmatr)[1])
        #ErrorCV[r,lamb]<-cvErrorOut[[1]]
        cat("Numb to rteplace:", length(R.all)*(lamb-1)+ r)
        #CVConfig[[length(R.all)*(lamb-1)+ r]]<-as.matrix(cvErrorOut[[4]])
        #write.table(as.matrix(cvErrorOut[[4]]),paste(outname,outname1,"/CVContribution_M",M,"_R",R,"lambda",lambda.all[lamb], 
        #                                             "SNR",SNR,sep=""),row.names = FALSE,col.names = FALSE)
        #ErrorContribution
      }} #lambda and R
    
    write.table(El1ErrorRaw, paste(outname,outname1,"/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",M,sep=""),row.names = FALSE,col.names = FALSE)
    
  }#Dlength
  #return(list(ErrorResults, CVConfig))
}


korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,16)
korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,32)
korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,64)

korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,128)

korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,256)
korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,512)

korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,1024) #behalt
#korte<-JustL1Reg(SimPath,locationsKernel,locationsData,SNR,M)

M<-128
outname<-paste0(SimPath, locationsData[1])


hibak16<-as.matrix(read.table(paste(outname,outname1,
                                    "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",16,sep="")))

hibak32<-as.matrix(read.table(paste(outname,outname1,
                                    "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",32,sep="")))

hibak64<-as.matrix(read.table(paste(outname,outname1,
                                     "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",64,sep="")))

hibak128<-as.matrix(read.table(paste(outname,outname1,
        "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",128,sep="")))
hibak256<-as.matrix(read.table(paste(outname,outname1,
                                     "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",256,sep="")))
hibak512<-as.matrix(read.table(paste(outname,outname1,
                                     "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",512,sep="")))
hibak1024<-as.matrix(read.table(paste(outname,outname1,
                                     "/ErrorL1Smoothed_Noise","_SNR",SNR,"_M",1024,sep="")))

zrange<-range(hibak32,hibak128,hibak512,hibak1024)
png("YsimmErrorM.png",un="in", width=5*1.5,  height=2*1.5, res=300  )

par(mfrow=c(1,5),  mar=c(4,4,2,0),oma=c(0,0,0,2))

image(1:5, 1:5,hibak32, xaxt="n",yaxt="n",xlab="R",ylab="Lambda", zlim=zrange,col=tim.colors(200), main="M = 32")
axis(1,1:5, R.all)
axis(2,1:5, lambda.all)


image(1:5, 1:5,hibak128, xaxt="n",yaxt="n",xlab="R",ylab="", zlim=zrange,col=tim.colors(200), main="M = 128")
axis(1,1:5, R.all)
#axis(2,1:5, lambda.all)

image(1:5, 1:5,hibak512, xaxt="n",yaxt="n",xlab="R",ylab="", zlim=zrange,col=tim.colors(200), main="M = 512")
axis(1,1:5, R.all)
#axis(2,1:5, lambda.all)

image(1:5, 1:5,hibak1024, xaxt="n",yaxt="n",xlab="R",ylab="", zlim=zrange,col=tim.colors(200), main="M = 1024")
axis(1,1:5, R.all)
#axis(2,1:5, lambda.all)


plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim=c(0,1),ylim=c(0,1))
image.plot(hibak128,legend.only = TRUE, zlim=zrange)#,  smallplot=c(0.1,0.1,0.3,0.9))
dev.off()

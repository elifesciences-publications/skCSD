library('scatterplot3d')
library('foreach')
library('doMC')
require(gWidgets)
library('fields')
library('MASS')

#options("guiToolkit"="RGtk2")


kernel_calc <- function() {
  #Get the location of the current working directory
  wherearewenow<-getwd()
  #for calculating the kernel functions
  kernel_calculate <- function(h,...) {
    #where to save the outputs
    where2save<-svalue(simulation.location.name)
    dir.create(where2save)
    
    
    outnameSub<-"skCSDreconsruct/"
    
    outname<-paste(where2save,"/",outnameSub,sep="")
    dir.create(outname)
    ##Parameters
    #The method calculates the solutions in case of different parameters (width and number of Gaussian functions) it is possible to set the minimum and maximum value of these parameter and the step sizes
    basis.width.min<-svalue(basis.width.min.g)
    basis.width.max<-svalue(basis.width.max.g)
    basis.width.step<-svalue(basis.width.step.g)
    
    basis.number.min<-svalue(basis.number.min.g)
    basis.number.max<-svalue(basis.number.max.g)
    basis.number.step<-svalue(basis.number.step.g)
    

    
    
    
    #setwd(where2save)
    #reading in electrode coordinates
    
    elec.kord <- as.matrix(read.table(svalue(elec.cord.location)))
    elec.kord<-matrix(elec.kord,ncol=3)
    
    #number of electrodes
    #LFP
    #lfp<-read.table('myLFP',colClasses='numeric' )
    LFP<-as.matrix(read.table(svalue(lfp.location )))
    lfpmatr<-dim(LFP)
    LFP<-matrix(as.numeric(LFP),nrow=lfpmatr[1])
  
    #Don't forget to write out also these parameters



    #m(lfpmatr)
    #morphology
    morpho <-read.table(svalue(membrane.current.location))
    #connection information
    connections <- read.table(svalue(connection.info.location))
    #cat(connections)
    #membrane currents
    memb.currents <- as.matrix(read.table(svalue(membraneCurrent)))
    cmatr<-dim(memb.currents )
    memb.currents <-matrix(as.numeric(memb.currents),nrow=cmatr[1])
    #rm(lfpmatr)
    #length of segments
    seg.length <- as.numeric(as.matrix(read.table(svalue(SegLength))))
    #################################
    
    source('utils/kernel_basisfunctions_Regularizal.R',local=TRUE,verbose=TRUE)
   # source('alprogik/errorCalcFun.R')
    #source('/media/BA0ED4600ED416EB/agy/ksCSD_SVN/trunk/alprogik/egyesitett_ksCSDguihoz.R',local=TRUE)
    
    #Rprof('profilekernel_sim.out',memory.profiling=TRUE)
    sCSD_currents<-ksCSD_all( basis.width.min, basis.width.max, basis.width.step, basis.number.min, 
                              basis.number.max, basis.number.step, LFP , elec.kord ,memb.currents,seg.length  ,where2save) 
    #Rprof(NULL)
    image(t(sCSD_currents))
  }
  
  win <- gwindow("Calculation of Kernel Functions")
  gp <- ggroup(horizontal=FALSE, cont=win)
  
  ###
  #PLotting
  #ggraphics(ps=6, container=win)
  #     hist(rnorm(100))
  
  
  #Setting the parameters
  tmp <- gframe("Width of basis functions", container=gp, expand=TRUE)
  
  
  glabel("R min:", container=tmp)
  basis.width.min.g<-gedit("10", container=tmp, coerce.with=as.numeric)
  
  glabel("R max:", container=tmp)
  basis.width.max.g<-gedit("100", container=tmp, coerce.with=as.numeric)
  
  glabel("R step:", container=tmp)
  basis.width.step.g<-gedit("10", container=tmp, coerce.with=as.numeric)
  
  ##M number of basis function
  #M=2^k....
  tmp <- gframe("Number of basis functions M=2^k", container=gp, expand=TRUE)
  glabel("k min:", container=tmp)
  basis.number.min.g<-gedit("4", container=tmp, coerce.with=as.numeric)
  
  glabel("k max:", container=tmp)
  basis.number.max.g<-gedit("10", container=tmp, coerce.with=as.numeric)
  
  glabel("k step:", container=tmp)
  basis.number.step.g<-gedit("1", container=tmp, coerce.with=as.numeric)


  
  
  
  ###############################################
  ###########################################
  #LOading files
  #tmp <- gframe("Loading files", container=gp)#, expand=TRUE)
  #Electrode coordinates
  glabel("Electrode coordinates", container=gp)
  elec.cord.location <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/elcoord_x_y_z'), quote=FALSE, cont=gp)
  #cb <- gcombobox("", cont=g)
  
  #addHandlerChanged(elec.cord.location, handler=function(h,...) {
  #  elec.kord <- as.matrix(read.table(svalue(elec.cord.location )))
  #elec.kord<-matrix(elec.kord,ncol=3)
  #el.nb<-dim(elec.kord)[1] #number of electrodes
  #})
  #LFP
  glabel("LFP", container=gp)
  lfp.location <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/myLFP'), quote=FALSE, cont=gp)
  #cb <- gcombobox("", cont=g)
  
  #addHandlerChanged(lfp.location, handler=function(h,...) {
  #  LFP <- as.matrix(read.table(svalue(lfp.location )))
  #image(LFP)
  #})
  
  
  ##Location of membrane currents
  #Membrane currents
  glabel("Morphology or membrane currents location", container=gp)
  membrane.current.location <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/segcoordinates.txt'), quote=FALSE, cont=gp)
  #cb <- gcombobox("", cont=g)
  
  #addHandlerChanged(membrane.current.location , handler=function(h,...) {
  # morpho <- as.matrix(read.table(svalue(membrane.current.location)))
  #  cb[] <- colnames(x)
  #})
  ##
  ##Location of membrane currents
  #Membrane currents
  glabel("Connection information", container=gp)
  connection.info.location <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/connections.txt'), quote=FALSE, cont=gp)
  #cb <- gcombobox("", cont=g)
  
  #addHandlerChanged(connection.info.location , handler=function(h,...) {
  # connections <- as.matrix(read.table(svalue(connection.info.location )))
  #  cb[] <- colnames(x)
  #cat(connections)
  #})
  
  
  #Membrane currents, if known
  glabel("Membrane currents", container=gp)
  membraneCurrent <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/membcurr'), quote=FALSE, cont=gp)
  #membraneCurrent <- gfilebrowse("NA", quote=FALSE, cont=gp)
  
  #Length of the segments,if known
  
  glabel("Length of segments", container=gp)
  SegLength <- gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/seglength'), quote=FALSE, cont=gp)
  
  
  #cb <- gcombobox("", cont=g)
  
  #addHandlerChanged(membraneCurrent, handler=function(h,...) {
  # memb.currents <- as.matrix(read.table(svalue(membraneCurrent  )))
  #  cb[] <- colnames(x)
  #image(memb.currents)
  #})
  
  ######################################
  
  
  glabel("Where to save the results", container=gp)
  #simulation.location.name<- gedit("/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/simulation/proba",width = 50, container=gp)
  simulation.location.name<-gfilebrowse(paste0(wherearewenow, '/simulation/gang_9x9_200/'), 
                                        type='selectdir', quote=FALSE, cont=gp)
  
  #gfilebrowse (text = "Simulation location...", type = "selectdir", quote = TRUE, 
  #        container =gp) 
  
  
  
  #run simulation
  startSimulation<-gbutton("Start simulation", container=gp,
                           handler =kernel_calculate )
  
  
  
  f_exit <-gbutton("Close simulation", container=gp,
                   handler = function( h,...) dispose( win ))
  
  
  setwd(wherearewenow)
  
}
############
kernel_calc() 

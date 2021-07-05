plot_nifti<-function(name.long.temp,nodal_metrics,name.nii){
  vol.template<-readNIfTI(name.long.temp,reorient=F)
  new.vol<-vol.template
  location<-unique(as.vector(vol.template))
  location<-location[!is.na(location)]
  location<-location[location!=0]
  location<-sort(location)
  print(location)
  n.regions<-length(location)
  for(i in 1:n.regions){
    
    index<-which(vol.template==location[i],arr.ind=TRUE)
    new.vol[index]<-nodal_metrics[i]
    
  }
  #convert.datatype(datatype.code = 8)
  #browser()
  writeNIfTI(new.vol,name.nii,onefile=TRUE,gzipped=FALSE)
  
}

plot_niftilib<-function(name.long.temp,nodal_metrics,name.nii){
  vol.template<-nifti.image.read(name.long.temp)
  new.vol<-vol.template
  location<-unique(as.vector(vol.template))
  location<-location[!is.na(location)]
  location<-location[location!=0]
  location<-sort(location)
  print(location)
  n.regions<-length(location)
  for(i in 1:n.regions){
    
    index<-which(vol.template==location[i],arr.ind=TRUE)
    new.vol[index]<-nodal_metrics[i]
    
  }
  nifti.image.setdatatype(new.vol,16)
  nifti.set.filenames(new.vol,name.nii)
  
  nifti.image.write(new.vol)
  
}

plot_graphs <- function(adj.mat,name.dir,file_coord,regions_selected = NULL, dim.edges,name.png){
  
  coord<-paste(file_coord)
  
  set2 <- read.table(coord)
  set2 <- as.matrix(set2)
  n.regions<- dim(set2)[1]
  
  if(is.null(regions_selected)){
    regions_selected = seq_len(n.regions)
  }
  
  set2 <- set2[regions_selected,]
  n.regions<- dim(set2)[1]

  euclid <- 2 * dist(set2, method = "euclidean")
  
  x.euclid <- as.matrix(euclid)
  
  index <- c(1:(n.regions/2))*2 # regions gauches et droites
  
  set2[index,c(2,3)] <- set2[index,c(2,3)] + 1.5 #	decale les pts superposes sur le graphe
  
  CairoPNG(file=paste(name.dir,'/','graph_','_dim',dim.edges,'_',name.png,'_sagittal.png',sep=''),width = 1000, height = 800,family='times')
  par(mar=c(4.5, 0.1, 5, 0.3))
  
  x.coord<-2 #i
  y.coord<-3 # j
  
  plot(set2[,x.coord], set2[,y.coord], type = "p",xlab= "", ylab="",cex.lab=2,asp=1,ylim=c(min(set2[,y.coord]),max(set2[,y.coord])),xlim=c(min(set2[,x.coord]),max(set2[,x.coord])),xaxt='n',yaxt='n',bty='n',main=paste(dim.edges),cex=0.5,pch=16,cex.main=2) #points noirs
  
  ### text(labels.coord[1,],labels.coord[2,],labels.names,cex = 1)
  
  for(kk in 2:(n.regions)){
    for(q in 1:(kk-1)){
      
      if(adj.mat[kk,q]==1)
      {
        if(x.euclid[kk,q]>85) visu <- "blue" #longues/courtes arretes, A CHANGER
        else visu <- "red"
        lines(c(set2[kk,x.coord], set2[q,x.coord]), c(set2[kk,y.coord], set2[q,y.coord]), col = visu,lw=2)
      }
      
    }
  }
  dev.off()
  
  CairoPNG(file=paste(name.dir,'/','graph_','_dim',dim.edges,'_',name.png,'_transverse.png',sep=''),width = 1000, height = 800,family='times')
  par(mar=c(4.5, 0.1, 5, 0.3))
  
  x.coord<-1 #j
  y.coord<-2 #k
  
  plot(set2[,x.coord], set2[,y.coord], type = "p",xlab= "", ylab="",cex.lab=2,asp=1,ylim=c(min(set2[,y.coord]),max(set2[,y.coord])),xlim=c(min(set2[,x.coord]),max(set2[,x.coord])),xaxt='n',yaxt='n',bty='n',main=paste(dim.edges),cex=0.5,pch=16,cex.main=2) #points noirs
  
  
  #text(labels.coord[1,],labels.coord[2,],labels.names,cex = 2)
  
  for(kk in 2:(n.regions)){
    for(q in 1:(kk-1)){
      
      if(adj.mat[kk,q]==1)
      {
        if(x.euclid[kk,q]>85) visu <- "blue" 
        else visu <- "red"
        lines(c(set2[kk,x.coord], set2[q,x.coord]), c(set2[kk,y.coord], set2[q,y.coord]), col = visu,lw=2)
      }
    }
  }
  dev.off()
  
}

#############################################################################################################################################

plot_mvt <- function(path,spm_path){
  # from plot_mvt.R
  
  library('Cairo')
  library('R2HTML')
  
  
  Sujet <- ""
  ref<- ""
  
  print(paste(path,'/mvt_zoom',Sujet,ref,'.pdf',sep=''))
  pdf(paste(path,'/mvt_zoom',Sujet,ref,'.pdf',sep=''))
  par(mfrow=c(2,1))	
  
  
  path.mov<-list.files(path=spm_path,pattern=glob2rx("rp*.txt"))
  mvt <- read.table(paste(spm_path,path.mov,sep='/'), header=FALSE)# fichiers des param de mvts
  
  ymin<-min(mvt[,1:3])-0.05
  ymax<-max(mvt[,1:3])+0.05
  
  plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
  lines(mvt[,2],type='l',col='green')
  lines(mvt[,3],type='l',col='red')
  legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
  
  ymin<-min(mvt[,4:6]*180/pi)-0.05
  ymax<-max(mvt[,4:6]*180/pi)+0.05
  plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
  lines(mvt[,5]*180/pi,type='l',col='green')
  lines(mvt[,6]*180/pi,type='l',col='red')
  legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
  dev.off()
  
  
  png(paste(path,'/mvt_zoom',Sujet,ref,'.png',sep=''))
  par(mfrow=c(2,1))	
  
  
  ymin<-min(mvt[,1:3])-0.05
  ymax<-max(mvt[,1:3])+0.05
  
  plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
  lines(mvt[,2],type='l',col='green')
  lines(mvt[,3],type='l',col='red')
  legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
  
  ymin<-min(mvt[,4:6]*180/pi)-0.05
  ymax<-max(mvt[,4:6]*180/pi)+0.05
  plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
  lines(mvt[,5]*180/pi,type='l',col='green')
  lines(mvt[,6]*180/pi,type='l',col='red')
  legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
  dev.off()
  
  
  
  CairoPDF(file=paste(path,'/mvt',ref,'.pdf',sep=''),width = 15, height = 11)
  par(mfrow=c(1,2))
  print(paste(path,'/mvt',Sujet,ref,'.pdf',sep=''))
  
  ymin<- -20
  ymax<-20
  
  plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
  lines(mvt[,2],type='l',col='green')
  lines(mvt[,3],type='l',col='red')
  vec<-rep(1,400)
  lines(2*vec,type='l',col='cyan')
  lines(-2*vec,type='l',col='cyan')
  legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
  
  ymin<- -20
  ymax<- 20
  
  plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
  lines(mvt[,5]*180/pi,type='l',col='green')
  lines(mvt[,6]*180/pi,type='l',col='red')
  lines(2*vec,type='l',col='cyan')
  lines(-2*vec,type='l',col='cyan')
  legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
  dev.off()
  
  
  CairoPNG(file=paste(path,'/mvt',ref,'.png',sep=''),width = 600, height = 500)
  par(mfrow=c(1,2))
  print(paste(path,'/mvt',Sujet,ref,'.png',sep=''))
  
  ymin<- -10
  ymax<-10
  
  plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
  lines(mvt[,2],type='l',col='green')
  lines(mvt[,3],type='l',col='red')
  vec<-rep(1,400)
  lines(2*vec,type='l',col='cyan')
  lines(-2*vec,type='l',col='cyan')
  legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
  
  ymin<- -10
  ymax<- 10
  
  plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
  lines(mvt[,5]*180/pi,type='l',col='green')
  lines(mvt[,6]*180/pi,type='l',col='red')
  lines(2*vec,type='l',col='cyan')
  lines(-2*vec,type='l',col='cyan')
  legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
  dev.off()
}

#############################################################################################################################################

plot_graphs_old <- function(adj.mat,path,templBaseName,file_coord,regions_selected = NULL, dim.edges,name.png){
  
  coord<-paste(file_coord)
  
  set2 <- read.table(coord)
  set2 <- as.matrix(set2)
  n.regions<- dim(set2)[1]
  
  if(is.null(regions_selected)){
    regions_selected = seq_len(n.regions)
  }
  
  set2 <- set2[regions_selected,]
  n.regions<- dim(set2)[1]
  
  euclid <- 2 * dist(set2, method = "euclidean")
  
  x.euclid <- as.matrix(euclid)
  
  index <- c(1:(n.regions/2))*2 # regions gauches et droites
  
  set2[index,c(2,3)] <- set2[index,c(2,3)] + 1.5 #	decale les pts superposes sur le graphe
  
  name.dir<-paste(path,'Graph_Measures',templBaseName,'',sep='/') # sortie
  
  CairoPNG(file=paste(name.dir,'/','graph_','_dim',dim.edges,'_',name.png,'_sagittal.png',sep=''),width = 1000, height = 800,family='times')
  par(mar=c(4.5, 0.1, 5, 0.3))
  
  x.coord<-2 #i
  y.coord<-3 # j
  
  plot(set2[,x.coord], set2[,y.coord], type = "p",xlab= "", ylab="",cex.lab=2,asp=1,ylim=c(min(set2[,y.coord]),max(set2[,y.coord])),xlim=c(min(set2[,x.coord]),max(set2[,x.coord])),xaxt='n',yaxt='n',bty='n',main=paste(dim.edges),cex=0.5,pch=16,cex.main=2) #points noirs
  
  ### text(labels.coord[1,],labels.coord[2,],labels.names,cex = 1)
  
  for(kk in 2:(n.regions)){
    for(q in 1:(kk-1)){
      
      if(adj.mat[kk,q]==1)
      {
        if(x.euclid[kk,q]>85) visu <- "blue" #longues/courtes arretes, A CHANGER
        else visu <- "red"
        lines(c(set2[kk,x.coord], set2[q,x.coord]), c(set2[kk,y.coord], set2[q,y.coord]), col = visu,lw=2)
      }
      
    }
  }
  dev.off()
  
  CairoPNG(file=paste(name.dir,'/','graph_','_dim',dim.edges,'_',name.png,'_transverse.png',sep=''),width = 1000, height = 800,family='times')
  par(mar=c(4.5, 0.1, 5, 0.3))
  
  x.coord<-1 #j
  y.coord<-2 #k
  
  plot(set2[,x.coord], set2[,y.coord], type = "p",xlab= "", ylab="",cex.lab=2,asp=1,ylim=c(min(set2[,y.coord]),max(set2[,y.coord])),xlim=c(min(set2[,x.coord]),max(set2[,x.coord])),xaxt='n',yaxt='n',bty='n',main=paste(dim.edges),cex=0.5,pch=16,cex.main=2) #points noirs
  
  
  #text(labels.coord[1,],labels.coord[2,],labels.names,cex = 2)
  
  for(kk in 2:(n.regions)){
    for(q in 1:(kk-1)){
      
      if(adj.mat[kk,q]==1)
      {
        if(x.euclid[kk,q]>85) visu <- "blue" 
        else visu <- "red"
        lines(c(set2[kk,x.coord], set2[q,x.coord]), c(set2[kk,y.coord], set2[q,y.coord]), col = visu,lw=2)
      }
    }
  }
  dev.off()
  
}



plot_mvt_old <- function(path,templBaseName,spm_path){
# from plot_mvt.R
	
  library('Cairo')
  library('R2HTML')
	
	
	Sujet <- ""
	ref<- ""

	name.dir<-paste(path,'Graph_Measures',templBaseName,'',sep='/') 

	print(paste(path,'/Graph_Measures/',templBaseName,'/mvt_zoom',Sujet,ref,'.pdf',sep=''))
	pdf(paste(path,'/Graph_Measures/',templBaseName,'/mvt_zoom',Sujet,ref,'.pdf',sep=''))
	par(mfrow=c(2,1))	
	

	path.mov<-list.files(path=spm_path,pattern=glob2rx("rp*.txt"))
	mvt <- read.table(paste(spm_path,path.mov,sep='/'), header=FALSE)# fichiers des param de mvts
	
	ymin<-min(mvt[,1:3])-0.05
	ymax<-max(mvt[,1:3])+0.05
	
	plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
	lines(mvt[,2],type='l',col='green')
	lines(mvt[,3],type='l',col='red')
	legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
	
	ymin<-min(mvt[,4:6]*180/pi)-0.05
	ymax<-max(mvt[,4:6]*180/pi)+0.05
	plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
	lines(mvt[,5]*180/pi,type='l',col='green')
	lines(mvt[,6]*180/pi,type='l',col='red')
	legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
	dev.off()
	

	png(paste(path,'/Graph_Measures/',templBaseName,'/mvt_zoom',Sujet,ref,'.png',sep=''))
	par(mfrow=c(2,1))	
	
		
	ymin<-min(mvt[,1:3])-0.05
	ymax<-max(mvt[,1:3])+0.05
	
	plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
	lines(mvt[,2],type='l',col='green')
	lines(mvt[,3],type='l',col='red')
	legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
	
	ymin<-min(mvt[,4:6]*180/pi)-0.05
	ymax<-max(mvt[,4:6]*180/pi)+0.05
	plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
	lines(mvt[,5]*180/pi,type='l',col='green')
	lines(mvt[,6]*180/pi,type='l',col='red')
	legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
	dev.off()
	


	CairoPDF(file=paste(path,'/Graph_Measures/',templBaseName,'/mvt',ref,'.pdf',sep=''),width = 15, height = 11)
	par(mfrow=c(1,2))
	print(paste(path,'/Graph_Measures/',templBaseName,'/mvt',Sujet,ref,'.pdf',sep=''))
	
	ymin<- -20
	ymax<-20
	
	plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
	lines(mvt[,2],type='l',col='green')
	lines(mvt[,3],type='l',col='red')
	vec<-rep(1,400)
	lines(2*vec,type='l',col='cyan')
	lines(-2*vec,type='l',col='cyan')
	legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
	
	ymin<- -20
	ymax<- 20
	
	plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
	lines(mvt[,5]*180/pi,type='l',col='green')
	lines(mvt[,6]*180/pi,type='l',col='red')
	lines(2*vec,type='l',col='cyan')
	lines(-2*vec,type='l',col='cyan')
	legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
	dev.off()


	CairoPNG(file=paste(path,'/Graph_Measures/',templBaseName,'/mvt',ref,'.png',sep=''),width = 600, height = 500)
	par(mfrow=c(1,2))
	print(paste(path,'/Graph_Measures/',templBaseName,'/mvt',Sujet,ref,'.png',sep=''))
	
	ymin<- -10
	ymax<-10
	
	plot(mvt[,1],type='l',col='blue',ylim=c(ymin,ymax),ylab='mm',xlab='image',main=paste(ref,'\nTranslation',sep=''))
	lines(mvt[,2],type='l',col='green')
	lines(mvt[,3],type='l',col='red')
	vec<-rep(1,400)
	lines(2*vec,type='l',col='cyan')
	lines(-2*vec,type='l',col='cyan')
	legend(x='topleft',legend=c('x','y','z'),col=c('blue','green','red'),lty=1)
	
	ymin<- -10
	ymax<- 10
	
	plot(mvt[,4]*180/pi,type='l',col='blue',ylim=c(ymin,ymax),ylab='degrees',xlab='image',main=paste('Rotation'))
	lines(mvt[,5]*180/pi,type='l',col='green')
	lines(mvt[,6]*180/pi,type='l',col='red')
	lines(2*vec,type='l',col='cyan')
	lines(-2*vec,type='l',col='cyan')
	legend(x='topleft',legend=c('pitch','roll','yaw'),col=c('blue','green','red'),lty=1)
	dev.off()



  }

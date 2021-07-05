#!/usr/bin/env Rscript
computeTS_annotated <- function(rs_path, 
                                r_dyn_name, 
                                art_path, 
                                atlas_file, 
                                gm_file, 
                                graph_path){
  
  vol.template<-readNIfTI(atlas_file,reorient=F)
  
  if(length(dim(vol.template))==4) vol.template<-vol.template[,,,1]
  
  location<-unique(as.vector(vol.template))
  location<-location[!is.na(location)]
  location<-location[location!=0]
  location<-sort(location)
  print(location)
  n.regions<-length(location)	
  
  #name.grey.matter<-list.files(path=T1_path,pattern=glob2rx(paste('^rc1*.nii', sep='')))
  
  vol.grey.matter<-readNIfTI(gm_file,reorient=F)
  
  #list.of.in.files<-list.files(path=spm_path,pattern=glob2rx(paste(Func_BaseName, '*.nii', sep='')))
  list.of.in.files <- list.files(rs_path, pattern=glob2rx(r_dyn_name))
  cat(list.of.in.files[1],'\n')
  
  length.proc<-length(list.of.in.files)
  
  # initialization
  print('Initialisation folders')
  templBaseName <- file_path_sans_ext(basename(atlas_file))
  name.ts <- templBaseName
  
  #Create directories
  dir.create(file.path(graph_path,"Functional"), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/',sep='/'), "data"), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/data',sep='/'), templBaseName), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/',sep='/'), "grey_matter_data"), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/grey_matter_data',sep='/'), templBaseName), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/',sep='/'), "index"), showWarnings = FALSE)
  dir.create(file.path(paste(graph_path,'Functional/index',sep='/'), templBaseName), showWarnings = FALSE)
  
  print('Extracting Time Series...')
  data.ts.gm<-matrix(0,length.proc,n.regions)
  
  
  for(i in 1:n.regions){
    
    name.txt<-paste(paste(graph_path,'Functional/data',templBaseName,'',sep='/'),name.ts,'_voxels_time_series_region_',i,'.txt',sep='')
    
    write.table(paste("Time series of each voxel for region",i,' ',sep=''),name.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=" ")
    write.table(" ",name.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    name.matter.txt<-paste(paste(graph_path,'Functional/grey_matter_data',templBaseName,'',sep='/'),'grey_matter_',name.ts,'_voxels_time_series_region_',i,'.txt',sep='')
    
    write.table(paste("Grey matter coefficients for region",i,' ',sep=''),name.matter.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=" ")
    write.table(" ",name.matter.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    
    index<-which(vol.template==location[i],arr.ind=TRUE)
    size.r<-dim(index)[1]
    tmp<-vol.grey.matter[index]
    tmp[is.na(tmp)]<-0
    write.table(c(tmp),name.matter.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,eol=" ")
    write.table(" ",name.matter.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    
    name.index.txt<-paste(paste(graph_path,'Functional/index',templBaseName,'',sep='/'),'index_',name.ts,'_voxels_time_series_region_',i,'.txt',sep='')
    
    write.table(paste("Voxels coordinates for region",i,sep=''),name.index.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=" ")
    write.table(" ",name.index.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    
    for(j in 1:3){
      write.table(c(index[,j]),name.index.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=" ",append=TRUE)
      write.table(" ",name.index.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    }    
  }
  
  # time series
  print('ts writting')
  
  for(i in 1:length.proc){	# loop on time
    cat(i,'\n')
    
    vol<-readNIfTI(file.path(rs_path,list.of.in.files[i]),reorient=F)
    
    for(j in 1:n.regions){	# loop on region
      
      index<-which(vol.template==location[j],arr.ind=TRUE)
      size.r<-dim(index)[1]
      
      name.txt<-paste(paste(graph_path,'Functional/data',templBaseName,'',sep='/'),name.ts,'_voxels_time_series_region_',j,'.txt',sep='')
      ##cat(vol[index],'/n')
      # TS for each voxel of each region
      write.table(c(i,vol[index]),name.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE,eol=" ")
      write.table(" ",name.txt,quote=FALSE,col.names=FALSE,row.names=FALSE,eol="\n",append=TRUE)
    }
  }
  
  data.ts<-matrix(0,length.proc,n.regions)
  data.ts.gm<-matrix(0,length.proc,n.regions)
  
  for(i in 1:n.regions){	# loop on regions
    cat('Regions ',i,'\n')
    name.txt<-paste(paste(graph_path,'Functional/data',templBaseName,'',sep='/'),name.ts,'_voxels_time_series_region_',i,'.txt',sep='')
    tmp<-read.table(name.txt,header=FALSE,skip=1)
    data<-as.matrix(tmp)[,-1]
    nb.voxels<-dim(data)[2]
    
    name.grey.matter<-paste(paste(graph_path,'Functional/grey_matter_data',templBaseName,'',sep='/'),'grey_matter_',name.ts,'_voxels_time_series_region_',i,'.txt',sep='')
    
    coef.grey.matter<-read.table(name.grey.matter,header=FALSE,skip=1)
    
    if(is.vector(data)){
      nb.voxels<-length(data)
      data.ts[,i]<-sum(data)/nb.voxels
    }
    if(is.matrix(data)){
      nb.voxels<-dim(data)[2]
      
      data.ts[,i]<-rowSums(data)/nb.voxels
    }	 
    data.ts.gm[,i]<-(data%*%t(coef.grey.matter))/nb.voxels
    
  }
  
  write.table(data.ts.gm,paste(graph_path,'/Functional/data/',templBaseName,'/',name.ts,'_ts_raw.txt',sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # REGRESSION	
  set1<-data.ts.gm  # time series changed to data.ts.gm, data.ts in the original
  N <- length(set1[,1])   # length of the time series
  print(c("ts length :",N))
  
  
  print('Regression: Movement correction...')
  # Movement correction
  graph_path.mov<-list.files(path=art_path,pattern=glob2rx("art_regression_outliers_and_movement_*.mat"))
  cat(c(graph_path),'\n')
  cat(c(graph_path.mov),'\n')
  set2 <- readMat(paste(art_path,graph_path.mov,sep='/'), header=FALSE)
  set2<-set2[[1]]
  set2 <- as.data.frame(set2) #art mat values
  
  set3 <- rbind(0, set2[-N,]) #delete last row and and first 0 row
  set3 <- cbind(set2, set2 - set3) #set3 is a reorganisation of data in set2
  names(set3) <- c(paste0("p",seq_len(dim(set2)[2])),paste0("p",seq_len(dim(set2)[2]),"d"))
  regressors<-cbind(set3)
  write.table(as.matrix(regressors),paste(graph_path,'/Functional/data/',templBaseName,'/Regressors.txt',sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  set1[is.na(set1)] <- 0 #set for the moment
  set4 <- resid(lm(as.matrix(set1) ~ as.matrix(regressors), na.action=na.omit)) #set4 gets the residuals of the linear regression between the ts and the movement data, WM, CSF and outliers
  set4 <- as.data.frame(set4)	
  
  data.correct<-set4
  write.table(data.correct,paste(graph_path,'/Functional/data/',templBaseName,'/',name.ts,'_ts.txt',sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE)
}
#computeTS(graph_path,spm_path,T1_path,templBaseName, name.long.temp)
#computeTS <- function(graph_path,spm_path,T1_path,templBaseName,name.long.temp)

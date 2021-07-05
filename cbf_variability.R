library("pracma")
library("neurobase")
library("oro.nifti")
library("ggplot2")

controlPath<-"/media/veronica/DATAPART2/EmoPark/Data/Controls/"
patientPath<-"/media/veronica/DATAPART2/EmoPark/Data/Patients/"
atlasPath <- "/usr/local/MATLAB/spm12/atlas/AAL1_brain.nii"

controlList <- dir(controlPath)
controlList <- controlList[controlList != "213"]
patientList <- dir(patientPath)
patientList <- patientList[patientList != "PAT15"]

atlas <- readNIfTI(atlasPath
  , reorient = FALSE
  , read_data = TRUE
)
atlas <- atlas@.Data
regions <- unique(as.vector(atlas))
regions <- regions[order(regions)][-1]
n.regions.max <- length(regions)

cbf_c <- vector(mode = "list", length = n.regions.max)
names(cbf_c) <- paste0("R",as.character(1:n.regions.max))
cbf_p <- vector(mode = "list", length = n.regions.max)
names(cbf_p) <- paste0("R",as.character(1:n.regions.max))
  #c("Frontal","Parietal","Temporal","Occipital","WM","CingIns","BG","RN","SN","STN","Cau","Put","GPe","GPi","Th","CS","CI")

for (i in 1:length(controlList)){
  if (file.exists(paste0(controlPath, controlList[i],"/pCASL/rCBF_resampled.nii"))){
    clust <- readNIfTI(
      paste0(controlPath, controlList[i],"/pCASL/rCBF_resampled.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Clust<-clust@.Data
    Clust[is.na(Clust)]<-0
    Clust <- Clust[is.finite(Clust)]
    
    atlas <- readNIfTI(
      paste0(controlPath, controlList[i],"/wAAL1_brain.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Atlas<-atlas@.Data
    
    for ( r in 1:n.regions.max ){
      cbf_c[[r]] <- c( cbf_c[[r]], Clust[Atlas == regions[r]] )
      #cbf_c[[r]][[i]] <- Clust[Atlas == regions[r]]
    }
  }
}

for (i in 1:length(patientList)){
  if (file.exists(paste0(patientPath, patientList[i],"/pCASL/rCBF_resampled.nii"))){
    clust <- readNIfTI(
      paste0(patientPath, patientList[i],"/pCASL/rCBF_resampled.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Clust<-clust@.Data
    Clust[is.na(Clust)]<-0
    Clust <- Clust[is.finite(Clust)]
    
    atlas <- readNIfTI(
      paste0(patientPath, patientList[i],"/wAAL1_brain.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Atlas<-atlas@.Data
    
    for ( r in 1:n.regions.max ){
      cbf_p[[r]][[i]] <- Clust[Atlas == regions[r]] 
    }
  }
}

p_tab <- data.frame(array(NA,dim = c(90,length(patientList))), row.names = regions)
colnames(p_tab) <- patientList

for ( r in 1:n.regions.max ){
  for ( i in 1:length(patientList) ){
    p_tab[r,i] <- wilcox.test(unlist(cbf_c[[r]]),unlist(cbf_p[[r]][i]))$p.value
  }
}

p_tab_cor <- apply(p_tab,MARGIN = 1, FUN = function(x) p.adjust(x,'bonferroni'))
p_tab_cor<-t(p_tab_cor)
p_tab_cor <- round(p_tab_cor,2)

a <- p_tab_cor
a[a>0.05]<-NA
#a<-a[rowSums(is.na(a)) != ncol(a), ]

library(plot.matrix)
par(mar=c(5.1, 4.1, 4.1, 4.1), las=1, cex.axis = 0.8)
plot(a, main = "P_values < 0.05 (corrected by # of subjects)", xlab = "", ylab="Atlas regions")




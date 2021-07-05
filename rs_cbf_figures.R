library("pracma")
library("neurobase")
library("oro.nifti")
library("ggplot2")

controlPath<-"/media/veronica/DATAPART2/EmoPark/Data/Controls/"
patientPath<-"/media/veronica/DATAPART2/EmoPark/Data/Patients/"
atlasPath <- "/usr/local/MATLAB/spm12/atlas/AAL1_brain.nii"
regions_selected <- NULL

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
n.regions <- length(regions)

if (!is.null(regions_selected)){
  regions <- regions[regions_selected]
  n.regions <- length(regions)
}

# --------------- CONTROLS ---------------
cbf_c2 <- data.frame(group = factor(),
                    subject = factor(),
                    region = factor(),
                    cbf = integer())

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
    
    for ( r in 1:n.regions ){
      cbf <- Clust[Atlas == regions[r]]
      group <- rep("Control", length(cbf))
      subject <- rep(controlList[i], length(cbf))
      region <- rep( paste0("R",as.character(regions[r])), length(cbf))
      cbf_c2 <- rbind(cbf_c2, data.frame(group,subject,region,cbf))
    }
  }
}


# -------------------- PATIENTS -------------------------

cbf_p2 <- data.frame(group = factor(),
                    subject = factor(),
                    region = factor(),
                    cbf = integer())

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
    
    for ( r in 1:n.regions ){
      cbf <- Clust[Atlas == regions[r]]
      group <- rep("Patient", length(cbf))
      subject <- rep(patientList[i], length(cbf))
      region <- rep( paste0("R",as.character(regions[r])), length(cbf))
      cbf_p2 <- rbind(cbf_p2, data.frame(group,subject,region,cbf))
    }
  }
}

cbf_all2 <- rbind(cbf_c2,cbf_p2) 

## ------------------------------ PLOTS -------------------------------------

# Distribution of CBF values from controls ----------------------------------
ggplot(data = cbf_c[cbf_c$subject!="213",], aes(x= cbf, fill=subject)) +
  #ggplot(data = cbf_c, aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  labs(x="CBF",
       title="CBF distribution values in controls (without 213)") +
  theme_classic()


# Distribution of CBF values from controls ----------------------------------
ggplot(data = cbf_p[!(cbf_p$subject %in% c("PAT15")),], aes(x= cbf, fill=subject)) +
#ggplot(data = cbf_p, aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  scale_fill_brewer(palette="Set3") +
  labs(x="CBF",
       title="CBF distribution values in patients (without PAT15)") +
  theme_classic()

# Distribution of CBF values from controls and patients ----------------------------------
ggplot(data = cbf_all[!(cbf_all$subject %in% c("213","PAT15")),], aes(x= cbf, fill=subject)) +
  #ggplot(data = cbf_p, aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  labs(x="CBF",
       title="CBF distribution values all subjects") +
  theme_classic()

# Distribution of CBF values controls vs patients ----------------------------------
ggplot(data = cbf_all[!(cbf_all$subject %in% c("213","PAT15")),], aes(x= cbf, fill=group)) +
  #ggplot(data = cbf_p, aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  labs(x="CBF",
       title="CBF distribution values in controls and patients") +
  theme_classic()


cbf_test <- cbf_all[!(cbf_all$subject %in% c("213","PAT15")),]
levels(cbf_test$subject) <- c(levels(cbf_test$subject), "All_controls")
cbf_test$subject[cbf_test$group == "Control"] <- "All_controls"

# Distribution of CBF values all controls vs individual patients ----------------------------------
library(RColorBrewer)
mycolors <- c(brewer.pal(12,"Set3"),"#333333")
ggplot(data = cbf_test, aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values = mycolors) +
  labs(x="CBF",
       title="CBF distribution values all controls versus individual subjects") +
  theme_classic()


# Distribution of CBF values all controls vs individual patients in a specific region ----------------------------------
library(RColorBrewer)
mycolors <- c(brewer.pal(12,"Set3"),"#333333")
ggplot(data = cbf_test[cbf_test$region == "R2001",], aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values = mycolors) +
  labs(x="CBF",
       title="PRECENTRAL_L - CBF distribution values all controls versus individual subjects") +
  theme_classic()



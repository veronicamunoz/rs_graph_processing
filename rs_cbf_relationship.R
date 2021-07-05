library("pracma")
library("neurobase")
library("oro.nifti")
library("ggplot2")

controlDataPath<-"/media/veronica/DATAPART2/EmoPark/Data/Controls/"
patientDataPath<-"/media/veronica/DATAPART2/EmoPark/Data/Patients/"
controlGraphPath<-"/media/veronica/DATAPART2/EmoPark/Graphs/Controls/"
patientGraphPath<-"/media/veronica/DATAPART2/EmoPark/Graphs/Patients/"
atlasPath <- "/usr/local/MATLAB/spm12/atlas/AAL1_brain.nii"

controlList <- dir(controlDataPath)
controlList <- controlList[controlList != "213"]
patientList <- dir(patientDataPath)
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
controls.df <- data.frame(group = factor(),
                    subject = factor(),
                    region = factor(),
                    cbf_mean = integer(),
                    cbf_var = integer(),
                    cbf_std = integer(),
                    rs_mean = integer(),
                    rs_var = integer(),
                    rs_std = integer())

for (i in 1:length(controlList)){
  cbf_path <- paste0(controlDataPath, controlList[i],"/pCASL/rCBF_resampled.nii")
  #rs_corr_path <- paste0(controlGraphPath, controlList[i],"/Graph_Brut/wAAL1_brain/wave.cor.mat_n.levels_3_n.regions_90_start_time_series_1proc_length500.grey.matter.txt")
  rs_corr_path <-paste0(controlGraphPath, controlList[i], "/Graph_Brut/wAAL3_mod_v1/wave.cor.mat_n.regions_100_start_time_series_1proc_length500.grey.matter.txt")
  
  if (file.exists(cbf_path) & file.exists(rs_corr_path)){
    clust <- readNIfTI(
      cbf_path
      , reorient = FALSE
      , read_data = TRUE
    )
    Clust<-clust@.Data
    Clust[is.na(Clust)]<-0
    Clust <- Clust[is.finite(Clust)]
    
    atlas <- readNIfTI(
      paste0(controlDataPath, controlList[i],"/wAAL1_brain.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Atlas<-atlas@.Data
    
    rs_corr <- read.table(rs_corr_path)
    rs_corr <- as.matrix(rs_corr)
    
    for ( r in 1:n.regions ){
      r_data <- data.frame(group = "Control",
                           subject = controlList[i],
                           region = paste0("R",as.character(regions[r])),
                           cbf_mean = mean(Clust[Atlas == regions[r]]),
                           cbf_var = var(Clust[Atlas == regions[r]]),
                           cbf_std = std(Clust[Atlas == regions[r]]),
                           rs_mean = mean(rs_corr[r,]),
                           rs_var = var(rs_corr[r,]),
                           rs_std = std(rs_corr[r,]))
                           
      controls.df <- rbind(controls.df, r_data)
    }
  }
}


# -------------------- PATIENTS -------------------------
patients.df <- data.frame(group = factor(),
                          subject = factor(),
                          region = factor(),
                          cbf_mean = integer(),
                          cbf_var = integer(),
                          cbf_std = integer(),
                          rs_mean = integer(),
                          rs_var = integer(),
                          rs_std = integer())

for (i in 1:length(patientList)){
  cbf_path <- paste0(patientDataPath, patientList[i],"/pCASL/rCBF_resampled.nii")
  #rs_corr_path <- paste0(patientGraphPath, patientList[i],"/Graph_Measures/wAAL1_brain/wave.cor.mat_n.levels_3_n.regions_90_start_time_series_1proc_length500.grey.matter.txt")
  rs_corr_path <-paste0(patientGraphPath, patientList[i], "/Graph_Brut/wAAL3_mod_v1/wave.cor.mat_n.regions_100_start_time_series_1proc_length500.grey.matter.txt")
  
  if (file.exists(cbf_path) & file.exists(rs_corr_path)){
    clust <- readNIfTI(
      cbf_path
      , reorient = FALSE
      , read_data = TRUE
    )
    Clust<-clust@.Data
    Clust[is.na(Clust)]<-0
    Clust <- Clust[is.finite(Clust)]
    
    atlas <- readNIfTI(
      paste0(patientDataPath, patientList[i],"/wAAL1_brain.nii")
      , reorient = FALSE
      , read_data = TRUE
    )
    Atlas<-atlas@.Data
    
    rs_corr <- read.table(rs_corr_path)
    rs_corr <- as.matrix(rs_corr)
    
    for ( r in 1:n.regions ){
      r_data <- data.frame(group = "Patient",
                           subject = patientList[i],
                           region = paste0("R",as.character(regions[r])),
                           cbf_mean = mean(Clust[Atlas == regions[r]]),
                           cbf_var = var(Clust[Atlas == regions[r]]),
                           cbf_std = std(Clust[Atlas == regions[r]]),
                           rs_mean = mean(rs_corr[r,]),
                           rs_var = var(rs_corr[r,]),
                           rs_std = std(rs_corr[r,]))
      
      patients.df <- rbind(patients.df, r_data)
    }
  }
}

DF <- rbind(controls.df,patients.df)
## --------------------PLOTS --------------------------------------

library(ggplot2)

ggplot(data = DF, aes(x = cbf_var, y = rs_mean)) + 
  geom_point(aes(colour=region, shape=group)) + 
  theme_classic()


ggplot(data = controls.df, aes(x = cbf_std, y = rs_mean, colour = region)) + 
  geom_point(shape=group) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  theme_classic()

ggplot(data = controls.df, aes(x = cbf_std, y = rs_mean)) + 
  geom_line(aes(colour=region)) + 
  theme_classic()


r_list <- (seq_len(45)*2)+1

ggplot(data = patients.df[patients.df$region %in% patients.df$region[r_list],], aes(x = cbf_std, y = rs_mean, colour = region)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  labs(title = "AAL3 Right brain regions - mean RS_corr =  f(std CBF)") + 
  theme_classic()


ggplot(data = DF[DF$region %in% c("R7022"),], aes(x = cbf_std, y = rs_mean, colour = group)) + 
  geom_point() + 
  geom_text(aes(label=subject),hjust=0, vjust=0) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  labs(title = "PALLIDIUM_R - mean RS_corr =  f(var CBF)") + 
  theme_classic()


library(RColorBrewer)
mycolors <- c(brewer.pal(12,"Set3"),"#333333")
ggplot(data = cbf_test[cbf_test$region == "R2001",], aes(x= cbf, fill=subject)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values = mycolors) +
  labs(x="CBF",
       title="PRECENTRAL_L - CBF distribution values all controls versus individual subjects") +
  theme_classic()




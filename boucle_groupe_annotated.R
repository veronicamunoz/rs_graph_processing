library('bitops')
library('tools')
library('oro.nifti')
library('R.matlab')
source('Figure.R')
source('computeTS_annotated.R')
source('computeGraphs_annotated.R')
source('computeGraphs_brut.R')

#list_control<-dir(path =paste("/media/achards/Samsung_T5/DATA_REF_PATIENTS/",sep=''),pattern=glob2rx('^patient*'), full.names = FALSE, recursive = FALSE)
data_path <- "/media/veronica/DATAPART2/EmoPark/Data/Controls"
list_control <- dir(path = data_path)
list_control <- list_control[list_control!="PAT06"]
## list of subjects in order to create corresponding directory and compute tim series 
## you should mention the directory where you put the list of subjects and the pattern to make an automatic list

res_path <-'/media/veronica/DATAPART2/EmoPark/Graphs/Controls'
##this should contain the path where you want to write the results of the graphs 

file_coord<- '/media/veronica/DATAPART2/EmoPark/Data/coord_AAL3_mod_v1.txt' ### coordinate of the parcels  should contain same number of regions as the time series 

for(name.control in list_control){
  cat(name.control,'\n')
  
  Func_BaseName <- 'rra' # to change corresponding to the fMRI data to select
  
  graph_path <- file.path(res_path, name.control)#results
  dir.create(graph_path, showWarnings = FALSE) ### create the directory to write the results for each subject 
  
  #spm_path<-dir(path =paste("/Volumes/DATASAVE/DATA_REF_PATIENTS/",name.control,"/Processed/",sep=''),pattern=glob2rx('*Function*'), full.names = TRUE, recursive = FALSE) #find path for functional data
  
  rs_path <- file.path(data_path, name.control, 'RS') ## path for realigned data
  #r_dyn_files <- dir(path = rs_path ,pattern=glob2rx('*rRS_MB3_3x3x3_135mmFH_500dyn_*'), full.names = TRUE, recursive = FALSE)
  r_dyn_name <- '*rRS_MB3_3x3x3_135mmFH_500dyn_*'
  
  art_path <- file.path(data_path, name.control, 'RS') ## path for art outputs
  
  T1_path <- file.path(data_path, name.control, 'Anat') # path to anatomical data
  
  atlas_file <- file.path(data_path, name.control, 'wAAL3_mod_v1.nii')
  gm_file <- file.path(data_path, name.control, 'rmwp1T1_3D.nii')

  # templBaseName <- 'AAL3v1' ### name of the template used to process SPM, used in computeTS_annotated to read the nifti file for parcels
  # name.long.temp<- 'AAL3v1' ### 
  
  
  #### Step for computation of time series 
  #computeTS_annotated(rs_path, r_dyn_name, art_path, atlas_file, gm_file, graph_path) ### pour calculer les series temporelles 
  #graph_file <- file.path(graph_path, 'Graph_Measures', 'wAAL3_mod_v1_5per','BC_n.levels_4_n.regions_108_start_time_series_1proc_length500.mst.grey.matter.txt')
  
  #if (!file.exists(graph_file)) {
  #### Step for computation of graphs
  res<-compute_Graph_brut(name.control,
                             rs_path,
                              atlas_file,
                              graph_path,
                              file_coord,
                              graphs=TRUE,
                              regions_selected = c(1:90, 97:106),
                              percentage_selected = 3,
                              num.levels=3,
                              num.nb.edges = seq(200,3000,200),
                              length_time_series=NULL
                              )
  ##browser()

  }

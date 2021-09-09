%--------------------------------------------------------------------------
%
%              rs-fMRI PROCESSING PIPELINES FOR GRAPH ANALYSIS
%                           (Grenoble, FR, 2021)
%
% Pre-processing by Veronica Munoz Ramirez
%                   Michel Dojat
%                   Chantal Delon-Martin  
% Time-series & graph extraction by Sophie Achard
%
%
% The code assumes that the Path folder contains one folder per subject
% with, at least, an anatomical (3D) and resting state (4D) MR nifti files
% corresponding to the individual. In practive the anat and resting state
% files have the same name for all subjects.
% 
%--------------------------------------------------------------------------

global spm_path
spm_path = '/usr/local/MATLAB/spm12';
addpath(spm_path);

art_toolbox_path = '/usr/local/MATLAB/spm12/toolbox/conn'; % artifact detection toolbox
addpath(art_toolbox_path);

Path = '/media/veronica/DATAPART2/EmoPark/Data/Controls_test/'; 

%% SEGMENTATION
% Brain segmentation using SPM/CAT12 functions
anat_file = 'Anat/T1_3D.nii'; 
anat_segmentation(Path,anat_file,2); % 1st : Path to files, 2th: How many simultaneous processes ? 
disp('We need to wait until all segmentations have finished.');
disp('Click ENTER when the process is finished.');
pause();

%% QUALITY CONTROL
% Display of all subjects from the population for quality control in one
% image
file_to_control='Anat/mri/wmT1_3D.nii';
slice = 125; % Slice to display
cat12quality_control(Path,file_to_control,slice);

%% REALIGN and 4D to 3D
% 
rs_file = 'RS/RS_MB3_3x3x3_135mmFH_500dyn.nii';
dyn_num = 500;
realign_n_separate(Path,rs_file,dyn_num);

%% REGISTER atlas in rs native space
% 
rs_file = 'RS/meanRS_MB3_3x3x3_135mmFH_500dyn.nii';
atlas_file = '/usr/local/MATLAB/spm12/atlas/AAL3_mod_v1.nii';
atlas2rs(Path,rs_file,atlas_file);

%% DETECT ARTIFACTS
realigned_rs = 'RS/rRS_MB3_3x3x3_135mmFH_500dyn.nii';
motion_params = 'RS/rp_RS_MB3_3x3x3_135mmFH_500dyn.txt';
cfg_file = '/media/veronica/DATAPART2/EmoPark/Data/art_detect.cfg';
artifact_detect(Path,realigned_rs,motion_params,cfg_file);

%% REGISTER other file in rs native space 
% (created for perfusion maps but can be used for any nifti 
% 
rs_file = 'RS/meanRS_MB3_3x3x3_135mmFH_500dyn.nii';
file_to_coreg = 'pCASL/CBF.nii';
vox_dims = [3 3 3];
coreg2rs(Path,rs_file,file_to_coreg,vox_dims);

disp('Is everything ok ?.');
disp('Click ENTER when you want to continue graph analysis.');
pause();

% %% CREATE ATLAS COORD
% Atlas_path = '/usr/local/MATLAB/spm12/atlas/AAL3_mod_v1.nii';
% Save_coord_path = '/media/veronica/DATAPART2/EmoPark/Data/coord_AAL3_mod_v1.txt';
% I = niftiread(Atlas_path);
% stats = regionprops(I);
% centroid = reshape([stats.Centroid],3,112);
% centroid = centroid.';
% centroid(any(isnan(centroid), 2), :) = [];
% fid = fopen(Save_coord_path,'w'); 
% fprintf(fid,'%5.12f %5.12f %5.12f\n', centroid);
% fclose(fid);

%% TIME SERIES AND GRAPH COMPUTATION
Subj_dir = dir([Path '/*']);
Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));

for i = 1 : size(Subj_dir,1)
    disp(Subj_dir(i,1).name);
    
    % Compute time series ------------------------------------------------
    script_path = '/home/veronica/Code/rs_graph_processing';
    rs_path = fullfile(Path,Subj_dir(i,1).name, 'RS');
    r_dyn_name = '*rRS_MB3_3x3x3_135mmFH_500dyn_*';
    art_path = rs_path;
    atlas_file = fullfile(Path,Subj_dir(i,1).name, 'wAAL3_mod_v1.nii');
    gm_file = fullfile(Path,Subj_dir(i,1).name, 'rmwp1T1_3D.nii');
    graph_path = fullfile('/media/veronica/DATAPART2/EmoPark/Graphs/Controls_test',Subj_dir(i,1).name);
    %mkdir(graph_path)
    
    cmd_ts = [script_path '/compute_TS_graph_Matlab.R ' script_path ' ' rs_path ' ' r_dyn_name ' ' art_path ' ' atlas_file ' ' gm_file ' ' graph_path];
    system(cmd_ts) % Comment to skip
    
    % Compute correlations and graphs ------------------------------------
    file_coord = '/media/veronica/DATAPART2/EmoPark/Data/coord_AAL3_mod_v1.txt'; 
    graphs = 'TRUE';
    regions_selected = '"c(1:90, 97:106)"'; %'NULL' to choose all regions;
    percentage_selected = '3';
    num.levels = '3';
    num.nb.edges = '"seq(200,3000,200)"';
    length_time_series = 'NULL';
    
    cmd_graph = [script_path '/compute_TS_graph_Matlab.R ' script_path ' ' rs_path ' ' atlas_file ' ' graph_path ' ' file_coord ' ' graphs ' ' regions_selected ' ' percentage_selected ' ' num.levels ' ' num.nb.edges ' ' length_time_series];
    system(cmd_graph) % Comment to skip
end
  
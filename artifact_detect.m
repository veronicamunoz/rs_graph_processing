function artifact_detect(Path,file,motion_file,cfg_file)

Subj_dir = dir([Path '/*']);
Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));

for i = 1 : size(Subj_dir,1)
    disp(Subj_dir(i,1).name);
    
    rs_path=fullfile(Path, Subj_dir(i,1).name, file);
    motion_path=fullfile(Path, Subj_dir(i,1).name, motion_file);
    
    if (Subj_dir(i,1).isdir==1 && exist(rs_path, 'file') ~= 0 && exist(motion_path, 'file') ~= 0)
        cd(fullfile(Path, Subj_dir(i,1).name, 'RS'))
        art('sess_file', cfg_file)
        delete('*temp*')
    end  
   
    cd(Path)
end

% Subj_dir = dir([Path1 '/*']);
% Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));
% 
% i = 1;
% 
% disp(Subj_dir(i,1).name);
%     
% cd(fullfile(Path1, Subj_dir(i,1).name, 'RS'))
% rs_path=fullfile(Path1, Subj_dir(i,1).name, realigned_rs);
% motion_path=fullfile(Path1, Subj_dir(i,1).name, motion_params);
% liste_rs = cell(dyn_num,1);
% 
% if (Subj_dir(i,1).isdir==1 && exist(rs_path, 'file') ~= 0 && exist(motion_path, 'file') ~= 0)
%     for j = 1:dyn_num
%         liste_rs{j} = [rs_path ',' num2str(j)];
%     end
% end
% 
% clear B
% B.P = {rs_path};
% B.M = {motion_path};
% B.global_threshold = 3;
% B.motion_threshold = 3;
% B.use_diff_motion =1;    %: 1/0 use scan-to-scan differences in motion parameters
% B.use_diff_global=1;     %: 1/0 use scan-to-scan differences in global BOLD signal
% B.use_norms=0;           %: 1/0 use motion composite measure
% B.drop_flag=0;           %: number of initial scans to flag as outliers (removal of initial scans)
% B.motion_file_type=0;    %: indicates type of realignment file (0: SPM rp_*.txt file; 1: FSL .par file; 2: Siemens .txt file; 3: .txt SPM-format but rotation parameters in degrees)
% B.close = 1;             %: 1/0 close gui
% B.print =1;             %: 1/0 print gui
% B.output_dir =  fullfile(Path1, Subj_dir(i,1).name);         %: directory for output files (default same folder as first-session functional files)

end
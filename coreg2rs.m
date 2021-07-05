function coreg2rs(Path, rs_file, file_to_coreg, vox_dims)

Subj_dir = dir([Path '/*']);
Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));
C3Dcommand='/home/veronica/Downloads/Programs/c3d/bin/c3d';

for i = 1 : size(Subj_dir,1)
    disp(Subj_dir(i,1).name);
    
    rs_path=fullfile(Path, Subj_dir(i,1).name, rs_file);
    file_path = fullfile(Path, Subj_dir(i,1).name, file_to_coreg);
    
    if (Subj_dir(i,1).isdir==1 && exist(rs_path, 'file') ~= 0 && exist(file_path, 'file') ~= 0)
        cd(fullfile(Path, Subj_dir(i,1).name, 'RS'))
        
        if ~isempty(vox_dims) && exist(C3Dcommand) ~= 0
            system([C3Dcommand ' ' file_path ' -resample-mm ' num2str(vox_dims(1)) 'x' num2str(vox_dims(2)) 'x' num2str(vox_dims(3)) 'mm -o ' strrep(file_path,'.nii','_resampled.nii')]);

            if exist(strrep(file_path,'.nii','_resampled.nii'),'file') ~= 0
                file_path = strrep(file_path,'.nii','_resampled.nii');
            end
        end
        
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {rs_path};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {file_path};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end  
end

end
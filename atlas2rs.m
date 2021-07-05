function atlas2rs(Path, rs_file, atlas_path)

global spm_path
Subj_dir = dir([Path '/*']);
Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));

for i = 1 : size(Subj_dir,1)
    disp(Subj_dir(i,1).name);
    
    rs_path=fullfile(Path, Subj_dir(i,1).name, rs_file);
    
    if (Subj_dir(i,1).isdir==1 && exist(rs_path, 'file') ~= 0)
        cd(fullfile(Path, Subj_dir(i,1).name, 'RS'))
        
        clear matlabbatch
        spm_jobman('initcfg');

        matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {rs_path};
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {[spm_path '/tpm/TPM.nii']};
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;
        % The output from the normalize function is used as input for the
        % extraction of the inverse deformation field, this is called
        % dependency in SPM
        matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.def(1) = cfg_dep('Normalise: Estimate: Deformation (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
        matlabbatch{2}.spm.util.defs.comp{1}.inv.space = {rs_path};
        matlabbatch{2}.spm.util.defs.out{1}.savedef.ofname = 'inverse';
        matlabbatch{2}.spm.util.defs.out{1}.savedef.savedir.savepwd = 1;
        
        cd .. % The subject general folder is now the current folder 
        matlabbatch{3}.spm.util.defs.comp{1}.def(1) = cfg_dep('Deformations: Deformation', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def'));
        matlabbatch{3}.spm.util.defs.out{1}.pull.fnames = {atlas_path};
        matlabbatch{3}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
        matlabbatch{3}.spm.util.defs.out{1}.pull.interp = 0;
        matlabbatch{3}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{3}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{3}.spm.util.defs.out{1}.pull.prefix = 'w';

        % Deform gray matter map from MNI space to the rs native space
        gm_path=fullfile(Path, Subj_dir(i,1).name, 'Anat/mri/mwp1T1_3D.nii');
        if (exist(gm_path, 'file') ~= 0)
            matlabbatch{4}.spm.util.defs.comp{1}.def = cfg_dep('Deformations: Deformation', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def'));
            matlabbatch{4}.spm.util.defs.out{1}.pull.fnames = {gm_path};
            matlabbatch{4}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
            matlabbatch{4}.spm.util.defs.out{1}.pull.interp = 4;
            matlabbatch{4}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{4}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            matlabbatch{4}.spm.util.defs.out{1}.pull.prefix = 'r';
        end

        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end  
end

end
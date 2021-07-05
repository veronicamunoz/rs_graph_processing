function realign_n_separate(Path,file,dyns)

Subj_dir = dir([Path '/*']);
Subj_dir = Subj_dir(arrayfun(@(x) ~strcmp(x.name(1),'.'),Subj_dir));

for i = 1 : size(Subj_dir,1)
    disp(Subj_dir(i,1).name);
    
    % R E A L I G N
    
    cd(fullfile(Path, Subj_dir(i,1).name, 'RS'))
    rs_path=fullfile(Path, Subj_dir(i,1).name, file);
    liste_rs = cell(dyns,1);
    
    if (Subj_dir(i,1).isdir==1 && exist(rs_path, 'file') ~= 0)
        
        for j = 1:dyns
            liste_rs{j} = [fullfile(Path,Subj_dir(i,1).name,file) ',' num2str(j)];
        end
        
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {liste_rs};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 0.001;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end  
    
    % S E P A R A T E ( 4 D -> 3 D )
    file_parts = strsplit(file,'/');
    new_path = [Path Subj_dir(i,1).name '/' file_parts{1} '/r' file_parts{2}];
        
    if (Subj_dir(i,1).isdir==1 && exist(new_path, 'file') ~= 0)
        clear matlabbatch
        spm_jobman('initcfg');
        matlabbatch{1}.spm.util.split.vol = {new_path};
        matlabbatch{1}.spm.util.split.outdir = {''};
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
   
    cd(Path)
end


end
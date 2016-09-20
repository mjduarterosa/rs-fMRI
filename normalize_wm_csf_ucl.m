% -------------------------------------------------------------------------
% Normalize WM and CSF images (UCL)
% -------------------------------------------------------------------------
clear
close all

% Data folders
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');
rootdir = '/Users/maria/Documents/NSPN/code';

% Set FSL env
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE', 'NIFTI');

% Sites
site = {'ucl','cbsu','wbic'};
sitedir{1} = fullfile(fmridir,site{1});
sitedir{2} = fullfile(fmridir,site{2});
sitedir{3} = fullfile(fmridir,site{3});

% find directories of scans
for s=1:3
    
    dirs{s} = dir(sitedir{s});
    dirs{s}(1:2) = []; % removes . ..
    dirs{s}(~[dirs{s}.isdir]) = []; % removes stuff that is not a directory
    
    % write full directory path
    for j = 1:numel(dirs{s})
        fulldirs{s}{j} = fullfile(sitedir{s},dirs{s}(j).name);
    end
    
end

% Batch options
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

% UCL
ndir = 1;
for i = 1:numel(fulldirs{ndir})
    
    disp(sprintf('Scan %d out of %d', i, numel(fulldirs{ndir})));
    
    subdir = dir(fulldirs{ndir}{i});
    subdir(1:2) = []; % removes . ..
    subdir(~[subdir.isdir]) = [];
    
    % remove .svn stuff
    for j=1:numel(subdir)
        in(j) =  strcmp(subdir(j).name,'.svn');
    end
    subdir(find(in)) = [];
    
    [pathstr, name, ext] = fileparts(fulldirs{ndir}{i});
    
    t1_flag = 1;
    p = 1;
    ntim = 263;
    
    for s = 1:numel(subdir)
        
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
        
        if length(subdir(s).name) > 9 && t1_flag

            forwd = [subfulldirs,'/y_T1w.nii'];
            
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwd};
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[subfulldirs,'/c2T1w.nii']};
            matlabbatch{2}.spm.spatial.normalise.write.subj.def = {forwd};
            matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {[subfulldirs,'/c3T1w.nii']};
            
            t1_flag = 0;
            
        end
        
    end
    try
        spm('defaults', 'PET');
        spm_jobman('run', matlabbatch); 
    catch err
        error(err,'Error: data could not be preprocessed!!!');
    end  
end








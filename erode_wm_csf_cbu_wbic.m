% -------------------------------------------------------------------------
% Normalize WM and CSF images (UBC and WBIC)
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
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{1}.spm.util.imcalc.expression = 'i1 > 0.99';
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.util.imcalc.expression = 'i1 > 0.99';

% Preprocess CBU or WBIC
% % ndir = 2; % CBU
ndir = 3; % WBIC
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
    
    for s = 1:numel(subdir)
        
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
        [pathstr, name_sub, ext] = fileparts(subfulldirs);
        
        subsubdir = dir(subfulldirs);
        subsubdir(1:2) = []; % removes . ..
        subsubdir(~[subsubdir.isdir]) = [];
        
        for j = 1:numel(subsubdir)
            in(j) =  strcmp(subsubdir(j).name,'.svn');
        end
        subsubdir(find(in)) = [];
        
        t1_flag = 1;
        p = 1;
        ntim = 263;
        
        % Find data
        for q = 1:numel(subsubdir)
            strmepi = strfind(subsubdir(q).name,'_t1_'); % Find T1s
            if ~isempty(strmepi) && t1_flag
                
                matlabbatch{1}.spm.util.imcalc.input = {[fullfile(subfulldirs,subsubdir(q).name),'/wc2T1w.nii']};
                matlabbatch{2}.spm.util.imcalc.input = {[fullfile(subfulldirs,subsubdir(q).name),'/wc3T1w.nii']};
                matlabbatch{1}.spm.util.imcalc.output = 'mask_wc2T1w';
                matlabbatch{2}.spm.util.imcalc.output = 'mask_wc3T1w';
                matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(subfulldirs,subsubdir(q).name)};
                matlabbatch{2}.spm.util.imcalc.outdir = {fullfile(subfulldirs,subsubdir(q).name)};
               
                t1_flag = 0;
            end 
        end
        try
            spm('defaults', 'PET');
            spm_jobman('run', matlabbatch);
        catch err
            error(err,'Error: data could not be preprocessed!!!');
        end
        % Find data
        t1_flag = 1;
        for q = 1:numel(subsubdir)
            strmepi = strfind(subsubdir(q).name,'_t1_'); % Find T1s
            if ~isempty(strmepi) && t1_flag
                
                system(['/usr/local/fsl/bin/fslmaths ',fullfile(subfulldirs,subsubdir(q).name),'/mask_wc2T1w.nii -ero ' ,fullfile(subfulldirs,subsubdir(q).name),'/mask_wc2T1w.nii']);
                system(['/usr/local/fsl/bin/fslmaths ',fullfile(subfulldirs,subsubdir(q).name),'/mask_wc3T1w.nii -ero ' ,fullfile(subfulldirs,subsubdir(q).name),'/mask_wc3T1w.nii']);               
                t1_flag = 0;
            end 
        end
    end    
 end
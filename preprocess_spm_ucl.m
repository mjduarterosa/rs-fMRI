% -------------------------------------------------------------------------
% Main script for preprocessing fMRI and T1 data for UCL
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
matlabbatch{1}.spm.temporal.st.nslices = 34;
matlabbatch{1}.spm.temporal.st.tr = 2.42;
matlabbatch{1}.spm.temporal.st.ta = 2.34882352941176;
matlabbatch{1}.spm.temporal.st.so = [34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
matlabbatch{1}.spm.temporal.st.refslice = 17;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,1'};
matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,2'};
matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,3'};
matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,4'};
matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,5'};
matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {'/Users/maria/Documents/spm12/tpm/TPM.nii,6'};
matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];


matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [4 4 4];
matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';

matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w';

matlabbatch{7}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 0;
matlabbatch{7}.spm.spatial.smooth.prefix = 's';

matlabbatch{8}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{8}.spm.spatial.smooth.dtype = 0;
matlabbatch{8}.spm.spatial.smooth.im = 0;
matlabbatch{8}.spm.spatial.smooth.prefix = 's';

matlabbatch{9}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{9}.spm.spatial.smooth.dtype = 0;
matlabbatch{9}.spm.spatial.smooth.im = 0;
matlabbatch{9}.spm.spatial.smooth.prefix = 's';

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
            
            system(['gunzip ',subfulldirs,'/T1w.nii.gz']);
            
            matlabbatch{3}.spm.spatial.coreg.estimate.source = {[subfulldirs,'/T1w.nii']};
            matlabbatch{4}.spm.spatial.preproc.channel.vols = {[subfulldirs,'/T1w.nii']};
            
            forwd = [subfulldirs,'/y_T1w.nii'];
            
            matlabbatch{6}.spm.spatial.normalise.write.subj.def = {forwd};
            matlabbatch{6}.spm.spatial.normalise.write.subj.resample = {[subfulldirs,'/mT1w.nii']};
            
            t1_flag = 0;
            
        end
        
    end
    
    for s = 1:numel(subdir)
        
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
        
        if length(subdir(s).name) < 10
            
            system(['rm ',subfulldirs,'/BOLD.nii']);
            system(['rm ',subfulldirs,'/cBOLD.nii']);
            system(['rm ',subfulldirs,'/acBOLD.nii']);
            system(['rm ',subfulldirs,'/cBOLD.mat']);
            system(['rm ',subfulldirs,'/acBOLD.mat']);
            system(['rm ',subfulldirs,'/racBOLD.mat']);
            system(['rm ',subfulldirs,'/rp_cBOLD.txt']);
            system(['rm ',subfulldirs,'/rp_acBOLD.txt']);
            system(['rm ',subfulldirs,'/meanacBOLD.nii']);
            
            system(['/usr/local/fsl/bin/fslroi ',subfulldirs,'/BOLD.nii.gz ' ,subfulldirs,'/cBOLD.nii 6 263']);
            
            for sc = 1:ntim,
                matlabbatch{1}.spm.temporal.st.scans{p}{sc,1} = [subfulldirs,'/cBOLD.nii,',num2str(sc)];
                matlabbatch{2}.spm.spatial.realign.estwrite.data{p}{sc,1} = [subfulldirs,'/acBOLD.nii,',num2str(sc)];
                matlabbatch{5}.spm.spatial.normalise.write.subj(p).resample{sc,1} = [subfulldirs,'/racBOLD.nii,',num2str(sc)];
                matlabbatch{6+p}.spm.spatial.smooth.data{sc,1} = [subfulldirs,'/wracBOLD.nii,',num2str(sc)];
            end
            
            if p == 2
                matlabbatch{3}.spm.spatial.coreg.estimate.ref = {[subfulldirs,'/racBOLD.nii,1']};
            end
            
            matlabbatch{5}.spm.spatial.normalise.write.subj(p).def = {forwd};
            
            if p==3
                tmp = matlabbatch{2}.spm.spatial.realign.estwrite.data{1};
                matlabbatch{2}.spm.spatial.realign.estwrite.data{1} = matlabbatch{2}.spm.spatial.realign.estwrite.data{2};
                matlabbatch{2}.spm.spatial.realign.estwrite.data{2} = tmp;
            end
            
            p = p+1; 
            
        end
        
    end
    
    try
        spm('defaults', 'PET');
        spm_jobman('run', matlabbatch);
        
    catch err
        
        error(err,'Error: data could not be preprocessed!!!');
        
    end
    
    
end








% -------------------------------------------------------------------------
% Main script for extracting ROIs from sulci atlas
% -------------------------------------------------------------------------
clear
close all
clc

% ROIs
nregions = 1;

% Output file
f_out    = '/Users/maria/Documents/NSPN/analysis/ROI_mask_nspn_ucl.mat';

% Directories
dir_img  = '/Users/maria/Documents/mask/';
% mask_name = 'ENTER NAME HERE';
ex_scan  = spm_vol('/Users/maria/Documents/NSPN/data/fMRI/ucl/20130325.14506.MQ01614/MQ01614.3/swracBOLD.nii,1');

% Get filenames
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');

% Site
ntim = 263;
nTEs = 3;

% Sites
site = {'ucl'};
sitedir = fullfile(fmridir,site{1});

% find directories of scans
dirs = dir(sitedir);
dirs(1:2) = []; % removes . ..
dirs(~[dirs.isdir]) = []; % removes stuff that is not a directory

% write full directory path
for j = 1:numel(dirs)
    fulldirs{j} = fullfile(sitedir,dirs(j).name);
end

data = {};
% Run for UCL
disp('Finding data >>>>');
nsubs = numel(fulldirs);
for i = 1:nsubs
    subdir = dir(fulldirs{i});
    subdir(1:2) = []; % removes . ..
    subdir(~[subdir.isdir]) = [];
    % remove .svn stuff
    for j=1:numel(subdir)
        in(j) =  strcmp(subdir(j).name,'.svn');
    end
    subdir(find(in)) = [];
    [pathstr, name, ext] = fileparts(fulldirs{i});
    p = 1;
    for s = 1:numel(subdir) 
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{i},subdir(s).name);
        if length(subdir(s).name) < 10 && p<=nTEs
            for sc=1:ntim,
                data{i,p}{sc,1} = [subfulldirs,'/swracBOLD.nii,',num2str(sc)];
            end
            p = p+1;
        end
    end
end
disp('Finding data: done!');


% Get ROI regions
disp('Loading ROIs >>>>');
for r = 1:nregions
   
    finput   = [dir_img,mask_name,'.nii'];

    Vol = spm_vol(finput);
    
    % Image dimensions
    % ---------------------------------------------------------------------
    V             = Vol(1);
    M             = V.mat;
    DIM           = V.dim(1:3)';
    xdim          = DIM(1); ydim  = DIM(2); zdim  = DIM(3);
    [xords,yords] = ndgrid(1:xdim,1:ydim);
    xords         = xords(:)';  yords = yords(:)';
    I             = 1:xdim*ydim;
    zords_init    = ones(1,xdim*ydim);
    
    % Get image values above zero for each fold and all folds
    % ---------------------------------------------------------------------
    xyz_above = [];
    z_above   = [];
    
    for z = 1:zdim,
        zords = z*zords_init;
        xyz   = [xords(I); yords(I); zords(I)];
        nVox  = size(xyz,2);
        mask_xyz = Vol.mat\ex_scan.mat*[xyz(:,1:nVox);ones(1,nVox)];     
        zvals = spm_get_data(V,mask_xyz);
        above = find(~isnan(zvals) & zvals > 0);
        if ~isempty(above)
            xyz_above = [xyz_above,xyz(:,above)];
        end
    end
    XYZ{r}   = xyz_above(1:3,:);
    
end
disp('Loading ROIs: done.');
 
clear R
disp('Loading ROI data >>>>');
ROI = cell(nsubs,nTEs,nregions);
for r = 1:nregions
    fprintf('Region %d of %d------------->>\n',r,nregions);
    for s = 1:nsubs
        fprintf('Sub %d of %d\n',s,nsubs);
        for t=1:nTEs % Three TE time-series 
            p = data{s,t};
            Vol  = spm_vol(p);
            R    = spm_get_data(Vol,XYZ{r});
            ROI{s,t,r} = mean(R,2);
        end
    end
end
save(f_out,'ROI','-v7.3');
disp('Finished!');

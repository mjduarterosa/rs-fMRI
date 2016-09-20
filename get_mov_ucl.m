% -------------------------------------------------------------------------
% Main script for QCing movement
% -------------------------------------------------------------------------
clear
close all
clc

% Get filenames
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');

% Output directory
f_dir   = '/Users/maria/Documents/NSPN/analysis/';

% Site
ntim = 263;
nTE = 3;    % Change echo-time

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

% Run for UCL
ndir = 1;
nsubs = numel(fulldirs{ndir});
k = 1;
mov_param = zeros(ntim,6,nsubs);
clear fnames_ucl
for i = 1:nsubs
    subdir = dir(fulldirs{ndir}{i});
    subdir(1:2) = []; % removes . ..
    subdir(~[subdir.isdir]) = [];
    % remove .svn stuff
    for j=1:numel(subdir)
        in(j) =  strcmp(subdir(j).name,'.svn');
    end
    subdir(find(in)) = [];
    [pathstr, name, ext] = fileparts(fulldirs{ndir}{i});
    p = 1;
    for s = 1:numel(subdir)
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
        if length(subdir(s).name) < 10 && p<=nTEs
            if p == nTE
                cd(subfulldirs)
                [PATHSTR,NAME,EXT] = fileparts(subfulldirs);
                fnames_ucl(k:k+2,:) = repmat(NAME,3,1);
                k = k+3;
                if exist('rp_acBOLD.txt','file') == 2
                    mov_param(:,:,i) = load('rp_acBOLD.txt');
                end
                cd('/Users/maria/Documents/NSPN/code')
            end
            p = p+1;
        end
    end
end

save(fullfile(f_dir,sprintf('ROIs_move_nspn_ucl_%s.mat',num2str(nTE))),'mov_param'); 

clear mean_mov_param std_mov_param
mean_mov_param(:,:) = mean(mov_param);
std_mov_param(:,:) = std(mov_param);

% Make plots
y = reshape(repmat([1;2;3],1,size(mean_mov_param(1:3,:),2)),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xm1 = reshape(mean_mov_param(1:3,:),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xd1 = reshape(std_mov_param(1:3,:),1,size(std_mov_param(1:3,:),1)*size(std_mov_param(1:3,:),2));
xm2 = reshape(mean_mov_param(4:6,:),1,size(mean_mov_param(4:6,:),1)*size(mean_mov_param(4:6,:),2));
xd2 = reshape(std_mov_param(4:6,:),1,size(std_mov_param(4:6,:),1)*size(std_mov_param(4:6,:),2));
figure;subplot(2,2,1);boxplot(mean_mov_param(1:3,:)'); title('UCL');
for i=1:3*nsubs, text(y(i),xm1(i),fnames_ucl(i,:)); end
subplot(2,2,2);boxplot(mean_mov_param(4:6,:)'); title('UCL');
for i=1:3*nsubs, text(y(i),xm2(i),fnames_ucl(i,:)); end
subplot(2,2,3);boxplot(std_mov_param(1:3,:)'); title('UCL');
for i=1:3*nsubs, text(y(i),xd1(i),fnames_ucl(i,:)); end
subplot(2,2,4);boxplot(std_mov_param(4:6,:)'); title('UCL');
for i=1:3*nsubs, text(y(i),xd2(i),fnames_ucl(i,:)); end



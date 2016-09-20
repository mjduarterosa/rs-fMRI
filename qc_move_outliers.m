% -------------------------------------------------------------------------
% Main script for QCing movement
% -------------------------------------------------------------------------
clear
close all
clc

% Get filenames
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');

% Site
ntim = 263;
nTEs = 3;

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

% Run for CBU or WBIC
ndir = 3; % CBU
nsubs = numel(fulldirs{ndir});
k = 1;
mov_param = zeros(ntim,6,nsubs);
clear fnames_wbic
for i = 1:nsubs
    
    fprintf('Subject %d of %d------------->>\n',i,nsubs);
    
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
        
        p = 1;
        
        % Find data
        for q = 1:numel(subsubdir)
            
            strmepi = strfind(subsubdir(q).name,'_MEPI2_'); % Find MEPIs
            
            if ~isempty(strmepi) && p<=nTEs
                if p == 3
                    cd(fullfile(subfulldirs,subsubdir(q).name))
                    [PATHSTR,NAME,EXT] = fileparts(fullfile(subfulldirs,subsubdir(q).name));
                    fnames_wbic(k:k+2,:) = repmat(PATHSTR(end-23:end-19),3,1);
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
end

f_dir   = '/Users/maria/Documents/NSPN/analysis/';
save(fullfile(f_dir,'ROIs_move_nspn_wbic_3.mat'),'mov_param'); 

clear mean_mov_param std_mov_param
mean_mov_param(:,:) = mean(mov_param);
std_mov_param(:,:) = std(mov_param);

% Make plots
y = reshape(repmat([1;2;3],1,size(mean_mov_param(1:3,:),2)),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
l = reshape(repmat(1:nsubs,size(mean_mov_param(1:3,:),1),1),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xm1 = reshape(mean_mov_param(1:3,:),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xd1 = reshape(std_mov_param(1:3,:),1,size(std_mov_param(1:3,:),1)*size(std_mov_param(1:3,:),2));
xm2 = reshape(mean_mov_param(4:6,:),1,size(mean_mov_param(4:6,:),1)*size(mean_mov_param(4:6,:),2));
xd2 = reshape(std_mov_param(4:6,:),1,size(std_mov_param(4:6,:),1)*size(std_mov_param(4:6,:),2));
figure;subplot(2,2,1);boxplot(mean_mov_param(1:3,:)'); title('WBIC');
for i=1:3*nsubs, text(y(i),xm1(i),fnames_wbic(i,:)); end
subplot(2,2,2);boxplot(mean_mov_param(4:6,:)'); title('WBIC');
for i=1:3*nsubs, text(y(i),xm2(i),fnames_wbic(i,:)); end
subplot(2,2,3);boxplot(std_mov_param(1:3,:)'); title('WBIC');
for i=1:3*nsubs, text(y(i),xd1(i),fnames_wbic(i,:)); end
subplot(2,2,4);boxplot(std_mov_param(4:6,:)'); title('WBIC');
for i=1:3*nsubs, text(y(i),xd2(i),fnames_wbic(i,:)); end

% Run for CBU or WBIC
ndir = 2; % CBU
nsubs = numel(fulldirs{ndir});
k = 1;
mov_param = zeros(ntim,6,nsubs);
clear fnames_cbu
for i = 1:nsubs
    
    fprintf('Subject %d of %d------------->>\n',i,nsubs);
    
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
        
        p = 1;
        
        % Find data
        for q = 1:numel(subsubdir)
            
            strmepi = strfind(subsubdir(q).name,'_MEPI2_'); % Find MEPIs
            
            if ~isempty(strmepi) && p<=nTEs 
                if p == 3
                    cd(fullfile(subfulldirs,subsubdir(q).name))
                    [PATHSTR,NAME,EXT] = fileparts(fullfile(subfulldirs,subsubdir(q).name));
                    fnames_cbu(k:k+2,:) = repmat(PATHSTR(end-37:end-32),3,1);
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
end

f_dir   = '/Users/maria/Documents/NSPN/analysis/';
save(fullfile(f_dir,'ROIs_move_nspn_cbu_3.mat'),'mov_param'); 

clear mean_mov_param std_mov_param
mean_mov_param(:,:) = mean(mov_param);
std_mov_param(:,:) = std(mov_param);

% Make plots
y = reshape(repmat([1;2;3],1,size(mean_mov_param(1:3,:),2)),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
l = reshape(repmat(1:nsubs,size(mean_mov_param(1:3,:),1),1),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xm1 = reshape(mean_mov_param(1:3,:),1,size(mean_mov_param(1:3,:),1)*size(mean_mov_param(1:3,:),2));
xd1 = reshape(std_mov_param(1:3,:),1,size(std_mov_param(1:3,:),1)*size(std_mov_param(1:3,:),2));
xm2 = reshape(mean_mov_param(4:6,:),1,size(mean_mov_param(4:6,:),1)*size(mean_mov_param(4:6,:),2));
xd2 = reshape(std_mov_param(4:6,:),1,size(std_mov_param(4:6,:),1)*size(std_mov_param(4:6,:),2));
figure;subplot(2,2,1);boxplot(mean_mov_param(1:3,:)'); title('CBU');
for i=1:3*nsubs, text(y(i),xm1(i),fnames_cbu(i,:)); end
subplot(2,2,2);boxplot(mean_mov_param(4:6,:)'); title('CBU');
for i=1:3*nsubs, text(y(i),xm2(i),fnames_cbu(i,:)); end
subplot(2,2,3);boxplot(std_mov_param(1:3,:)'); title('CBU');
for i=1:3*nsubs, text(y(i),xd1(i),fnames_cbu(i,:)); end
subplot(2,2,4);boxplot(std_mov_param(4:6,:)'); title('CBU');
for i=1:3*nsubs, text(y(i),xd2(i),fnames_cbu(i,:)); end

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
            if p == 3
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

f_dir   = '/Users/maria/Documents/NSPN/analysis/';
save(fullfile(f_dir,'ROIs_move_nspn_ucl_3.mat'),'mov_param'); 

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



 clear 
close all

% New folders for data copies
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');
rootdir = '/Users/maria/Documents/NSPN/code';

% Set FSL env
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');

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
    
    % only keep directories with defined length (to remove the rest)
    for j = 1:numel(dirs{s})
        len{s}(j) = length(dirs{s}(j).name);
    end
    
    % find unique lenghts and frequencies
    u = unique(len{s});
    for j = 1:length(u)
        fq(j) = sum(u(j)==len{s});
    end
    
    % remove directories with smallest frequency
    [ma,in] = max(fq);
    dirs{s}(~[u(in)==len{s}])=[];
    
    % write full directory path
    for j = 1:numel(dirs{s})
        fulldirs{s}{j} = fullfile(sitedir{s},dirs{s}(j).name);
    end
    
end
% 
% % WBIC
% ndir = 2;
% for i = 1:numel(fulldirs{ndir})
% 
%     disp(sprintf('Scan %d out of %d', i, numel(fulldirs{ndir})));
%     
%     subdir = dir(fulldirs{ndir}{i});
%     subdir(1:2) = []; % removes . ..
%     subdir(~[subdir.isdir]) = [];
%     
%     % remove .svn stuff
%     for j=1:numel(subdir)
%         in(j) =  strcmp(subdir(j).name,'.svn');
%     end
%     subdir(find(in)) = [];
%     
%     [pathstr, name, ext] = fileparts(fulldirs{ndir}{i});
%     
%     for s = 1:numel(subdir)
%         
%         % only the first subdirectory
%         subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
%         [pathstr, name_sub, ext] = fileparts(subfulldirs);
%         
%         subsubdir = dir(subfulldirs);
%         subsubdir(1:2) = []; % removes . ..
%         subsubdir(~[subsubdir.isdir]) = [];
%         
%         for j = 1:numel(subsubdir)
%             in(j) =  strcmp(subsubdir(j).name,'.svn');
%         end
%         subsubdir(find(in)) = [];
%         
%         % Find MEPI data
%         % Copy *.dcm to local directory
%         % Move to new directory
%         for q = 1:numel(subsubdir)
%             
%             strmepi = strfind(subsubdir(q).name,'_MEPI2_'); % Find MEPI2
%             
%             if ~isempty(strmepi)
%                 
% %                 images = dir(fullfile(subfulldirs,subsubdir(q).name));
% %                 images(1:2) = [];
%                 cd(fullfile(subfulldirs,subsubdir(q).name))
%                 
% 
%                 try
% %                     system('/usr/local/fsl/bin/fslmerge -tr BOLD *.nii 2.429719');
%                     
%                     system('rm *.nii');
%                     cd(rootdir)
%                     
% 
%                 catch err
%                     
%                     error(err,'Error: data could not be averaged!!!');
%                     
%                 end
%                 
%             end
%            
%         end
%         
%     end
%     
% end

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
    
    for s = 1:numel(subdir)
        
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{i},subdir(s).name);
        
        if length(subdir(s).name) < 10
            
            cd(fullfile(subfulldirs))
            
            
            try
                system('/usr/local/fsl/bin/fslmerge -tr BOLD *.nii 2.429719');
                
                system('rm *.nii');
                cd(rootdir)
                
                
            catch err
                
                error(err,'Error: data could not be averaged!!!');
                
            end
            
        end
        
    end
    
end



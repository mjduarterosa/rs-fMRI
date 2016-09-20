clear 
close all

% First used wget --spider -r --level=inf url to create directory structure locally

% Directory sctructure saved here
rawdir = '/Users/maria/Documents/NSPN/omd.fil.ion.ucl.ac.uk';
url = 'http://omd.fil.ion.ucl.ac.uk/';

% New folders for data copies
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');

% Sites
site = {'ucl','cbsu','wbic'};
sitedir{1} = fullfile(rawdir,site{1},'uchange');
sitedir{2} = fullfile(rawdir,site{2},'uchange','mri');
sitedir{3} = fullfile(rawdir,site{3},'uchange','dicom');
urldir{1} = fullfile(url,site{1},'uchange');
urldir{2} = fullfile(url,site{2},'uchange','mri');
urldir{3} = fullfile(url,site{3},'uchange','dicom');

% Create site directories
% for s=1:3
%      mkdir(fullfile(fmridir,site{s}));
% end

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
        fulldirsurl{s}{j} = fullfile(urldir{s},dirs{s}(j).name);
    end
    
end

% Find CBU scan lists
cbu_scans = [];
for j = 1:numel(dirs{2})
    cbu_scans = [cbu_scans; dirs{2}(j).name(1:9)];
end

% Find UCL scan lists
ucl_scans = [];
for j = 1:numel(dirs{1})
    if (j == 100) ||  (j == 98)
        ucl_scans = [ucl_scans; dirs{1}(j).name(end-12:end-6)];
    else
        ucl_scans= [ucl_scans; dirs{1}(j).name(end-6:end)];
    end
end

% Find WBIC scan lists
wbic_scans = [];
for j = 1:numel(dirs{3})
    wbic_scans = [wbic_scans; dirs{3}(j).name];
end

% Read all baseline files
NSPNscans = csvread('nspn_subjects_baseline.csv');

% Get all baseline scans in NSPN database
cbu_scans_nspn = [];
ucl_scans_nspn = [];
wbic_scans_nspn = [];
for i = 1:size(NSPNscans,1)
   if NSPNscans(i,2) == 1
       cbu_scans_nspn = [cbu_scans_nspn; char(['CBU',num2str(NSPNscans(i,3))])];
   end
   if NSPNscans(i,2) == 2
       ucl_scans_nspn = [ucl_scans_nspn; char(['MQ0',num2str(NSPNscans(i,3))])];
   end
   if NSPNscans(i,2) == 3
       wbic_scans_nspn = [wbic_scans_nspn; num2str(NSPNscans(i,3))];
   end 
end

% Find ones available for download
fid=fopen('NSPN_scans_all.csv','wt');

cbu_scans_in = [];
flag = 0;
for i = 1:size(cbu_scans,1)
    if ~isempty(find((str2num(cbu_scans_nspn(:,4:end)) - repmat(str2num(cbu_scans(i,4:end)),size(cbu_scans_nspn,1),1))==0))
        if strcmp(cbu_scans(i,4:end),'140313')
            if flag == 0
                cbu_scans_in = [cbu_scans_in; char(cbu_scans(i,:))];
                fprintf(fid,[cbu_scans(i,4:end),'\n']);
%                 cbu_scans_in = [cbu_scans_in; i];
            end
            flag = flag+1;
        else
            cbu_scans_in = [cbu_scans_in; char(cbu_scans(i,:))];
            fprintf(fid,[cbu_scans(i,4:end),'\n']);
%             cbu_scans_in = [cbu_scans_in; i];
        end
    end
    
end

wbic_scans_in = [];
for i = 1:size(wbic_scans,1)
    if i ~= 204
        if ~isempty(find((str2num(wbic_scans_nspn(:,:)) - repmat(str2num(wbic_scans(i,:)),size(wbic_scans_nspn,1),1))==0))
            wbic_scans_in = [wbic_scans_in; wbic_scans(i,:)];
            fprintf(fid,[wbic_scans(i,:),'\n']);
            %         wbic_scans_in = [wbic_scans_in; i];
            
        end
    end
end


ucl_scans_in = [];
flag = 0;
for i = 1:size(ucl_scans,1)
    if i ~= 43
        if ~isempty(find((str2num(ucl_scans_nspn(:,4:end)) - repmat(str2num(ucl_scans(i,4:end)),size(ucl_scans_nspn,1),1))==0))
            if ((ucl_scans(i,3) == num2str(0)) && (ucl_scans(i,2) == 'Q'))
                if strcmp(ucl_scans(i,4:end),'1947')
                    if flag == 0
                        ucl_scans_in = [ucl_scans_in; ucl_scans(i,:)];
                        fprintf(fid,[ucl_scans(i,3:end),'\n']);
                        %                     ucl_scans_in = [ucl_scans_in; i];
                    end
                    flag = flag+1;
                else
                    ucl_scans_in = [ucl_scans_in; ucl_scans(i,:)];
                    fprintf(fid,[ucl_scans(i,3:end),'\n']);
                    %                     ucl_scans_in = [ucl_scans_in; i];
                end
            end
        end
    end
    
end
fclose(fid)





keyboard

ndir = 2; % import cbsu
fid = fopen('CBU_scans_err.txt', 'a');
for i = 1:size(cbu_scans_in,1)
    subdir = dir(fulldirs{ndir}{cbu_scans_in(i)});
    subdir(1:2) = []; % removes . ..
    subdir(~[subdir.isdir]) = [];
    
    % remove .svn stuff
    for j=1:numel(subdir)
        in(j) =  strcmp(subdir(j).name,'.svn');
    end
    subdir(find(in)) = [];
    
    [pathstr, name, ext] = fileparts(fulldirs{ndir}{cbu_scans_in(i)});
    
    for s = 1:numel(subdir)
        
        if numel(subdir) > 2
            fprintf(fid, ['\n',num2str(cbu_scans_in(i)),'\n']);
        end
        
        % only the first subdirectory
        subfulldirs  = fullfile(fulldirs{ndir}{cbu_scans_in(i)},subdir(s).name);
        subfulldirsurl = fullfile(fulldirsurl{ndir}{cbu_scans_in(i)},subdir(s).name);
        [pathstr, name_sub, ext] = fileparts(subfulldirs);
        
        subsubdir = dir(subfulldirs);
        subsubdir(1:2) = []; % removes . ..
        subsubdir(~[subsubdir.isdir]) = [];
        
        for j = 1:numel(subsubdir)
            in(j) =  strcmp(subsubdir(j).name,'.svn');
        end
        subsubdir(find(in)) = [];
        
        % Find MEPI data
        % Copy *.dcm to local directory
        % Move to new directory
        for q = 1:numel(subsubdir)
            
            strmepi = strfind(subsubdir(q).name,'_t1_'); % Find T1s
            %strmepi = strfind(subsubdir(q).name,'_MEPI2_'); % Find MEPI2
            
            if ~isempty(strmepi)
                
                copydir = fullfile(subfulldirsurl,subsubdir(q).name);
                new_dir = fullfile(fmridir,site{ndir},name,name_sub,subsubdir(q).name);

                mkdir(new_dir);
                
                try
                    
                    system(['/usr/local/bin/wget -A dcm -r -l 1 -nd ', copydir]);
                    system(['mv *.dcm ',new_dir])

                catch err
                    
                    error(err,'Error: data could not be copied!!!');
                    
                end
                
                dicomfiles = spm_select('FPList',new_dir,'.*');
                
                try
                    
                    matlabbatch{1}.spm.util.import.dicom.data = cellstr(dicomfiles);
                    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
                    matlabbatch{1}.spm.util.import.dicom.outdir = {new_dir};
                    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
                    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
                    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
                    spm('defaults', 'PET');
                    spm_jobman('run', matlabbatch);
                    
                    system(['find ', new_dir, ' -name ''*.dcm'' | xargs rm'])
                    system(['rm ',new_dir,'/*.dcm']);
                    
                catch err
                    
                    error(err,'Error: data could not be converted!!!');
                    
                    
                end
                
            end
           
        end
        
    end
    
end

fclose(fid);

% ndir = 1; % import ucl
% for i = 1:size(ucl_scans_in,1)
%    
%    % Create new directory
%    [pathstr, name, ext] = fileparts(fulldirs{ndir}{ucl_scans_in(i)});
%    
%    try
%        
% %        system(['/usr/local/bin/wget -A tar -r -l 1 -nd ', fulldirsurl{ndir}{ucl_scans_in(i)}]);
%         system(['/usr/local/bin/wget -A 12.tar -r -l 1 -nd ', fulldirsurl{ndir}{ucl_scans_in(i)}]);
%         system(['/usr/local/bin/wget -A 13.tar -r -l 1 -nd ', fulldirsurl{ndir}{ucl_scans_in(i)}]);
%        
% %        tarfiles = spm_select('FPList',pwd,'.tar');
%        tarfiles = [spm_select('FPList',pwd,'.12.tar'); spm_select('FPList',pwd,'.13.tar')];
%        
%        for j = 1:size(tarfiles,1)
%            
%            [tmp,dname,tmp]=fileparts(tarfiles(j,:));
%            untar(strtrim(tarfiles(j,:)),dname);
%            imafiles = spm_select('FPList',fullfile(pwd,dname),'.ima');
%            
%            matlabbatch{1}.spm.util.import.dicom.data = cellstr(imafiles);
%            matlabbatch{1}.spm.util.import.dicom.root = 'flat';
%            matlabbatch{1}.spm.util.import.dicom.outdir = {fullfile(pwd,dname)};
%            matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
%            matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
%            matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
%            spm('defaults', 'PET');
%            spm_jobman('run', matlabbatch);
%            
%            niifiles = spm_select('FPList',fullfile(pwd,dname),'.nii');
%            
%            %            if strfind(niifiles(1,:),'fMQ')
%            
%            [tmp,dirname,dirname2] = fileparts(fulldirs{ndir}{ucl_scans_in(i)});
%            new_dir = fullfile(fmridir,site{ndir},[dirname,dirname2],dname);
%            
%            mkdir(new_dir)
%            
%            system(['rm ', dname, '/*.ima']);
%            system(['cp -r ', dname, '/ ', new_dir]);
%            system(['rm -r ',dname]);
%            
% %            else
% %                
% %                system(['rm -r ',dname]);
% %                
% %            end 
%            
%        end
%        
%        system('rm *.tar');
%        
%    catch err
%        
%        error(err,'Error: data could not be copied!!!');
%        
%    end
% 
% end

% ndir = 3; % import wbic
% for i = 1:size(wbic_scans_in,1)
% 
%     subdir = dir(fulldirs{ndir}{wbic_scans_in(i)});
%     subdir(1:2) = []; % removes . ..
%     subdir(~[subdir.isdir]) = [];
%     
%     % remove .svn stuff
%     for j=1:numel(subdir)
%         in(j) =  strcmp(subdir(j).name,'.svn');
%     end
%     subdir(find(in)) = [];
%     
%     % Create new directory
%     [pathstr, name, ext] = fileparts(fulldirs{ndir}{wbic_scans_in(i)});
%     
% %     for s = 1:numel(subdir)
%     for s = 1:1 % Get only first scan (baseline data)
%         
%         % only the first subdirectory
%         subfulldirs  = fullfile(fulldirs{ndir}{wbic_scans_in(i)},subdir(s).name);
%         subfulldirsurl = fullfile(fulldirsurl{ndir}{wbic_scans_in(i)},subdir(s).name);
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
%             if ~isempty(strfind(subsubdir(q).name,'_MEPI2_')) || ~isempty(strfind(subsubdir(q).name,'_t1_'))
%                 copydir = fullfile(subfulldirsurl,subsubdir(q).name);
%                 new_dir = fullfile(fmridir,site{ndir},name,name_sub,subsubdir(q).name);
%                 mkdir(new_dir);
%                 
%                 try
%                     
%                     system(['/usr/local/bin/wget -A dcm -r -l 1 -nd ', copydir]);
%                     system(['mv *.dcm ',new_dir])
%                     
%                 catch err
%                     
%                     error(err,'Error: data could not be copied!!!');
%                     
%                 end
%                 
%                 dicomfiles = spm_select('FPList',new_dir,'.*');
%                 
%                 try
%                     
%                     matlabbatch{1}.spm.util.import.dicom.data = cellstr(dicomfiles);
%                     matlabbatch{1}.spm.util.import.dicom.root = 'flat';
%                     matlabbatch{1}.spm.util.import.dicom.outdir = {new_dir};
%                     matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
%                     matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
%                     matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
%                     spm('defaults', 'PET');
%                     spm_jobman('run', matlabbatch);
%                     
%                     system(['find ', new_dir, ' -name ''*.dcm'' | xargs rm'])
%                     
%                 catch err
%                     
%                     error(err,'Error: data could not be converted!!!');
%                     
%                     
%                 end
%                 
%             end
%         end
%         
%     end
% 
% end
% 
% 

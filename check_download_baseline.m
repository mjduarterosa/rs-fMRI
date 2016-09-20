clear 
close all

% New folders for data copies
fmridir = fullfile('/Users/maria/Documents/NSPN/data/fMRI');

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
ucl_scans_in = [];
flag = 0;
for i = 1:size(ucl_scans,1)
    if ~isempty(find((str2num(ucl_scans_nspn(:,4:end)) - repmat(str2num(ucl_scans(i,4:end)),size(ucl_scans_nspn,1),1))==0))
        if ((ucl_scans(i,3) == num2str(0)) && (ucl_scans(i,2) == 'Q'))
            if strcmp(ucl_scans(i,4:end),'1947')
                if flag == 0
                    % ucl_scans_in = [ucl_scans_in; ucl_scans(i,:)];
                    ucl_scans_in = [ucl_scans_in; i];
                end
                flag = flag+1;
            else
                    % ucl_scans_in = [ucl_scans_in; ucl_scans(i,:)];
                    ucl_scans_in = [ucl_scans_in; i];
            end
        end
    end
    
end

cbu_scans_in = [];
flag = 0;
for i = 1:size(cbu_scans,1)
    if ~isempty(find((str2num(cbu_scans_nspn(:,4:end)) - repmat(str2num(cbu_scans(i,4:end)),size(cbu_scans_nspn,1),1))==0))
        if strcmp(cbu_scans(i,4:end),'140313')
            if flag == 0
                % cbu_scans_in = [cbu_scans_in; cbu_scans(i,:)];
                cbu_scans_in = [cbu_scans_in; i];
            end
            flag = flag+1;
        else
            % cbu_scans_in = [cbu_scans_in; cbu_scans(i,:)];
            cbu_scans_in = [cbu_scans_in; i];
        end
    end
    
end

wbic_scans_in = [];
for i = 1:size(wbic_scans,1)
    if ~isempty(find((str2num(wbic_scans_nspn(:,:)) - repmat(str2num(wbic_scans(i,:)),size(wbic_scans_nspn,1),1))==0))
        % wbic_scans_in = [wbic_scans_in; wbic_scans(i,:)];
        wbic_scans_in = [wbic_scans_in; i];

    end
    
end


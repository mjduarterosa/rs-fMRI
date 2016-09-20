% Create partial correlation matrices using FSL netmats code

% Load clean data
files{1} = 'ROIs_clean_nspn_cbu.mat';
files{2} = 'ROIs_clean_nspn_wbic.mat';
files{3} = 'ROIs_clean_nspn_ucl.mat';

fdir = '/Users/maria/Documents/NSPN/analysis';

nTE = 2;
nreg = 137;

% NET (all netmats)
NET = [];

% Calculate netmats for CBU
load(fullfile(fdir,files{1}))
netmats_cbu = [];
nsub = size(ROIcbu,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    t = [];
    for j = 1:nreg
        t(:,j) = ROIcbu{s,nTE,j};
    end
    t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
    tmp = nets_netmats(t,0,'corr');
    netmats_cbu(s,:,:) = tmp;
    ri = tril(ones(size(tmp,1), size(tmp,1)),-1);
    r = tmp(ri==1)';
    NET = [NET; r];
end

% Calculate netmats for WBIC
load(fullfile(fdir,files{2}))
netmats_wbic = [];
nsub = size(ROIwbic,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    t = [];
    for j = 1:nreg
        t(:,j) = ROIwbic{s,nTE,j};
    end
    t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
    tmp = nets_netmats(t,0,'corr');
    netmats_wbic(s,:,:) = tmp;
    ri = tril(ones(size(tmp,1), size(tmp,1)),-1);
    r = tmp(ri==1)';
    NET = [NET; r];
end

% Calculate netmats for UCL
load(fullfile(fdir,files{3}))
netmats_ucl = [];
nsub = size(ROIucl,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    t = [];
    for j = 1:nreg
        t(:,j) = ROIucl{s,nTE,j};
    end
    t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
    tmp = nets_netmats(t,0,'corr');
    netmats_ucl(s,:,:) = tmp;
    ri = tril(ones(size(tmp,1), size(tmp,1)),-1);
    r = tmp(ri==1)';
    NET = [NET; r];
end

save(fullfile(fdir,'netmats_NET_corr_TE30.mat'),'NET','-v7.3');




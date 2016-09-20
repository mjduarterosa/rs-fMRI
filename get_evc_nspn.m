% Create partial correlation matrices using FSL netmats code

% Load clean data
files{1} = 'ROIs_clean_nspn_cbu.mat';
files{2} = 'ROIs_clean_nspn_wbic.mat';
files{3} = 'ROIs_clean_nspn_ucl.mat';

fdir = '/Users/maria/Documents/NSPN/analysis';

nTEs = 3;
nreg = 137;

% NET (all netmats)
NET = [];

% Calculate netmats for CBU
load(fullfile(fdir,files{1}))
netmats_cbu = [];
nsub = size(ROIcbu,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    netmats = zeros(nreg,nreg);
    for i = 1:nTEs
        t = [];
        for j = 1:nreg
            t(:,j) = ROIcbu{s,i,j};
        end
        t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
        netmats = netmats + nets_netmats(t,0,'corr');
    end
    tmp = abs(netmats/nTEs);
    r = eigenvector_centrality_und(tmp)';
    NET = [NET; r]; 
end

% Calculate netmats for WBIC
load(fullfile(fdir,files{2}))
netmats_wbic = [];
nsub = size(ROIwbic,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    netmats = zeros(nreg,nreg);
    for i = 1:nTEs
        t = [];
        for j = 1:nreg
            t(:,j) = ROIwbic{s,i,j};
        end
        t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
        netmats = netmats + nets_netmats(t,0,'corr');
    end
    tmp = abs(netmats/nTEs);
    r = eigenvector_centrality_und(tmp)';
    NET = [NET; r]; 
end

% Calculate netmats for UCL
load(fullfile(fdir,files{3}))
netmats_ucl = [];
nsub = size(ROIucl,1);
for s = 1: nsub,
    fprintf('\nSubject %d out of %d >>>>\n',s,nsub);
    netmats = zeros(nreg,nreg);
    for i = 1:nTEs
        t = [];
        for j = 1:nreg
            t(:,j) = ROIucl{s,i,j};
        end
        t = (t-repmat(mean(t),size(t,1),1))./repmat(std(t),size(t,1),1);
        netmats = netmats + nets_netmats(t,0,'corr');
    end
    tmp = abs(netmats/nTEs);
    r = eigenvector_centrality_und(tmp)';
    NET = [NET; r]; 
end

save(fullfile(fdir,'netmats_EVC.mat'),'NET','-v7.3');




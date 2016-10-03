% -------------------------------------------------------------------------
% Permute CCA - NSPN data (all data - healthy and depressed)
%
% Based on Smith et al. Nature Neuroscience 2016
% Additional matlab toolboxes required (will need to addpath for each of these)
% FSLNets     http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
% PALM        http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
% nearestSPD  http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
% -------------------------------------------------------------------------

% Clear
% -------------------------------------------------------------------------
close all
clear all
clc

% Options
% -------------------------------------------------------------------------
cohort = 'healthy';
% cohort = 'all';
inc_age_gender = 0; % include age and gender as confound?

% Load data
% -------------------------------------------------------------------------
load_data_cca_nspn;

% Process confounds
% -------------------------------------------------------------------------

% Choose which confounds to use
if inc_age_gender
    varsQconf = [varsQconf, fd, age]; % remove age
else
    varsQconf = [varsQconf, fd]; % remove age
end
varsQconftmp = varsQconf;

% Setup confounds matrix
varsQconftmp(varsQconftmp==999)=0;  % impute missing data as zeros
conf = palm_inormal(varsQconftmp);  % Gaussianise
conf(isnan(conf))=0;                % impute missing data as zeros
conf = nets_normalise(conf);        % normalise

% Choose other confounds to use
if inc_age_gender
    conf = [conf mri_centre gender];
else
    conf = [conf mri_centre];
end

% Process cognitive/clinical variables
% -------------------------------------------------------------------------

% Treat variables (missing data)
if strcmp(cohort,'healthy')
    nn = 269;
else
    nn = 312;
end
vars(vars==999) = NaN;
vars(:,sum(~isnan(vars))<nn) = [];  % pre-delete vars with LOADS of missing data
vars = vars(:,2:end);
med = repmat(median(vars,'omitnan'),size(vars,1),1); % missing data imputation using the median
vars(isnan(vars)) = med(isnan(vars));

% Remove WASI (intellegence scores)
% vars = vars(:,1:end-9);
N = size(vars,1);

% prepare main netmat matrix
NET1 = nets_demean(NET);  NET1=NET1/std(NET1(:));   % no norm
grot = NET1;
NETd = nets_demean(grot-conf*(pinv(conf)*grot));    % deconfound and demean

% Prepare permutations
% -------------------------------------------------------------------------

% prepare permutation scheme using PALM - for more details see:
Nperm = 10000;                                                                       
EB = ones(N,1);
PAPset = palm_quickperms([ ], EB, Nperm);

% Prepare cognitive/clinical variables
% -------------------------------------------------------------------------

% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values
badvars=[];
nkeep = 250; % Smith default
for i=1:size(vars,2)
  Y=vars(:,i); grotKEEP=~isnan(Y);  
  grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?
  if (sum(grotKEEP)>nkeep) & (std(Y(grotKEEP))>0) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95) & (grot<100)
  else
%     [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
    badvars=[badvars i];
  end
end

% get list of which SMs to feed into CCA
varskeep=setdiff([1:size(vars,2)],badvars);                                                                

% "impute" missing vars data - actually this avoids any imputation
varsd=palm_inormal(vars(:,varskeep)); % Gaussianise
for i=1:size(varsd,2) % deconfound ignoring missing data
   grot=(isnan(varsd(:,i))==0);  
   grotconf=nets_demean(conf(grot,:)); 
   varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsdCOV=zeros(size(varsd,1));
for i=1:size(varsd,1) % estimate "pairwise" covariance, ignoring missing data
  for j=1:size(varsd,1)
    grot=varsd([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); varsdCOV(i,j)=grot(1,2);
  end
end
varsdCOV2=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix


% Get different PCA datasets
% -------------------------------------------------------------------------
pca_comp = [10, 25, 50, 75, 100, 125, 150, 200, 225, 250, 275, 300];
for Nkeep = 1:length(pca_comp),
    
    disp(sprintf('PCA dataset %d / Number of components %d >>>>>', Nkeep, pca_comp(Nkeep)));
    
    [uu1,ss1,vv1]=nets_svds(NETd,pca_comp(Nkeep));     % SVD reduction
    [uu,dd]=eigs(varsdCOV2,pca_comp(Nkeep));           % SVD (eigs actually)
    uu2=uu-conf*(pinv(conf)*uu);             % deconfound again just to be safe
    
    % Run CCA
    % -------------------------------------------------------------------------
    [grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,uu2);
    
    % Run CCA - permutation tests
    % -------------------------------------------------------------------------
    grotRp=zeros(Nperm,pca_comp(Nkeep)); clear grotRpval;
    for j=1:Nperm
        [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,uu2(PAPset(:,j),:));
    end
    for i=1:pca_comp(Nkeep);  % get FWE-corrected pvalues
        grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
        %         grotRpval(i) = (grotR(1)-sum(grotRp(2:end,1))/Nperm)/std(grotRp(2:end,1));
    end
    
    Ncca=sum(grotRpval<0.05);  % number of FWE-significant CCA components
    
    disp(sprintf('Number of FWE-significant CCA components: %d / Canonical Correlation: %f / P-value: %f',Ncca,grotR(1),grotRpval(1)));
    
    perm_pvals(Nkeep) = grotRpval(1);
    
end

% Get best PCA
% -------------------------------------------------------------------------
[m,mi]=min(perm_pvals);
[uu1,ss1,vv1]=nets_svds(NETd,pca_comp(mi));     % SVD reduction
[uu,dd]=eigs(varsdCOV2,pca_comp(mi));           % SVD (eigs actually)
uu2=uu-conf*(pinv(conf)*uu);             % deconfound again just to be safe

% Get CCA modes
% -------------------------------------------------------------------------
% Re-run CCA
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,uu2);

% netmat weights for CCA mode 1
grotAA = corr(grotU(:,1),NET)';
% or
grotAAd = corr(grotU(:,1),NETd(:,1:size(NET,2)))'; % weights after deconfounding

% variable weights for CCA mode 1
grotBB = corr(grotV(:,1),palm_inormal(vars),'rows','pairwise');
% or 
varsgrot=palm_inormal(vars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end

grotBBd = corr(grotV(:,1),varsgrot,'rows','pairwise')'; % weights after deconfounding

% Plot results
% -------------------------------------------------------------------------
make_cca_mode_plot;

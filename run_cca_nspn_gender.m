% RUN CCA - NSPN data (gender)

%%% additional matlab toolboxes required (will need to addpath for each of these)
% FSLNets     http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
% PALM        http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
% nearestSPD  http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

%%% input data
vars = csvread('/Users/maria/Documents/NSPN/docs/NSPN_vars_baseline_110816v5.csv');     
vars(:,sum(vars~=999)<269) = [];  % pre-delete vars with LOADS of missing data
vars = vars(:,2:end);
med = repmat(median(vars),size(vars,1),1); % missing data imputation using the median
vars(vars==999) = med(vars==999);
N = size(vars,1);

%%% prepare permutation scheme using PALM - for more details see:
Nperm = 10000;                                                                       % in the paper we used 100000 but 10000 should be enough
EB = ones(N,1);
PAPset = palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

%%% load netmats
load('/Users/maria/Documents/NSPN/analysis/cca_analysis/netmats_NET_corr.mat');

% load controls
ids = csvread('/Users/maria/Documents/NSPN/docs/NSPN_MRIids_baseline_new.csv'); 
ids = ids(:,1);
NET = NET(ids,:);

%%% Load confounds
varsQconf = csvread('/Users/maria/Documents/NSPN/docs/NSPN_weight_height_baseline.csv');     % weight and height
varsQconf = varsQconf(:,2:end);
age = csvread('/Users/maria/Documents/NSPN/docs/NSPN_age_baseline_cambridge.csv');     % age
age = age(:,2:end);
load('/Users/maria/Documents/NSPN/code/mean_frame_displacement.mat');      % Frame displacement
fd = mean(frame_disp_sub(ids,:),2); 
% varsQconf = [varsQconf, age, fd];
varsQconf = [varsQconf, fd]; % remove age
mri_centre = csvread('/Users/maria/Documents/NSPN/docs/NSPN_MRIcentre_baseline.csv');      % MRI centre
mri_centre = mri_centre(:,2:end);

%%% Choose one gender
gender = csvread('/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_baseline.csv');      % gender
gender = gender(:,2:end);

% g = 1; % male = 1 / female = 0
% NET = NET(gender==g,:);
% vars = vars(gender==g,:);
% varsQconf = varsQconf(gender==g,:);
% mri_centre = mri_centre(gender==g,:);

%%% setup confounds matrix
conf = palm_inormal(varsQconf);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf = nets_normalise(conf);  % add on squared terms and renormalise
% conf = [conf mri_centre gender];
conf = [conf mri_centre];

%%% number of components to feed into CCA
Nkeep = 100;

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1 = nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot = NET1;
NETd = nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu1,ss1,vv1]=nets_svds(NETd,Nkeep); % SVD reduction

%%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values
badvars=[];
for i=1:size(vars,2)
  Y=vars(:,i); grotKEEP=~isnan(Y);  
  grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?
  if (sum(grotKEEP)>250) & (std(Y(grotKEEP))>0) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95) & (grot<100)
  else
    [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
    badvars=[badvars i];
  end
end

%%% get list of which SMs to feed into CCA
varskeep=setdiff([1:size(vars,2)],badvars);                                                                

%%% "impute" missing vars data - actually this avoids any imputation
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
[uu,dd]=eigs(varsdCOV2,Nkeep);  % SVD (eigs actually)
uu2=uu-conf*(pinv(conf)*uu);    % deconfound again just to be safe

%%% CCA
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,uu2);

%%% CCA permutation testing
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,uu2(PAPset(:,j),:));
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%%% netmat weights for CCA mode 1
grotAA = corr(grotU(:,1),NET)';
%  % or
% grotAAd = corr(grotU(:,1),NETd(:,1:size(NET,2)))'; % weights after deconfounding

%%% SM weights for CCA mode 1
grotBB = corr(grotV(:,1),palm_inormal(vars),'rows','pairwise');
 % or 
varsgrot=palm_inormal(vars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grotBBd = corr(grotV(:,1),varsgrot,'rows','pairwise')'; % weights after deconfounding

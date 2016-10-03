% Input data - healthy
% -------------------------------------------------------------------------

% Cognitive/Clinical variables - healthy
% vars = csvread('/Users/maria/Documents/NSPN/docs/NSPN_vars_baseline_110816v5.csv');
load('/Users/maria/Documents/NSPN/docs/NSPNvarsbaseline110816v5.mat')
vars = NSPNvarsbaseline110816v5;

% Connectivity (netmats - full correlation, NET)
load('/Users/maria/Documents/NSPN/analysis/cca_analysis/netmats_NET_corr.mat');

% Load only healthy scans
ids = csvread('/Users/maria/Documents/NSPN/docs/NSPN_MRIids_baseline_new.csv'); 
ids = ids(:,1);
NET = NET(ids,:);

% Confounds - healthy
% -------------------------------------------------------------------------

% Weight, height and waist
varsQconf = csvread('/Users/maria/Documents/NSPN/docs/NSPN_weight_height_waist_baseline.csv');     
varsQconf = varsQconf(:,2:end);

% Age
age = csvread('/Users/maria/Documents/NSPN/docs/NSPN_age_baseline_cambridge.csv');    
age = age(:,2:end);

% Mean frame displacement
load('/Users/maria/Documents/NSPN/code/mean_frame_displacement.mat');      
fd = mean(frame_disp_sub(ids,:),2); 

% MRI centre
mri_centre = csvread('/Users/maria/Documents/NSPN/docs/NSPN_MRIcentre_baseline.csv');      
mri_centre = mri_centre(:,2:end);

% Gender
gender = csvread('/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_baseline.csv');      
gender = gender(:,2:end);

% Psychiatric scores
psy_scores = csvread('/Users/maria/Documents/NSPN/docs/NSPN_psy_baseline.csv');   
% psy_scores(250,:) = NaN;
med = repmat(median(psy_scores),size(psy_scores,1),1); % missing data imputation using the median
psy_scores(isnan(psy_scores)) = med(isnan(psy_scores));

% Load depressed 
% -------------------------------------------------------------------------
% varsdep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_vars_depressed.csv');
load('/Users/maria/Documents/NSPN/docs/NSPNvarsdepressed.mat')
varsdep = NSPNvarsdepressed;

% load depressed Ids
ids_d = csvread('/Users/maria/Documents/NSPN/docs/NSPN_IDs_depressed.csv'); 
ids_d = ids_d(:,2);
NETdep = NET(ids_d,:);

% Load confounds - depressed
% -------------------------------------------------------------------------

% Weight, height and waist
varsQconfdep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_weight_height_waist_depressed.csv');     % weight and height
varsQconfdep = varsQconfdep(:,2:end);

% Age
agedep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_age_depressed.csv');     % age
agedep = agedep(:,2:end);

% Mean frame displacement
fddep = mean(frame_disp_sub(ids_d,:),2); 

% MRI centre
mri_centred_dep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_MRIcentre_depressed.csv');      % MRI centre
mri_centred_dep = mri_centred_dep(:,2:end);

% Gender
genderdep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_depressed.csv');      % gender
genderdep = genderdep(:,2:end);

% Psychiatric scores
psy_scores_dep = csvread('/Users/maria/Documents/NSPN/docs/NSPN_psy_depressed.csv');   

% Concatenate data (comment for healthy only)
% -------------------------------------------------------------------------

switch cohort
    case 'healthy'
        disp('Healthy cohort!');
    case 'all'
        disp('Healthy + Depressed cohort!');
        % Concatenate nets
        NET = [NET; NETdep];
        
        % Concatenate variables
        vars = [vars; varsdep];
        
        % Concatenate confounds
        varsQconf = [varsQconf; varsQconfdep];
        mri_centre = [mri_centre; mri_centred_dep];
        fd = [fd; fddep];
        age = [age; agedep];
        gender = [gender; genderdep];
        
        psy_scores = [psy_scores; psy_scores_dep];
end

confraw = [varsQconf, fd, mri_centre, age, gender];
figure;subplot(1,3,1);imagesc(NET);subplot(1,3,2);imagesc(vars);subplot(1,3,3);imagesc(confraw);


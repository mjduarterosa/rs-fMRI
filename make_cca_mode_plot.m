load '/Users/maria/Documents/NSPN/docs/vars_label.mat'
load '/Users/maria/Documents/NSPN/analysis/cca_analysis/results_100pca_allconfs_wage_depressed.mat'

% Plot most important variables
ncomp = 30;
[s,si] = sort(grotBB);
figure;plot(ones(size(s(end-ncomp-1:end),1),1),s(end-ncomp-1:end),'r+');
hold on;
plot(ones(size(s(1:ncomp),1),1),s(1:ncomp),'b+');
for i = 1:ncomp, text(1.1,s(i+end-ncomp),char(varslabelsbaseline{si(i+end-ncomp)}),'FontSize',14); end
for i = 1:ncomp, text(1.1,s(i),char(varslabelsbaseline{si(i)}),'FontSize',14); end

% Plot projections CCA1 with color corresponding to variable
% gender = csvread('/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_baseline.csv');      % gender
% gender = gender(:,2:end);
% age = csvread('/Users/maria/Documents/NSPN/docs/NSPN_age_baseline_cambridge.csv');     % age
% age = age(:,2:end);
% c(:)=vars(:,si(end-1));

% When in the model
% c(:)=[ones(299,1);zeros(33,1)]; % controls and depressed
% figure;scatter(grotU(:,1),grotV(:,1),40,c(:),'filled');

% When projected
c(:)=[ones(299,1);zeros(33,1)]; % controls and depressed
figure;scatter([grotU(:,1);Udep(:,1)],[grotV(:,1);Vdep(:,1)],40,c(:),'filled');
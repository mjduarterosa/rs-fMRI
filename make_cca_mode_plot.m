% -------------------------------------------------------------------------
% Plot CCA results - NSPN data (all data - healthy and depressed)
%
% Based on Smith et al. Nature Neuroscience 2016
% Additional matlab toolboxes required (will need to addpath for each of these)
% FSLNets     http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
% PALM        http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
% nearestSPD  http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
% -------------------------------------------------------------------------

% Load variable labels
% -------------------------------------------------------------------------
load '/Users/maria/Documents/NSPN/docs/vars_label.mat'

% Plot most important variables (n = 20)
% -------------------------------------------------------------------------
ncomp = 20;
[s,si] = sort(grotBB);
[sa,sia] = sort(abs(grotBB));

% Make bar plot
scales = {'apsd','bis','cads','dasi','icu','k10','spq','wemwbs','ypq','ctq','wasi'};
clb = hsv(length(scales));
for i = 1:length(varslabelsbaseline)
    for j =1:length(scales)
        if strfind(varslabelsbaseline{i},scales{j})
            colorlabels(i,:) = j;
        end
    end
end
data = [s(1:ncomp),s(end-ncomp+1:end)];
data_low = [sa(1:2*ncomp)];
datalab_low = [sia(1:2*ncomp)];
datalab = [si(1:ncomp),si(end-ncomp+1:end)];
labelsw = cell(1,length(datalab));
labelsw_low = cell(1,length(datalab_low));
for i = 1:length(datalab)
    labelsw{i} = char(varslabelsbaseline{datalab(i)});
    labelsw_low{i} = char(varslabelsbaseline{datalab_low(i)});
end
grp = colorlabels(datalab);
grp_low = colorlabels(datalab_low);
for i = 1:length(scales), grp_labels{i,1} = scales{i}; grp_labels{i,2} = i; end
plot_w_joao(data, labelsw, grp, grp_labels);
plot_w_joao(data_low, labelsw_low, grp_low, grp_labels);

% Plot projections CCA1 with color corresponding to variable
% -------------------------------------------------------------------------
c(:) = age;
c(c<16)=1;
c((c>=16 & c<18))=2;
c((c>=18 & c<20))=3;
c((c>=20 & c<22))=4;
c((c>=22 & c<26))=5;
figure;scatter(grotU(:,1),grotV(:,1),100,c(:),'filled','o');
colorbar;
c(:)=gender;
figure;scatter(grotU(:,1),grotV(:,1),100,c(:),'filled','o');
colorbar
if length(c)>300,
    c(1:299)=1;
    c(300:332)=0;
figure;scatter(grotU(:,1),grotV(:,1),100,c(:),'filled','o');
end
colorbar
    
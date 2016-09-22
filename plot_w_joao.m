function plot_w_joao(w, labels, grp, grp_labels)
%
%   Function to plot weight vectors
%
%   Input: w - weight vector
%          labels - labels for the variables of the weight vector
%          grp - group that each varible belongs to (interger vector)
%          grp_labels - labels for the groups (cell with interger and
%          string) grp_labels(:,1) = str; grp_labels(:,2) = int
%
%
%   Version: 2016-06-01
%__________________________________________________________________________

% Written by Joao Matos Monteiro
% Email: joao.monteiro@ucl.ac.uk


%--- Set colors of the groups
colormap('jet');
label_colors = colormap;
if size(label_colors, 1) < length(unique(grp))
    error('Too many groups, not enough colors to plot them.')
end
% Get a subset of colors from the colorscale
label_colors = label_colors(round(linspace(1,size(label_colors,1),size(grp_labels,1))),:);

%--- Get colors and labels all in one group
for i = 1:size(grp_labels,1)
    grp_labels{i,3} = label_colors(i,:);
end

%--- Plot
hold on
for i = 1:size(grp_labels,1)
    grp_ind = grp_labels{i,2} == grp;
    
    dummy = w;
    dummy(~grp_ind) = 0;
    
    bar(dummy, 'facecolor', grp_labels{i,3});
end
hold off

%--- Cosmetic stuff
hl = legend(grp_labels(:,1), 'Location','NW');
set(hl, 'FontSize', 18);

% xlabel('Variable')
% ylabel('Correlation with CCA projection: variables')

if ~isempty(labels)
    ax = gca;
    tickLocations = 1:length(labels); % change to whatever you want
    set(ax,'xTick',tickLocations,'xTickLabel',labels,'FontSize',18)
    ax.XTickLabelRotation=45;
end

xlim([0, length(w)+1])
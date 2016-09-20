% -------------------------------------------------------------------------
% Main script for cleaning ROI signals from all sites
% -------------------------------------------------------------------------
clear
close all
clc

% Output file
f_dir   = '/Users/maria/Documents/NSPN/analysis/';

% Choose TE time-series
for ts_TE = 1:3,
    
    % ROI files
    files{1} = 'ROIs_conf_nspn_cbu.mat';
    files{2} = 'ROIs_conf_nspn_wbic.mat';
    files{3} = 'ROIs_conf_nspn_ucl.mat';
    files{4} = 'ROIs_sulci_nspn_cbu.mat';
    files{5} = 'ROIs_sulci_nspn_wbic.mat';
    files{6} = 'ROIs_sulci_nspn_ucl.mat';
    files{7} = sprintf('ROIs_move_nspn_cbu_%s.mat',num2str(ts_TE));
    files{8} = sprintf('ROIs_move_nspn_wbic_%s.mat',num2str(ts_TE));
    files{9} = sprintf('ROIs_move_nspn_ucl_%s.mat',num2str(ts_TE));
    
    % Get 2nd TE time-series and confounds
    fprintf('Obtaining data >>>>\n');
    kc = 1;
    kr = 1;
    km = 1;
    ks = 0;
    nsites = 3;
    tsub = 0;
    time_series_conf = [];
    time_series_rois = [];
    time_series_move = [];
    for i = 1:length(files)
        if i <=3
            load(fullfile(f_dir,files{i}));
            nsub = size(ROIc,1);
            nreg = size(ROIc,3);
            tsub = tsub + nsub;
            for j = 1:nsub,
                if j~=201 % Accomodate for bad subject
                    for r = 1:nreg,
                        time_series_conf(:,kc) = ROIc{j,ts_TE,r};
                        kc = kc+1;
                    end
                end
            end
        else
            if i <= 6
                load(fullfile(f_dir,files{i}));
                nsub = size(ROI,1);
                nreg = size(ROI,3);
                for j = 1:nsub,
                    if j~=201
                        for r = 1:nreg,
                            time_series_rois(:,kr) = ROI{j,ts_TE,r};
                            kr = kr+1;
                        end
                    end
                end
            else
                load(fullfile(f_dir,files{i}));
                nsub = size(mov_param,3);
                nreg = size(mov_param,2);
                for j = 1:nsub,
                    
                    frame_disp = zeros(size(mov_param,1),1);
                    if j~=201
                        ks = ks + 1;
                        for r = 1:nreg,
                            time_series_move(:,km) = mov_param(:,r,j);
                            frame_disp = frame_disp + [time_series_move(2:end,km); time_series_move(end,km);] - time_series_move(:,km);
                            km = km+1;
                        end
                    end
                    frame_disp_sub(ks,ts_TE) = mean(frame_disp);

                end
                
            end
        end
    end
    fprintf('Done.');
    
    
    
%     % Filter data
%     fprintf('\nFiltering data >>>>\n');
%     time_series_conf_filt = y_IdealFilter(time_series_conf, 2.42, [0.008 0.1]);
%     time_series_rois_filt = y_IdealFilter(time_series_rois, 2.42, [0.008 0.1]);
%     time_series_move_filt = y_IdealFilter(time_series_move, 2.42, [0.008 0.1]);
%     fprintf('Done.\n');
%     
%     % Remove confounds
%     nreg = 137;
%     k = 1;
%     time_series_rois_filt_corr = [];
%     Rcorr = [];
%     for i = 1:tsub-1
%         fprintf('\nSubject %d out of %d >>>>\n',i,tsub-1);
%         tmp_cov = [];
%         for r = 1:nreg
%             R = time_series_rois_filt(:,k);
%             X0 = [time_series_conf_filt(:,(i-1)*2+1:i*2) time_series_move_filt(:,(i-1)*6+1:i*6) ones(length(R),1)];
%             R0 = eye(length(R))-X0*pinv(X0);
%             time_series_rois_filt_corr(:,k) = R0*R;
%             if i <= 53,
%                 ROIcbu{i,ts_TE,r} = time_series_rois_filt_corr(:,k);
%             else
%                 if i <= (53 + 247)
%                     ROIwbic{i-53,ts_TE,r} = time_series_rois_filt_corr(:,k);
%                 else
%                     ROIucl{i-(53 + 247),ts_TE,r} = time_series_rois_filt_corr(:,k);
%                 end
%             end
%             k = k + 1;
%         end
%     end
    
end

% Save clean filesROI
% save(fullfile(f_dir,'ROIs_clean_nspn_cbu.mat'),'ROIcbu','-v7.3');
% save(fullfile(f_dir,'ROIs_clean_nspn_wbic.mat'),'ROIwbic','-v7.3');
% save(fullfile(f_dir,'ROIs_clean_nspn_ucl.mat'),'ROIucl','-v7.3');

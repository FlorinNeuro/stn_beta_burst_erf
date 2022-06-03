%% This script is for a t test against baseline on the cortical activity
% For this script Brainstorm needs to be started
% Matthias Sure
clear
path_brainstorm_data = 'E:\brainstorm_db\Cortical_activity_bipolar_first_z_75_200ms\data\';
%define the time points of interest
PoI = {'start','peak','random'};
%which averages are of interest
bursts_duration = {};
bursts_duration_name = {};
%Hemisphere
Medication = {'OFF','ON'};
%Was the data temporly smoothed?
smoothing = 5;%ms
for iSide = 1 : size(Side,2)
    for iMed = 1 : size(Medication,2)
        for iPoI = 1 : 2% size(PoI,2)
            for iBurst = 1 : size(bursts_duration,2)
                if (iPoI == 3 && iSide > 2) || (iPoI == 3 && iBurst > 2)
                    continue
                end                
                if iPoI == 3
                    target_folder = [PoI{iPoI} '_burts_' Medication{iMed} '_tsmoothed'];
                else
                    target_folder = [PoI{iPoI} '_' bursts_duration_name{iBurst} '_burts_' Side{iSide} '_STN_' Medication{iMed} '_tsmoothed'];
                end
                cd([path_brainstorm_data 'Group_analysis\' target_folder])
                Files = dir;
                sFiles = {};
                n_good_Files = 0;
                folder = ['Group_analysis\' target_folder '\'];
                for iFile = 1 : size(Files,1)
                    if contains(Files(iFile).name,'results') && ~contains(Files(iFile).name,'brainstorm') && ~contains(Files(iFile).name,'presults') && ~contains(Files(iFile).name,'pthresh')&& ~contains(Files(iFile).name,'base')
                        n_good_Files = n_good_Files + 1;
                        sFiles{1,n_good_Files} = [folder  Files(iFile).name];                              
                    end
                end
                if isfile(['presults_t_test_baseline_' target_folder '_' mat2str(n_good_Files) '_files.mat'])
                    continue
                end                
                % Process: t-test baseline [-200ms,100ms]          H0:(X=Baseline), H1:(X<>Baseline)
                sFiles = bst_process('CallProcess', 'process_test_baseline', sFiles, [], ...
                    'baseline',      [-0.2, -0.1], ...
                    'timewindow',    [-0.2, 0.1], ...
                    'scoutsel',      {}, ...
                    'scoutfunc',     3, ...  % PCA
                    'isnorm',        0, ...
                    'avgtime',       0, ...
                    'test_type',     'ttest_baseline', ...  % Student's t-test vs baseline        X~N(m,v)Y = mean_trials(X)        Y~N(m,v)t = (Y - mean_time(Y(baseline)) / std_time(Y(baseline)))df = Nbaseline - 1 = length(baseline) - 1
                    'tail',          'two');  % Two-tailed)
                t = load([path_brainstorm_data  sFiles.FileName]);
                t.Comment = ['t_test_baseline_' target_folder '_' mat2str(n_good_Files) '_files'];
                save([path_brainstorm_data  folder 'presults_' t.Comment '.mat'],'-struct','t')
                delete([path_brainstorm_data  sFiles.FileName]);                
            end
        end
    end
end
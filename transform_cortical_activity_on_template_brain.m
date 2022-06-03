%% This script is for transforming subject level data on a template brain
% For this script Brainstorm needs to be started
% Matthias Sure
clear
%define the subject ID
subject_ID = {'S001','S003','S004', ...};
%define the path for the data
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
                for iSubj = 12% : size(subjects,2)
                    if iPoI == 3
                        target_folder = [PoI{iPoI} '_burts_' Medication{iMed}];
                    else
                        target_folder = [PoI{iPoI} '_' bursts_duration_name{iBurst} '_burts_' Side{iSide} '_STN_' Medication{iMed}];
                    end     
                    if ~exist([path_brainstorm_data subject_ID{iSubj} '_PD_peri\' target_folder], 'dir')
                        continue
                    end
                    cd([path_brainstorm_data subject_ID{iSubj} '_PD_peri\' target_folder])
                    Files = dir;
                    for iFile = 1 : size(Files,1)
                        if contains(Files(iFile).name,'results') && ~contains(Files(iFile).name,'brainstorm') && ~contains(Files(iFile).name,'tsmooth') %&& contains(Files(iFile).name,'base')
                            result_file = [path_brainstorm_data 'Group_analysis\' target_folder '_tsmoothed\' Files(iFile).name(1:end-4) '_' subject_ID{iSubj} '_PD_peri.mat'];
                            if exist(result_file, 'file') == 2
                                continue
                            elseif iSide == 1 && iSubj == 1
                                continue
                            else
                                sFiles = {[subject_ID{iSubj} '_PD_peri/' target_folder '/' Files(iFile).name]};
                                % Process: Project on default anatomy: surface
                                sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                                    'headmodeltype', 'surface');  % Cortex surface                         
                            end                                
                        end
                    end
                end
            end
        end
    end
end
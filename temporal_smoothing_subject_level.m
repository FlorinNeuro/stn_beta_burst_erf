%% temporal smoothing with a gaussianfilter
%This script is for temporal smoothing of time series with a
%gaussian filter.
%It is based on the data structure of Brainstorm.

%Matthias Sure

clear
%define the subject ID
subject_ID = {'S001','S003','S004', ...};
%define the path for the data
path_brainstorm_data = ['.../brainstorm_db/Cortical_activity_bipolar/data/'];
%define the time points of interest
PoI = {'start','peak','random'};
%which averages are of interest
bursts_duration = {};
bursts_duration_name = {};
%Hemisphere
Side = {'left','right'};
Medication = {'OFF','ON'};
%What is the window length for the smoothing; Data was sampled with 1000 Hz 
smoothing = 5;%ms
for iSide = 2% : size(Side,2)
    for iMed = 2%1 : size(Medication,2)
        for iPoI = 2%1: 2%1 : size(PoI,2)
            for iBurst = 2%1 : size(bursts_duration,2)
                for iSubj = 12% : size(subjects,2)
                    if iPoI == 3
                        target_folder = [path_brainstorm_data subject_ID{iSubj} '_PD_peri\' PoI{iPoI} '_burts_' Medication{iMed}];
                    else
                        target_folder = [path_brainstorm_data subject_ID{iSubj} '_PD_peri\' PoI{iPoI} '_' bursts_duration_name{iBurst} '_burts_' Side{iSide} '_STN_' Medication{iMed}];
                    end
                    if exist(target_folder, 'dir')
                        cd(target_folder)
                        Files = dir;
                        for iFile = 1 : size(Files,1)
                            if contains(Files(iFile).name,'results') && ~contains(Files(iFile).name,'brainstorm') && ~contains(Files(iFile).name,'tsmooth') && ~contains(Files(iFile).name,'base')
                                data = load(Files(iFile).name);
                                % Process: Spatial smoothing (1.00)
                                data.ImageGridAmp = smoothdata(data.ImageGridAmp,2,'gaussian',smoothing);
                                data.Comment = [data.Comment '_tsmooth' mat2str(smoothing)];
                                save(['results_' data.Comment '.mat'],'-struct','data')                   
                            end            
                        end                        
                    end
                end
            end
        end
    end
end
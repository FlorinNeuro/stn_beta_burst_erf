function Average_cortical_source_activity_bipolar_fun(iSubject)
%% This Script should average the cortical activity time-locked to a event
%This script is used to average cortical activity around specific events. 
%It is assumed that the data is stored in a Brainstorm structure. 
%Brainstorm functions are used. The events must have been determined 
%beforehand. For each event, the cortical activity is averaged for a 
%certain period of time by determining the inverse solution for the MEG 
%sensors. This activity is corrected against a baseline. The activity at 
%the time of each event is averaged at the end.

%clear
warning('off','all')
Med_state = 'OFF';
%% Define the important paths
path_brainstorm_project = ['/gpfs/project/sumat100/brainstorm_db/STN_COR_' Med_state '/'];
path_brainsotrm_data=['/gpfs/project/sumat100/brainstorm_db/STN_COR_' Med_state '/data/'];
path_brainstorm_code = '/gpfs/project/sumat100/brainstorm3/';
addpath(genpath(path_brainstorm_code));

BrainstormDbDir = 'E:/brainstorm_db';

%%how to start brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
ProtocolName = ['STN_COR_' Med_state];
sProtocol.Comment = ProtocolName;
sProtocol.SUBJECTS = [path_brainstorm_project 'anat'];
sProtocol.STUDIES = [path_brainstorm_project 'data'];
db_edit_protocol('load',sProtocol);
% Get the protocol index
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

% %% Define the files you want to investigate
subject_ID = {'S001','S003','S004','S008','S009','S010','S011','S012','S013','S014','S016','S020','S021','S022','S023','S025','S027','S028','S029','S032','S033','S034','S036','S037','S038','S041','S043'};
switch Med_state
    case 'ON'
        % ON
        subject = {...
            'S001_PD_peri' {'S001PD_rest_peri_ON_001_notch_high', ...
            'S001PD_rest_peri_ON_002_notch_high', ...
            'S001PD_rest_peri_ON_003_notch_high'}; ...
            'S003_PD_peri' {'S003_PD_ON_peri_04052017_rest_001_notch_high_notch', ...
            'S003_PD_ON_peri_04052017_rest_002_notch_high_notch', ...
            'S003_PD_ON_peri_04052017_rest_003_notch_high_notch'}; ...
            'S004_PD_peri' {'S004_PD_ON_peri_27042017_rest_001_notch_high_notch', ...
            'S004_PD_ON_peri_27042017_rest_002_notch_high_notch', ...
            'S004_PD_ON_peri_27042017_rest_003_notch_high_notch'}; ...
            'S008_PD_peri' {'bk_PD_ON_peri_04092017_rest_001_notch_high_notch', ...
            'bk_PD_ON_peri_04092017_rest_002_notch_high_notch'}; ...
            'S009_PD_peri' {'S009_PD_ON_peri_16102017_rest_001_notch_high_notch', ...
            'S009_PD_ON_peri_16102017_rest_002_notch_high_notch', ...
            'S009_PD_ON_peri_16102017_rest_003_notch_high_notch'}; ...
            'S010_PD_peri' {'S010_PD_ON_peri_27102017_rest_001_notch_high', ...
            'S010_PD_ON_peri_27102017_rest_002_notch_high', ...
            'S010_PD_ON_peri_27102017_rest_003_notch_high'}; ...
            'S011_PD_peri' {'S011_PD_ON_peri_30112017_rest_001_notch_02_high', ...
            'S011_PD_ON_peri_30112017_rest_002_notch_02_high', ...
            'S011_PD_ON_peri_30112017_rest_003_notch_02_high'}; ...
            'S013_PD_peri' {}; ... 
            'S013_PD_peri' {'S013_PD_ON_peri_22032018_rest_001_notch_high_notch', ...
            'S013_PD_ON_peri_2203208_rest_002_notch_high_notch', ...
            'S013_PD_ON_peri_22032018_rest_003_notch_high_notch'}; ...    
            'S014_PD_peri' {'S014_PD_ON_peri_2303018_rest_001_notch_high_notch', ...
            'S014_PD_ON_peri_23032018_rest_002_notch_high_notch', ...
            'S014_PD_ON_peri_23032018_rest_003_notch_high_notch'}; ...
            'S016_PD_peri' {'S016_PD_ON_peri_22112018_rest_001_notch_high', ...
            'S016_PD_ON_peri_22112018_rest_002_notch_high', ...
            'S016_PD_ON_peri_22112018_rest_003_notch_high'};
            'S020_PD_peri' {'S020_PD_ON_peri_18102018_rest_001_notch_high', ...
            'S020_PD_ON_peri_18102018_rest_002_notch_high', ...
            'S020_PD_ON_peri_18102018_rest_003_notch_high'};
            'S021_PD_peri' {'S021_PD_ON_peri_01112018_rest_001_notch_high', ...
            'S021_PD_ON_peri_01112018_rest_002_notch_high', ...
            'S021_PD_ON_peri_01112018_rest_003_notch_high'};
            'S022_PD_peri' {'S022_PD_ON_peri_15112018_rest_001_notch_high', ...
            'S022_PD_ON_peri_15112018_rest_002_notch_high', ...
            'S022_PD_ON_peri_15112018_rest_003_notch_high'};
            'S023_PD_peri' {'S023_PD_ON_peri_29112018_rest_001_notch_high', ...
            'S023_PD_ON_peri_29112018_rest_002_notch_high', ...
            'S023_PD_ON_peri_29112018_rest_003_notch_high'};    
            'S025_PD_peri' {'S025_PD_ON_peri_06122018_rest_001_notch_high', ...
            'S025_PD_ON_peri_06122018_rest_002_notch_high', ...
            'S025_PD_ON_peri_06122018_rest_003_notch_high'}
            'S027_PD_peri' {'S027_PD_ON_peri_18012019_rest_001_notch_high', ...
            'S027_PD_ON_peri_18012019_rest_002_notch_high', ...
            'S027_PD_ON_peri_18012019_rest_003_notch_high'}
            'S028_PD_peri' {'S028_PD_ON_peri_31012019_rest_001_notch_high', ...
            'S028_PD_ON_peri_31012019_rest_002_notch_high', ...
            'S028_PD_ON_peri_31012019_rest_003_notch_high'}    
            'S029_PD_peri' {'S029_PD_ON_peri_14022019_rest_01_notch_high', ...
            'S029_PD_ON_peri_14022019_rest_002_notch_high', ...
            'S029_PD_ON_peri_14022019_rest_003_notch_high'}
            'S032_PD_peri' {'S032_PD_ON_peri_07042019_rest_001_notch_high', ...
            'S032_PD_ON_peri_07042019_rest_002_notch_high', ...
            'S032_PD_ON_peri_07042019_rest_003_notch_high'}
            'S033_PD_peri' {'S033_PD_ON_peri_28042019_rest_001_notch_high', ...
            'S033_PD_ON_peri_28042019_rest_002_notch_high', ...
            'S033_PD_ON_peri_28042019_rest_003_notch_high'} 
            'S034_PD_peri' {'S034_PD_ON_peri_03052019_rest_001_notch_high', ...
            'S034_PD_ON_peri_03052019_rest_002_notch_high', ...
            'S034_PD_ON_peri_03052019_rest_003_notch_high'} 
            'S036_PD_peri' {'S036_PD_ON_peri_04072019_rest_001_notch_high', ...
            'S036_PD_ON_peri_04072019_rest_002_notch_high'}
            'S037_PD_peri' {'S037_PD_ON_peri_09082019_rest_001_notch_high', ...
            'S037_PD_ON_peri_09082019_rest_002_notch_high'}     
            'S038_PD_peri' {'S038_PD_ON_peri_08082019_rest_001_notch_high', ...
            'S038_PD_ON_peri_08082019_rest_002_notch_high', ...
            'S038_PD_ON_peri_08082019_rest_003_notch_high'}
            'S041_PD_peri' {'S041_PD_ON_peri_05092019_rest_001_notch_high', ...
            'S041_PD_ON_peri_05092019_rest_002_notch_high', ...
            'S041_PD_ON_peri_05092019_rest_003_notch_high'}
            'S043_PD_peri' {'S043_PD_ON_peri_18102019_rest_001_notch_high', ...
            'S043_PD_ON_peri_18102019_rest_002_notch_high', ...
            'S043_PD_ON_peri_18102019_rest_003_notch_high'}};
    case 'OFF'
        %OFF
        subject = {...
            'S001_PD_peri' {'S001PD_rest_peri_OFF_003_notch_notch_high_notch', ...
            'S001PD_rest_peri_OFF_02_notch_notch_high_notch', ...
            'S001PD_rest_peri_OFF_1_notch_notch_high_notch'}; ...
            'S003_PD_peri' {'S003_PD_OFF_peri_04052017_rest_001_notch_notch_high_notch', ...
            'S003_PD_OFF_peri_04052017_rest_002_notch_notch_high_notch', ...
            'S003_PD_OFF_peri_04052017_rest_003_notch_notch_high_notch'}; ...
            'S004_PD_peri' {'S004_PD_OFF_peri_27042017_reset_003_notch_notch_high_notch', ...
            'S004_PD_OFF_peri_27042017_rest_001_notch_notch_high_notch', ...
            'S004_PD_OFF_peri_27042027_rest_002_notch_notch_high_notch'}; ...
            'S008_PD_peri' {'bk_PD_OFF_peri_04092017_rest_001_notch_notch_high_notch', ...
            'bk_PD_OFF_peri_04092017_rest_002_notch_notch_high_notch'}; ...
            'S009_PD_peri' {'S009_PD_OFF_peri_16102017_rest_001_notch_notch_high_notch', ...
            'S009_PD_OFF_peri_16102017_rest_002_notch_notch_high_notch', ...
            'S009_PD_OFF_peri_16102017_rest_003_notch_notch_high_notch'}; ...
            'S010_PD_peri' {'S010_PD_OFF_peri_27102017_rest_001_notch_notch_high_notch', ...
            'S010_PD_OFF_peri_27102017_rest_002_notch_notch_high_notch', ...
            'S010_PD_OFF_peri_27102017_rest_003_notch_notch_high_notch'}; ...
            'S011_PD_peri' {'S011_PD_OFF_peri_30112017_rest_001_notch_notch_high_notch', ...
            'S011_PD_OFF_peri_30112017_rest_002_notch_notch_high_notch', ...
            'S011_PD_OFF_peri_30112017_rest_003_notch_notch_high_notch'}; ...
            'S012_PD_peri' {'S012_PD_OFF_peri_09022018_rest_001_notch_notch_high_notch', ...
            'S012_PD_OFF_peri_09022018_rest_002_notch_notch_high_notch'}; ...
            'S013_PD_peri' {}; ... 
            'S014_PD_peri' {'S014_PD_OFF_peri_23032018_rest_001_notch_high_notch', ...
            'S014_PD_OFF_peri_23032018_rest_002_notch_high_notch', ...
            'S014_PD_OFF_peri_23032018_rest_003_notch_high_notch'}; ...
            'S016_PD_peri' {'S016_PD_OFF_peri_22112018_rest_001_notch_high', ...
            'S016_PD_OFF_peri_22112018_rest_002_notch_high', ...
            'S016_PD_OFF_peri_22112018_rest_003_notch_high'};
            'S020_PD_peri' {'S020_PD_OFF_peri_18102018_rest_001_notch_high', ...
            'S020_PD_OFF_peri_18102018_rest_002_notch_high', ...
            'S020_PD_OFF_peri_18102018_rest_003_notch_high'};
            'S021_PD_peri' {'S021_PD_OFF_peri_01112018_rest_001_notch_high', ...
            'S021_PD_OFF_peri_01112018_rest_002_notch_high', ...
            'S021_PD_OFF_peri_01112018_rest_003_notch_high'};
            'S022_PD_peri' {'S022_PD_OFF_peri_15112018_rest_001_notch_high', ...
            'S022_PD_OFF_peri_15112018_rest_002_notch_high', ...
            'S022_PD_OFF_peri_15112018_rest_003_notch_high'};
            'S023_PD_peri' {'S023_PD_OFF_peri_29112018_rest_001_notch_high', ...
            'S023_PD_OFF_peri_29112018_rest_002_notch_high', ...
            'S023_PD_OFF_peri_29112018_rest_003_notch_high'};
            'S025_PD_peri' {'S025_PD_OFF_peri_06122018_rest_001_notch_high', ...
            'S025_PD_OFF_peri_06122018_rest_002_notch_high', ...
            'S025_PD_OFF_peri_06122018_rest_003_notch_high'}
            'S027_PD_peri' {'S027_PD_OFF_peri_18012019_rest_001_notch_high', ...
            'S027_PD_OFF_peri_18012019_rest_002_notch_high', ...
            'S027_PD_OFF_peri_18012019_rest_003_notch_high', ...
            'S027_PD_OFF_peri_18012019_rest_004_notch_high'}
            'S028_PD_peri' {'S028_PD_OFF_peri_31012019_rest_001_notch_high', ...
            'S028_PD_OFF_peri_31012019_rest_002_notch_high', ...
            'S028_PD_OFF_peri_31012019_rest_003_notch_high'}    
            'S029_PD_peri' {'S029_PD_OFF_peri_14022019_rest_001_notch_high', ...
            'S029_PD_OFF_peri_14022019_rest_002_notch_high', ...
            'S029_PD_OFF_peri_14022019_rest_003_notch_high'}
            'S032_PD_peri' {'S032_PD_OFF_peri_07042019_rest_001_notch_high', ...
            'S032_PD_OFF_peri_07042019_rest_002_notch_high', ...
            'S032_PD_OFF_peri_07042019_rest_003_notch_high'}
            'S033_PD_peri' {'S033_PD_OFF_peri_28042019_rest_001_notch_high', ...
            'S033_PD_OFF_peri_28042019_rest_002_notch_high', ...
            'S033_PD_OFF_peri_28042019_rest_003_notch_high'} 
            'S034_PD_peri' {'S034_PD_OFF_peri_03052019_rest_001_notch_high', ...
            'S034_PD_OFF_peri_03052019_rest_002_notch_high', ...
            'S034_PD_OFF_peri_03052019_rest_003_notch_high'} 
            'S036_PD_peri' {'S036_PD_OFF_peri_04072019_rest_001_notch_high', ...
            'S036_PD_OFF_peri_04072019_rest_002_notch_high', ...
            'S036_PD_OFF_peri_04072019_rest_003_notch_high', ...
            'S036_PD_OFF_peri_04072019_rest_004_notch_high'}
            'S037_PD_peri' {'S037_PD_OFF_peri_09082019_rest_001_notch_high', ...
            'S037_PD_OFF_peri_09082019_rest_002_notch_high'}     
            'S038_PD_peri' {'S038_PD_OFF_peri_08082019_rest_001_notch_high', ...
            'S038_PD_OFF_peri_08082019_rest_002_notch_high', ...
            'S038_PD_OFF_peri_08082019_rest_003_notch_high'}
            'S041_PD_peri' {'S041_PD_OFF_peri_05092019_rest_001_notch_high', ...
            'S041_PD_OFF_peri_05092019_rest_002_notch_high', ...
            'S041_PD_OFF_peri_05092019_rest_003_notch_high'}
            'S043_PD_peri' {'S043_PD_OFF_peri_18102019_rest_001_notch_notch_high', ...
            'S043_PD_OFF_peri_18102019_rest_002_notch_high', ...
            'S043_PD_OFF_peri_18102019_rest_003_notch_high'}};
end


% Name the frequency bands for burst detection
freqband_name={'Morlet_peak';'Morlet_lBB'; 'Morlet_hBB'};
% starting, peak and end point of burst
event_pos = {'start','peak','end'};
%around wich time point should be averaged
Event_of_Interest = 2; % 2 = peak
% how much time should be between bursts and bad segments
time_intervall = 0.8;%800ms
% which time period before and after the burst event should considered
time_before_event = 0.100;% 100 ms
time_after_event = 0.100;% 100 ms
min_duration = 0.200;%minimum duration of the burst selected
max_duration = 100.0;%maximum duration of the burst selected
% length of baseline
baseline = 0.100;% 100 ms
% respective to which burst event the baseline should be determined
Baseline_of_Interest = 1;
% how much time before this event the baseline should begin
Basline_start = 0.200;% 200 ms
% create a string to get information about the averaged period
average_period = [mat2str(time_before_event*1000) '_' event_pos{Event_of_Interest} '_' mat2str(time_after_event*1000) 'ms_minDur_' mat2str(min_duration*1000) '_maxDur_' mat2str(max_duration*1000) 'ms'];

% determine the detection setting for the detected bursts
thresh_prc = 75;
detection_identifier = 'first_z_bipolar';
Detection_method = [detection_identifier '_' mat2str(thresh_prc)];
peak_freq = [21,15,16,17,24,18,24,20,0,0,0,17,23,27,0,0,0,0,0,0,15,16,15,20,14,15,30; ...
             21,23,16,15,24,21,24,20,23,27,23,17,0,27,22,20,0,27,0,0,15,0,0,18,0,0,15];
%get the index for the burst events
burst_start = 1;
burst_peak = 2;
burst_end = 3;
%Hemispheres of the bilateral LFP
Hemispheres = {'right_STN','left_STN'};



for iFreq = 1 %: length(freqband_name)
    for iEEG = 1 : 2%length(EEG)
        for iSubj= iSubject%loop over subjects
            %% if you want to average accross beta peak bursts and their is now peak than skip
            if peak_freq(iEEG,iSubj) == 0 && iFreq == 1
                continue
            end
            % reload subject folder
            [~, BST_Subject] = bst_get('Subject', subject{iSubj,1});
            db_reload_conditions(BST_Subject);
            %define frequency band
            freqband = [peak_freq(iEEG,iSubj)-3 peak_freq(iEEG,iSubj)+3];
            % check if the burst average already exist
            if exist([path_brainsotrm_data subject{iSubj,1} '/@intra/results_PNAI_avg_'  Detection_method '_' average_period '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{iSubj} '_Bipolar_' Hemispheres{iEEG} '_' Med_state '.mat'], 'file') == 2
                 continue
            end
            % get the path to the beta bursts
            Events_path = ['/gpfs/project/sumat100/Output/Event_detection/first_z_bipolar_75/' subject_ID{iSubj} '/'];
            %create cell arays to colect the names of created files which
            %can be deletet at the end
            sSource = [];
            s_etime = [];
            s_source = [];
            sFiles = [];
            %count how many files are created
            n = 0;
            for irun=1:length(subject{iSubj,2}) % loop over runs 
                iavg_n = 0; %count how many bursts are averaged
                %% Set up the arrays
                clear Events Time
                % get the time points of bad segments
                Bad_tmp = [];
                Bad_LFP_tmp = [];
                % get the burst time poins
                peak_Burst = [];
                start_Burst = [];
                end_Burst = [];
                Burst_tmp = [];                
                index_data = [];
                %matrix for the averaged data
                av_source = [];
                bs_files=dir(fullfile([path_brainsotrm_data subject{iSubj,1} '/'  subject{iSubj,2}{irun}  '/']));
                %% load the events
                % get the filename of the Brainstorm file with the
                % Information of the Bad segments
                for i=1:length(bs_files)
                    if ~isempty(strfind(bs_files(i).name, 'block'))
                        index_data(end+1) = i;
                    end
                end
                if size(index_data)>1
                    disp('There were more than one datafile')
                    quit
                end
                % load the Bad segments
                load([path_brainsotrm_data subject{iSubj,1} '/' subject{iSubj,2}{irun} '/' bs_files(index_data(1)).name],'Events','Time');
                % load the beta bursts
                try
                    load([Events_path 'Morlet_Bipolar_' Hemispheres{iEEG} '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{iSubj} '_run_' num2str(irun) '_' Detection_method '_' Med_state '.mat'],'DataMat');
                catch
                    continue
                end
                Bursts = DataMat.F.events(2:end);
                clear DataMat
                %% define the timepoints of the BAD segments and events
                for iEvent=1:length(Events)
                    if strcmp(Events(1,iEvent).label,'BAD')
                        %extract the epochs of the BAD-segments
                        Bad_tmp=(Events(1,iEvent).times);
                    end
                    if strcmp(Events(1,iEvent).label,'BAD_LFP')
                        %extract the epochs of the BAD-LFP-segments
                        Bad_LFP_tmp=(Events(1,iEvent).times);
                    end
                end
                clear iEvent
                for iBurst = 1:size(Bursts,2)
                    %find the correct Burst set
                    if strcmp(Bursts(1,iBurst).label,['Morlet_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' Hemispheres{iEEG} '_' event_pos{burst_peak}])
                        %get the timepoints
                        peak_Burst=(Bursts(1,iBurst).times);
                        Burst_tmp = ones(size(peak_Burst));
                    end
                    if strcmp(Bursts(1,iBurst).label,['Morlet_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' Hemispheres{iEEG} '_' event_pos{burst_start}])
                        %get the timepoints
                        start_Burst=(Bursts(1,iBurst).times);
                    end
                    if strcmp(Bursts(1,iBurst).label,['Morlet_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' Hemispheres{iEEG} '_' event_pos{burst_end}])
                        %get the timepoints
                        end_Burst=(Bursts(1,iBurst).times);
                    end                        
                end
                clear iBurst
                % if no BUrst are present than skip here
                if isempty(Burst_tmp)
                    continue
                end
                %% check if the Bursts are part of BAD epochs
                for iEvent = 1 : size(Burst_tmp,2)
                    %check if the event +- intervall is complete in the
                    %Timespan
                    if start_Burst(iEvent) - time_intervall <= Time(1) || end_Burst(iEvent) + time_intervall >= Time(end)
                        Burst_tmp(iEvent) = 0;
                    end
                    %check if the event +- intervall is colliding with
                    %the BAD epochs
                    %start with BAD later on BAD_LFP
                    if ~isempty(Bad_tmp)
                        for iBAD = 1 : size(Bad_tmp,2)
                            if Bad_tmp(1,iBAD) >= start_Burst(iEvent) - time_intervall && Bad_tmp(1,iBAD) <= end_Burst(iEvent) + time_intervall
                                %begin of BAD is inside Event +- intervall
                                Burst_tmp(iEvent) = 0;
                            elseif Bad_tmp(2,iBAD) >= start_Burst(iEvent) - time_intervall && Bad_tmp(2,iBAD) <= end_Burst(iEvent) + time_intervall
                                %end of BAD is inside Event +- inetvall
                                Burst_tmp(iEvent) = 0;
                            elseif Bad_tmp(1,iBAD) <= start_Burst(iEvent) - time_intervall && Bad_tmp(2,iBAD) >= end_Burst(iEvent) + time_intervall
                                %Event +- intervall is inside BAD
                                Burst_tmp(iEvent) = 0;
                            end                                
                        end
                    end
                    if ~isempty(Bad_LFP_tmp)
                        for iBAD_LFP = 1 : size(Bad_LFP_tmp,2)
                            if Bad_LFP_tmp(1,iBAD_LFP) >= start_Burst(iEvent) - time_intervall && Bad_LFP_tmp(1,iBAD_LFP) <= end_Burst(iEvent) + time_intervall
                                %begin of BAD is inside Event +- intervall
                                Burst_tmp(iEvent) = 0;
                            elseif Bad_LFP_tmp(2,iBAD_LFP) >= start_Burst(iEvent) - time_intervall && Bad_LFP_tmp(2,iBAD_LFP) <= end_Burst(iEvent) + time_intervall
                                %end of BAD is inside Event +- inetvall
                                Burst_tmp(iEvent) = 0;
                            elseif Bad_LFP_tmp(1,iBAD_LFP) <= start_Burst(iEvent) - time_intervall && Bad_LFP_tmp(2,iBAD_LFP) >= end_Burst(iEvent) + time_intervall
                                %Event +- intervall is inside BAD
                                Burst_tmp(iEvent) = 0;
                            end                                
                        end
                    end
                end
                %take only the events which are in good epochs
                peak_Burst = peak_Burst(Burst_tmp == 1);
                start_Burst = start_Burst(Burst_tmp == 1);
                end_Burst = end_Burst(Burst_tmp == 1);
                Burst_tmp = Burst_tmp(Burst_tmp == 1);
                % select the time points over which it should be averaged
                if Event_of_Interest == 1
                    eRipple = start_Burst;
                elseif Event_of_Interest == 2
                    eRipple = peak_Burst;
                elseif Event_of_Interest == 3
                    eRipple = end_Burst;
                end
                eBaseline = [];
                if Baseline_of_Interest == 1
                    eBaseline = start_Burst;
                elseif Baseline_of_Interest == 2
                    eBaseline = peak_Burst;
                elseif Baseline_of_Interest == 3
                    eBaseline = end_Burst;
                end
                  
                %% start with averaging part
                for iEvent = 1 : size(eRipple,2)
                    if Burst_tmp(iEvent) == 1
                        iavg_n = iavg_n + 1;
                        %to reduce computation time you can add a limit of
                        %burst to be averaged
                        if iavg_n >= 5000
                            break
                        end
                        % check if burst duration is colliding with a limit
                        if end_Burst(iEvent) - start_Burst(iEvent) < min_duration || end_Burst(iEvent) - start_Burst(iEvent) > max_duration
                            continue 
                        end
                        sFiles_etime = [];
                        sFiles_source = [];
                        %File which should be used to extract the timespan;
                        sFiles{1} = [subject{iSubj,1} '/' subject{iSubj,2}{irun} '/' bs_files(index_data(1)).name];
                        % Start a new report
                        bst_report('Start', sFiles);
                        %Extract the time window from baseline til end of
                        %burst
                        % Process: Extract time:
                        eRipple(iEvent) = round(eRipple(iEvent),3);%round so the number can be rea
                        sFiles_etime = bst_process('CallProcess', 'process_extract_time', sFiles, [], ...
                            'timewindow', [eBaseline(iEvent)-Basline_start, eRipple(iEvent)+time_after_event], ...
                            'overwrite',  0);
                        %adap the timeline for the extracted data;
                        %event is centerd at t = 0
                        load([ path_brainsotrm_data (sFiles_etime.FileName)],'Time','F')
                        % get the timepoint of the end of the baseline
                        end_episode_base = find(round(Time,4)  == round(eBaseline(iEvent)-Basline_start+baseline,4));
                        % get the timepoint of the start of the burst
                        start_episode_burst = find(round(Time,4)  == round(eRipple(iEvent)-time_before_event,4));
                        % remove time between baseline and burst
                        F = [F(:,1:end_episode_base-1) F(:,start_episode_burst:end)];
                        % save the updadet time line
                        save([ path_brainsotrm_data (sFiles_etime.FileName)],'F','-append')
                        % update the time vector
                        Time_step = round(Time(2)-Time(1),4);
                        Time = round([Time_step:Time_step:Time_step*(size(F,2))]-Time_step*size(F,2)+baseline,4);
                        save([ path_brainsotrm_data (sFiles_etime.FileName)],'Time','-append')
                        clear Time
                        % relaod the bs databse
                        db_reload_conditions(BST_Subject);
                        %mark the files which should be deleted
                        % Process: Add tag: _delnow_
                        sFiles_etime = bst_process('CallProcess', 'process_add_tag', sFiles_etime, [], ...
                            'tag',    ['_delnow_' int2str(iEvent)], ...
                            'output', 2);  % Add to file name
                        %compute the data in the source space
                        % Process: Compute sources [2018]
                        sFiles_inter = bst_process('CallProcess', 'process_inverse_2018', sFiles_etime, [], ...
                            'output',  2, ...  % Kernel only: one per file
                            'inverse', struct(...
                                 'Comment',        'PNAI: MEG ALL', ...
                                 'InverseMethod',  'lcmv', ...
                                 'InverseMeasure', 'nai', ...
                                 'SourceOrient',   {{'fixed'}}, ...
                                 'Loose',          0.2, ...
                                 'UseDepth',       1, ...
                                 'WeightExp',      0.5, ...
                                 'WeightLimit',    10, ...
                                 'NoiseMethod',    'median', ...
                                 'NoiseReg',       0.1, ...
                                 'SnrMethod',      'rms', ...
                                 'SnrRms',         1e-06, ...
                                 'SnrFixed',       3, ...
                                 'ComputeKernel',  1, ...
                                 'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));
                        % Process: DC offset correction: [-200ms,-50ms]
                        sFiles_bl = bst_process('CallProcess', 'process_baseline_norm', sFiles_inter, [], ...
                            'baseline',   [-(time_before_event+baseline), -time_before_event], ...
                            'source_abs', 0, ...
                            'method',     'bl', ...  % DC offset correction:    x_std = x - &mu;
                            'overwrite',  1);
                        %compute the data in the source space
                        % Process: Add tag: _avgnow_
                        sFiles_source = bst_process('CallProcess', 'process_add_tag', sFiles_bl, [], ...
                            'tag',    ['_avgnow_' int2str(iEvent)], ...
                            'output', 2);  % Add to file name
                        % relaod the bs databse
                        db_reload_conditions(BST_Subject);
                        
                        %% start to average
                        if iavg_n == 1
                            av_source = load([path_brainsotrm_data sFiles_source.FileName]);
                            av_source_filename = [path_brainsotrm_data subject{iSubj,1} '/' subject{iSubj,2}{irun} '/results_PNAI_avg_' average_period '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_Bipolar_' Hemispheres{iEEG} '_' subject_ID{iSubj} '_run_' mat2str(irun) '_' Med_state '.mat'];                                                                
                            av_source.Comment = ['PNAI_avg_' average_period '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_Bipolar_' Hemispheres{iEEG} '_' subject_ID{iSubj} '_run_' mat2str(irun) '_' Med_state];
                            save([av_source_filename],'-struct','av_source','-v7.3')
                        elseif iavg_n >= 2
                            av_source_new = load([path_brainsotrm_data sFiles_source.FileName]);
                            % sum the cortical maps and divide them at the end by the total number of bursts 
                            av_source.ImageGridAmp = (av_source.ImageGridAmp + av_source_new.ImageGridAmp);
                            save([av_source_filename],'-struct','av_source','-v7.3')
                        end

                        %colect the file names which should be deleted
                        n = n + 1;
                        s_bl{n} = sFiles_bl.FileName;%sFiles_morlet;
                        s_etime{n} = sFiles_etime;
                        s_source{n}= sFiles_source;

                        %delet the files which are not needed anymore
                        % Process: Delete selected files
                        if ~isempty(sFiles_etime)
                        sFiles_del = bst_process('CallProcess', 'process_delete', sFiles_etime, [], ...
                            'target', 1);  % Delete selected files
                        end
                        if ~isempty(sFiles_source)
                        sFiles_del = bst_process('CallProcess', 'process_delete', sFiles_source, [], ...
                            'target', 1);  % Delete selected files
                        end
                        if ~isempty(sFiles_bl)
                        sFiles_del = bst_process('CallProcess', 'process_delete', sFiles_bl, [], ...
                            'target', 1);  % Delete selected files
                        end
                        if ~isempty(sFiles_inter)
                        sFiles_del = bst_process('CallProcess', 'process_delete', sFiles_inter, [], ...
                            'target', 1);  % Delete selected files
                        end 

                    end
                end
                %keep the number of events included in the average for
                %each run
                Navg_N(irun) = iavg_n;
                if iavg_n ~= 0  
                    av_source.ImageGridAmp = (av_source.ImageGridAmp/iavg_n);
                    save([av_source_filename],'-struct','av_source','-v7.3')
                    sSource{irun} = [av_source_filename];
                elseif iavg_n == 0
                    sSource{irun} = {};
                end

            end
            %Determine the average over the runs for each patient
            if ~isempty(sSource)
                cd([path_brainsotrm_data subject{iSubj,1} '/@intra'])
                for irun=1:length(Navg_N)
                    if isempty(sSource{irun})
                         continue
                    end
                    if irun == 1
                        av_source_fin = load([sSource{irun}]);
                        delete([sSource{irun}])
                        av_source_filename_fin = [path_brainsotrm_data subject{iSubj,1} '/@intra/results_PNAI_avg_' Detection_method '_' average_period '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{iSubj} '_Bipolar_' Hemispheres{iEEG} '_' Med_state '.mat'];                            
                        save([av_source_filename_fin],'-struct','av_source_fin','-v7.3')
                        av_source_fin.Comment = (['Power_avg_' Detection_method '_' average_period '_' freqband_name{iFreq} '_Bipolar_' Hemispheres{iEEG}]);
                    elseif irun >= 2
                        if ~exist('av_source_fin','var')
                            av_source_fin = load([sSource{irun}]);
                            delete([sSource{irun}])
                            av_source_filename_fin = [path_brainsotrm_data subject{iSubj,1} '/@intra/results_PNAI_avg_'  Detection_method '_' average_period '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{iSubj} '_Bipolar_' Hemispheres{iEEG} '_' Med_state '.mat'];                            
                            av_source_fin.Comment = (['Power_avg_'  Detection_method '_' average_period '_' freqband_name{iFreq} '_Bipolar_' Hemispheres{iEEG}]);                                
                            save([av_source_filename_fin],'-struct','av_source_fin','-v7.3')
                        else
                            av_sourc_fin_new = load([sSource{irun}]);
                            delete([sSource{irun}])
                            av_source_fin.ImageGridAmp = (av_source_fin.ImageGridAmp .* (sum(Navg_N(1:irun-1))) + av_sourc_fin_new.ImageGridAmp)./(sum(Navg_N(1:irun)));
                            av_source_fin.Comment = (['Power_avg_'  Detection_method '_' average_period '_' freqband_name{iFreq} '_Bipolar_' Hemispheres{iEEG}]);
                            save([av_source_filename_fin],'-struct','av_source_fin','-v7.3')                                
                        end
                    end                       
                end
                clear av_morlet_fin av_source_fin

                % Save and display report
                ReportFile = bst_report('Save', sFiles);


            end
        clear Navg_N

        end
    end
end

brainstorm stop

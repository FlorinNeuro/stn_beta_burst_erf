%% Reorder burst information
% This script is used to write the detected bursts into a structure. 
% Previously, the start and end times were noted for each run of a patient.
% Here a structure is created which combines the burst duration and the 
% maximum amplitude of all passes of a patient.This script is applicable to
% detection with bipolar reref.
% Besides the detected bursts also the brainstorm database is needed.
% Matthias Sure

%%
clear
% determine the relevant time points of the bursts
time_point={'start'; 'peak'; 'end'};
% Hemispheres of the STN
Hemispheres = {'right_STN','left_STN'};
% path to the bet bursts
DataPath='path_to_bursts';
% Medication sate
Med_State = 'OFF';
% cell array for max number of runs per patient
runs = {'run_1','run_2','run_3','run_4'};
% Subject ID
subjects = {'S001','S003','S004', ...};
%name the frequency bands
freqband_name = {'iBP'};% can be more than one but than you have to adap the frequency band defintion in the freqband variable in this script
% if you use individual beta peak frequencies: % beta peak frequencies --> 0 if no freq; shoud be identical as in the burst detction script
peak_freq = [21,15,16,17, ...; ...
             21,23,16,0, ...];
% set which bursts you want to work with         
thresh_prc = 75;
detection_identifier = 'first_z_bipolar';
Data_of_Interest = [detection_identifier '_' mat2str(thresh_prc)];
%% get the filenames
%files are from a brainstorm database
%in the first column you have the subject ID and in the second column the
%corresponding runs
switch Med_state
    case 'ON'
        % ON
        subject = {...
            'S001' {'S001_ON_run001', ...
            'S001_ON_run002', ...
            'S001_ON_run003'}; ...
            'S002' {'S002_ON_run001', ...
            'S002_ON_run002', ...
            'S002_ON_run003'}; 
            'S003' {'S003_ON_run001', ...
            'S003_ON_run002', ...
            'S003_ON_run003'}; 
            'S004' {'S003_ON_run001', ...
            'S004_ON_run002', ...
            'S004_ON_run003'}};
    case 'OFF'
        %OFF
        subject = {...
            'S001' {'S001_OFF_run001', ...
            'S001_OFF_run002', ...
            'S001_OFF_run003'}; ...
            'S002' {'S002_OFF_run001', ...
            'S002_OFF_run002', ...
            'S002_OFF_run003'}; 
            'S003' {'S003_OFF_run001', ...
            'S003_OFF_run002', ...
            'S003_OFF_run003'}; 
            'S004' {'S003_OFF_run001', ...
            'S004_OFF_run002', ...
            'S004_OFF_run003'}};
end

%% reorder burst information
database = ['link_to_datbase'];
for iSubj = 1 : size(subject,1)
for iFreq = 1 : size(freqband_name,2)
    % variable for the total time
    Time = 0;
    % variable for the total time without bad segments
    Time_wo_BAD = 0;
    % check if dataset exist
    try
        cd([database subjects{iSubj} '_PD_peri'])
    catch
        continue
    end
    for iRun = 1 : size(subject{iSubj,2},2)
        %% Evaluate measurement time 
        %Determine for every run the time of measurement and the time
        %segments of bad events. Moreover test if events are part of bad
        %segments
        % load the measured data
        a=dir(fullfile([database subject{iSubj,1} '/'  subject{iSubj,2}{iRun}  '/']));
        for i=1:length(a)
            if ~isempty(strfind(a(i).name, 'block'))
                load([database subject{iSubj,1} '/'  subject{iSubj,2}{iRun}  '/' a(i).name],'Events')
            end
        end
        Bad_times = [];
        Bad_LFP_times = [];
        for iEvent = 1 : size(Events,2)
            %determine the measurement time
            if strfind(Events(iEvent).label,'transient')
                duration = Events(iEvent).times(1,2)-Events(iEvent).times(1,1);
            %determine the Bad segments
            elseif strcmp(Events(iEvent).label,'BAD')
                Bad_times = Events(iEvent).times;
            %determine the BAD_LFP_segments
            elseif strcmp(Events(iEvent).label,'BAD_LFP')
                Bad_LFP_times = Events(iEvent).times;
            end
        end
        %determine the duration off all Bad_segments_together
        Bad_duration = 0;
        for iBad = 1 : size(Bad_times,2)
            Bad_duration = Bad_duration + Bad_times(2,iBad)-Bad_times(1,iBad);
        end
        %determine the duration off all Bad_LFP_segments_together
        Bad_LFP_duration = 0;
        for iBad = 1 : size(Bad_LFP_times,2)
            Bad_LFP_duration = Bad_LFP_duration + Bad_LFP_times(2,iBad)-Bad_LFP_times(1,iBad);
        end
        %sum up the measurment time of all runs
        Time = Time + duration;
        Time_wo_BAD = Time_wo_BAD + duration - (Bad_duration+Bad_LFP_duration);
                  
        for iChannel = 1 : 2
        % define the frequency band defintion
        freqband = [peak_freq(iChannel,iSubj)-3 peak_freq(iChannel,iSubj)+3];
        %% load new burst set
        try
            temp = load([DataPath Data_of_Interest '/' subjects{iSubj} '/Morlet_Bipolar_' Hemispheres{iChannel} '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subjects{iSubj} '_run_' num2str(iRun) '_' Data_of_Interest '_' Med_State '.mat']);
            %% check if bursts part of Bad segments
            % do it for every threshold
            DataMat.(Data_of_Interest{1}).events = temp.DataMat.F.events;        
            Amplitude_cell.(Data_of_Interest{1}) = temp.Amplitude_cell;
            % clear all reused arrays
            start = [];peak = [];peak_end = [];Burst_duration = [];
            %%remove the empty first line of DataMat
            if isempty(DataMat.(Data_of_Interest{1}).events(1).label)
                DataMat_new.(Data_of_Interest{1}).events = DataMat.(Data_of_Interest{1}).events(2:end);
            end
            iEEG = find(reshape(contains({DataMat_new.(Data_of_Interest{1}).events.label},[Hemispheres{iChannel} '_peak']),size(DataMat_new.(Data_of_Interest{1}).events)))+1;
            % pick the burst time points
            start = DataMat.(Data_of_Interest{1}).events(iEEG-1).times;
            peak = DataMat.(Data_of_Interest{1}).events(iEEG).times;
            peak_end = DataMat.(Data_of_Interest{1}).events(iEEG+1).times;
            % pick the amplitudes
            Burst_amplitude = Amplitude_cell.(Data_of_Interest{1}){iChannel};
            % get the burst duration
            Burst_duration = peak_end - start;
            %% group the burst duration and zScore amplitude of all runs togehter
            try
                Duration.(Data_of_Interest{1}){iChannel}(end+1:end+size(Burst_duration,2)) = Burst_duration;
            catch
                Duration.(Data_of_Interest{1}){iChannel} = Burst_duration;
            end
            try
                Amplitude.(Data_of_Interest{1}){iChannel}(end+1:end+size(Burst_duration,2)) = Burst_amplitude;
            catch
                Amplitude.(Data_of_Interest{1}){iChannel} = Burst_amplitude;
            end               
        catch
            Duration.(Data_of_Interest{1}){iChannel} = [];
            Amplitude.(Data_of_Interest{1}){iChannel} = [];                    
        end
        clear iEEG
        end
    end
    clear DataMat DataMat2 DataMat_noBad Beta_Power Beta_Power_noBad Amplitude_cell Amplitude_cell_noBad Thresholds 
    Thresholds = fieldnames(Power);
    for iThresh = 1 : size(Thresholds,1)
        Results.(Data_of_Interest{1}).zScore_Amplitude = Amplitude.(Data_of_Interest{1});
        Results.(Data_of_Interest{1}).Duration = Duration.(Data_of_Interest{1});
        Results.(Data_of_Interest{1}).MeasurementDuration = Time;%Time_wo_BAD;%
    end
    save([DataPath Data_of_Interest '/' subjects{iSubj} '/' 'Results_Bursts_' freqband_name{iFreq} '_' Data_of_Interest{1} '_' Med_State '.mat'],'Results','-v7.3')
    clear Thresholds Power Amplitude Duration Time_wo_BAD Time
end
end


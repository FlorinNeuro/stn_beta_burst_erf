% function Burst_detection_z_scored_raw_data_rereference_l_r_LFP_fun(isubject)
%% Burst detection with z-Score normalized data
%This script is for detecting beta bursts in LFP channels. The method is 
%adapted from Tinkhauser et al (2017). In this modification, the raw data 
%are normalized as the first z-score. The LFP channel should be from 
%directional DBS electrodes from Boston or St.Jude. A bipolar re ref will 
%be calculated. The script is divided into four parts. The first two parts 
%are used to define the burst detection threshold. This is determined using
%OFF and ON data. Here, the OFF data is prepared in the first block and the
%ON data in the second block. In the third block the bursts are determined 
%in OFF and in the fourth block in ON. The script is designed so that the 
%data is stored in a Brainstorm readable structure and also uses Brainstorm
%functions.

%Matthias Sure


% clear -except isubject
% isubject = 1;
warning('off','all')
%% Define how these detected bursts should be named
thresh_prc = 75; %burst detection threshold
detection_identifier = 'first_z_bipolar';
final_file_name = [detection_identifier '_' mat2str(thresh_prc)]; % create a string to later on identify your detected burst if you use i.e. different thresholds
%% Define the important paths
Med_state = 'OFF';
path_brainstorm_project = ['.../brainstorm_db/STN_COR_' Med_state '/'];
path_brainstorm_code = '.../brainstorm3/';
addpath(genpath(path_brainstorm_code));
OutputPath='Where to store the bursts'; % in this folder should be a folder for each subject with the subject ID as folder name
BrainstormDbDir = '.../brainstorm_db';
path_brainsotrm_data=['.../brainstorm_db/STN_COR_' Med_state '/data/'];
% Name the frequency bands for burst detection
freqband_name = {'Morlet_beta_peak'};% can be more than one but than you have to adap the frequency band defintion in the four blocks

subject_ID = {'S001','S003','S004', ...};
% determine for each subject the used electrode: St.Jude SJ, Boston Scientific BS   
Electrode = {'SJ','SJ','BS','SJ', ...};
% determine for each subject the beta frequency --- each column one
% subject; first row right bipolar LFP; second row left bipolar LFP; enter 0 if no peak frequency could be find 
peak_freq = [21,15,16,17, ...; ...
             21,23,16,0, ...];
% Name the max number of runs per subject
Runs = {'Run1','Run2','Run3','Run4'};
% starting, peak and end point of burst
time_point={'start'; 'peak'; 'end'};
%create a matrix for the thersholds, for each subject, frequency band and
%channel
Subj_Thresholds = zeros(size(freqband_name,2),size(subject_ID,2),2);

Hemispheres = {'right_STN','left_STN'};


%%how to start brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
ProtocolName = ['STN_COR_' Med_state];
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

%% Determine bursts threshold
%OFF Part
Med_state = 'OFF'; % if you use a different Medication states
%get the filenames
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

for irun=1:length(subject{isubject,2})
    %% Determine new Re-Reference Signal
    % load the Channel information
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/channel_vectorview306_acc1.mat'], 'Channel')
    % load the data (Channels x Time), Time vector, and the information
    % which channels are good or bad
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/data_block001.mat'], 'F','Time','ChannelFlag')
    % get the position of the LFP Channels in the data matrix F
    SEEG_index = zeros(16,1);% index for the LFP=SEEG channel
    nn = 0;
    for n = 1 : size(Channel,2)
        if strcmp(Channel(n).Type,'SEEG')
            nn = nn +1;
            SEEG_index(nn) = n;
        end
    end
    % re ref the LFP signal based on the used electrode. It should top
    % vs bottom but it also possible to re ref top vs the three
    % directional contacts at the bottom or bottom vs the three
    % contacts at the top. If one of the three contacts is bad the re
    % ref channel is bad.
    LFP_data = zeros(2,size(F,2));
    if strcmp(Electrode{isubject},'SJ')
        %%right STN
        if ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(5)) == 1
            temp_LFP = F(SEEG_index(5),:) - F(SEEG_index(1),:); 
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(5)) == 1 && ChannelFlag(SEEG_index(2)) == 1 && ChannelFlag(SEEG_index(4)) == 1 && ChannelFlag(SEEG_index(7)) == 1
            temp_LFP = F(SEEG_index(5),:) - mean([F(SEEG_index(2),:);F(SEEG_index(4),:);F(SEEG_index(7),:)]);
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(3)) == 1 && ChannelFlag(SEEG_index(6)) == 1 && ChannelFlag(SEEG_index(8)) == 1
            temp_LFP = mean([F(SEEG_index(3),:);F(SEEG_index(6),:);F(SEEG_index(8),:)]) - F(SEEG_index(1),:);
            LFP_data(1,:) = temp_LFP;
        end
        %%left STN
        if ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(13)) == 1
            temp_LFP = F(SEEG_index(13),:) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(13)) == 1 && ChannelFlag(SEEG_index(10)) == 1 && ChannelFlag(SEEG_index(12)) == 1 && ChannelFlag(SEEG_index(15)) == 1
            temp_LFP = F(SEEG_index(13),:) - mean([F(SEEG_index(10),:);F(SEEG_index(12),:);F(SEEG_index(15),:)]);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(11)) == 1 && ChannelFlag(SEEG_index(14)) == 1 && ChannelFlag(SEEG_index(16)) == 1
            temp_LFP = mean([F(SEEG_index(11),:);F(SEEG_index(14),:);F(SEEG_index(16),:)]) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        end
    elseif strcmp(Electrode{isubject},'BS')
        %%right STN
        if ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(8)) == 1
            temp_LFP = F(SEEG_index(8),:) - F(SEEG_index(1),:); 
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(8)) == 1 && ChannelFlag(SEEG_index(2)) == 1 && ChannelFlag(SEEG_index(3)) == 1 && ChannelFlag(SEEG_index(4)) == 1
            temp_LFP = F(SEEG_index(8),:) - mean([F(SEEG_index(2),:);F(SEEG_index(3),:);F(SEEG_index(4),:)]);
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(5)) == 1 && ChannelFlag(SEEG_index(6)) == 1 && ChannelFlag(SEEG_index(7)) == 1
            temp_LFP = mean([F(SEEG_index(5),:);F(SEEG_index(6),:);F(SEEG_index(7),:)]) - F(SEEG_index(1),:);
            LFP_data(1,:) = temp_LFP;
        end
        %%left STN
        if ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(16)) == 1
            temp_LFP = F(SEEG_index(16),:) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(16)) == 1 && ChannelFlag(SEEG_index(10)) == 1 && ChannelFlag(SEEG_index(11)) == 1 && ChannelFlag(SEEG_index(12)) == 1
            temp_LFP = F(SEEG_index(16),:) - mean([F(SEEG_index(10),:);F(SEEG_index(11),:);F(SEEG_index(12),:)]);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(13)) == 1 && ChannelFlag(SEEG_index(14)) == 1 && ChannelFlag(SEEG_index(15)) == 1
            temp_LFP = mean([F(SEEG_index(13),:);F(SEEG_index(14),:);F(SEEG_index(15),:)]) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        end            
    end     
    %% create Data Strucutre for brainstorm
    data.Atlas.Name = 'process_extract_scout';
    data.Atlas.Scouts(1).Vertices = [1 2 3];
    data.Atlas.Scouts(1).Seed = [2];
    data.Atlas.Scouts(1).Color = [0.5 0.5 0.5];
    data.Atlas.Scouts(1).Label = 'trash1';
    data.Atlas.Scouts(1).Function = 'pca';
    data.Atlas.Scouts(1).Region = 'LT';
    data.Atlas.Scouts(1).Handles = [];
    data.Atlas.Scouts(2).Vertices = [5 6 7];
    data.Atlas.Scouts(2).Seed = [6];
    data.Atlas.Scouts(2).Color = [0.5 0.5 0.5];
    data.Atlas.Scouts(2).Label = 'trash2';
    data.Atlas.Scouts(2).Function = 'pca';
    data.Atlas.Scouts(2).Region = 'LT';
    data.Atlas.Scouts(2).Handles = [];        
    data.ChannelFlag = ChannelFlag;
    data.Comment = 'trash';
    data.Description = {'right Re-Ref LFP';'left Re-Ref LFP'};
    data.DisplayUnits = [];
    data.Events = struct([]);
    data.History{1,1} = date;
    data.History{1,2} = 'project';
    data.History{1,3} = 'create test template';
    data.Leff = 1;
    data.nAvg = 1;
    data.Std = [];
    data.SurfaceFile = [subject{isubject,1} 'tess_cortex_pial_low.mat'];
    data.Time = Time;
    data.Value = LFP_data;
    data_name_OFF = [path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/matrix_scout_201211_0000.mat'];
    save(data_name_OFF,'-struct','data')
    % reload the bs database
    [sSubject, iSubject] = bst_get('Subject', subject{isubject,1});
    db_reload_conditions(iSubject);

    % create the filenames to run the brainstorm functions
    sFiles{1} = [subject{isubject,1} '/' subject{isubject,2}{irun} '/matrix_scout_201211_0000.mat'];
    sFiles_int = sFiles;
        
    %% pre process the LFP signal for burst detection
    for iLFP = 1 : 2
        % define the frequency band defintion
        freqband = [peak_freq(iLFP,isubject)-3 peak_freq(iLFP,isubject)+3];
        % z-score normalise the data
        sFiles_z = bst_process('CallProcess', 'process_baseline_norm', sFiles_int, [], ...
            'baseline',  [], ...
            'method',    'zscore', ...  % Z-score transformation:    x_std = (x - &mu;) / &sigma;
            'overwrite', 0);  
        for iFreq = 1 : size(freqband_name,2) % if you want to use more than one frequency band definiton
            % if you want beta peak bursts but no beta peak --> no need to detect bursts.
            if peak_freq(iLFP,isubject) == 0 && iFreq == 1
                continue
            end
            %% Process: Time-frequency (Morlet wavelets)
            sFiles = bst_process('CallProcess', 'process_timefreq', sFiles_z, [], ...
                'sensortypes', 'SEEG', ...
                'edit',        struct(...
                     'Comment',         'Power,FreqBands', ...
                     'TimeBands',       [], ...
                     'Freqs',           {{'Freq', [num2str(freqband(iFreq,1)) ', ' num2str(freqband(iFreq,2))], 'mean'}}, ...
                     'MorletFc',        1, ...
                     'MorletFwhmTc',    10, ...
                     'ClusterFuncTime', 'none', ...
                     'Measure',         'power', ...
                     'Output',          'all', ...
                     'SaveKernel',      0), ...
                'normalize',   'multiply');  % 1/f compensation: Multiply output values by frequency
            %load the Time-frequency signal
            load ([path_brainsotrm_data sFiles.FileName],'TF');
            %create an empty entry for the preocesed data
            DC_Data_OFF.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP}) = [];
            % get the TF data from the channel of interest
            Morlet_data = TF(iLFP,:,1);
            %% temporal smoothing with 0.2s
            % sampling rate
            sRate = 1000;
            Smooth_time = 0.2;%s Smoothing time Window
            Smooth_timepoints = sRate*Smooth_time;
            Smooth_data = smooth(Morlet_data,Smooth_timepoints, 'moving')';
            %% DC-Correction with a time constant of 20 s  
            DC_time = 20;%s DC-correction time constant
            DC_timepoints = sRate*DC_time;
            DC_data = zeros(size(Smooth_data));
            if size(Smooth_data,2) > DC_timepoints
                for iData = 1 : size(Smooth_data,2)
                    if iData <= DC_timepoints/2
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(1:iData+DC_timepoints/2));
                    elseif iData >= size(Smooth_data,2) - DC_timepoints/2
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(iData-DC_timepoints/2:size(Smooth_data,2)));
                    else
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(iData-DC_timepoints/2:iData+DC_timepoints/2));
                    end
                end
            end
            % save pre-processed data for later use
            DC_Data_OFF.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP}) = DC_data;
            % determine Trshold based on OFF and ON data
            if isempty(Subj_Thresholds(iFreq,isubject,iLFP)) || Subj_Thresholds(iFreq,isubject,iLFP) == 0
                Subj_Thresholds(iFreq,isubject,iLFP) =  prctile(DC_data, thresh_prc);
            else
                Subj_Thresholds(iFreq,isubject,iLFP) = mean([Subj_Thresholds(iFreq,isubject,iLFP) prctile(DC_data, thresh_prc)]);
            end                
        end
        % Process: Delete selected files
        sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
        'target', 1);  % Delete selected files  
    end
    % Process: Delete selected files
    sFiles = bst_process('CallProcess', 'process_delete', sFiles_z, [], ...
    'target', 1);  % Delete selected files   
end
brainstorm stop   

%% ON part
%% Define the important paths
Med_state = 'ON';
path_brainstorm_project = ['.../brainstorm_db/STN_COR_' Med_state '/'];
path_brainsotrm_data=['.../brainstorm_db/STN_COR_' Med_state '/data/'];


%%how to start brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
ProtocolName = ['STN_COR_' Med_state];
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

%get the filenames
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


for irun=1:length(subject{isubject,2})
    %% Determine new Re-Reference Signal
    % load the Channel information
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/channel_vectorview306_acc1.mat'], 'Channel')
    % load the data (Channels x Time), Time vector, and the information
    % which channels are good or bad
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/data_block001.mat'], 'F','Time','ChannelFlag')
    % get the position of the LFP Channels in the data matrix F
    SEEG_index = zeros(16,1);% index for the LFP=SEEG channel
    nn = 0;
    for n = 1 : size(Channel,2)
        if strcmp(Channel(n).Type,'SEEG')
            nn = nn +1;
            SEEG_index(nn) = n;
        end
    end
    % re ref the LFP signal based on the used electrode. It should top
    % vs bottom but it also possible to re ref top vs the three
    % directional contacts at the bottom or bottom vs the three
    % contacts at the top. If one of the three contacts is bad the re
    % ref channel is bad.
    LFP_data = zeros(2,size(F,2));
    if strcmp(Electrode{isubject},'SJ')
        %%right STN
        if ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(5)) == 1
            temp_LFP = F(SEEG_index(5),:) - F(SEEG_index(1),:); 
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(5)) == 1 && ChannelFlag(SEEG_index(2)) == 1 && ChannelFlag(SEEG_index(4)) == 1 && ChannelFlag(SEEG_index(7)) == 1
            temp_LFP = F(SEEG_index(5),:) - mean([F(SEEG_index(2),:);F(SEEG_index(4),:);F(SEEG_index(7),:)]);
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(3)) == 1 && ChannelFlag(SEEG_index(6)) == 1 && ChannelFlag(SEEG_index(8)) == 1
            temp_LFP = mean([F(SEEG_index(3),:);F(SEEG_index(6),:);F(SEEG_index(8),:)]) - F(SEEG_index(1),:);
            LFP_data(1,:) = temp_LFP;
        end
        %%left STN
        if ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(13)) == 1
            temp_LFP = F(SEEG_index(13),:) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(13)) == 1 && ChannelFlag(SEEG_index(10)) == 1 && ChannelFlag(SEEG_index(12)) == 1 && ChannelFlag(SEEG_index(15)) == 1
            temp_LFP = F(SEEG_index(13),:) - mean([F(SEEG_index(10),:);F(SEEG_index(12),:);F(SEEG_index(15),:)]);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(11)) == 1 && ChannelFlag(SEEG_index(14)) == 1 && ChannelFlag(SEEG_index(16)) == 1
            temp_LFP = mean([F(SEEG_index(11),:);F(SEEG_index(14),:);F(SEEG_index(16),:)]) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        end
    elseif strcmp(Electrode{isubject},'BS')
        %%right STN
        if ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(8)) == 1
            temp_LFP = F(SEEG_index(8),:) - F(SEEG_index(1),:); 
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(8)) == 1 && ChannelFlag(SEEG_index(2)) == 1 && ChannelFlag(SEEG_index(3)) == 1 && ChannelFlag(SEEG_index(4)) == 1
            temp_LFP = F(SEEG_index(8),:) - mean([F(SEEG_index(2),:);F(SEEG_index(3),:);F(SEEG_index(4),:)]);
            LFP_data(1,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(1)) == 1 && ChannelFlag(SEEG_index(5)) == 1 && ChannelFlag(SEEG_index(6)) == 1 && ChannelFlag(SEEG_index(7)) == 1
            temp_LFP = mean([F(SEEG_index(5),:);F(SEEG_index(6),:);F(SEEG_index(7),:)]) - F(SEEG_index(1),:);
            LFP_data(1,:) = temp_LFP;
        end
        %%left STN
        if ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(16)) == 1
            temp_LFP = F(SEEG_index(16),:) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(16)) == 1 && ChannelFlag(SEEG_index(10)) == 1 && ChannelFlag(SEEG_index(11)) == 1 && ChannelFlag(SEEG_index(12)) == 1
            temp_LFP = F(SEEG_index(16),:) - mean([F(SEEG_index(10),:);F(SEEG_index(11),:);F(SEEG_index(12),:)]);
            LFP_data(2,:) = temp_LFP;
        elseif ChannelFlag(SEEG_index(9)) == 1 && ChannelFlag(SEEG_index(13)) == 1 && ChannelFlag(SEEG_index(14)) == 1 && ChannelFlag(SEEG_index(15)) == 1
            temp_LFP = mean([F(SEEG_index(13),:);F(SEEG_index(14),:);F(SEEG_index(15),:)]) - F(SEEG_index(9),:);
            LFP_data(2,:) = temp_LFP;
        end            
    end     
    %% create Data Strucutre for brainstorm
    data.Atlas.Name = 'process_extract_scout';
    data.Atlas.Scouts(1).Vertices = [1 2 3];
    data.Atlas.Scouts(1).Seed = [2];
    data.Atlas.Scouts(1).Color = [0.5 0.5 0.5];
    data.Atlas.Scouts(1).Label = 'trash1';
    data.Atlas.Scouts(1).Function = 'pca';
    data.Atlas.Scouts(1).Region = 'LT';
    data.Atlas.Scouts(1).Handles = [];
    data.Atlas.Scouts(2).Vertices = [5 6 7];
    data.Atlas.Scouts(2).Seed = [6];
    data.Atlas.Scouts(2).Color = [0.5 0.5 0.5];
    data.Atlas.Scouts(2).Label = 'trash2';
    data.Atlas.Scouts(2).Function = 'pca';
    data.Atlas.Scouts(2).Region = 'LT';
    data.Atlas.Scouts(2).Handles = [];        
    data.ChannelFlag = ChannelFlag;
    data.Comment = 'trash';
    data.Description = {'right Re-Ref LFP';'left Re-Ref LFP'};
    data.DisplayUnits = [];
    data.Events = struct([]);
    data.History{1,1} = date;
    data.History{1,2} = 'project';
    data.History{1,3} = 'create test template';
    data.Leff = 1;
    data.nAvg = 1;
    data.Std = [];
    data.SurfaceFile = [subject{isubject,1} 'tess_cortex_pial_low.mat'];
    data.Time = Time;
    data.Value = LFP_data;
    data_name_ON = [path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/matrix_scout_201211_0000.mat'];
    save(data_name_ON,'-struct','data')
    % reload the bs database
    [sSubject, iSubject] = bst_get('Subject', subject{isubject,1});
    db_reload_conditions(iSubject);

    % create the filenames to run the brainstorm functions
    sFiles{1} = [subject{isubject,1} '/' subject{isubject,2}{irun} '/matrix_scout_201211_0000.mat'];
    sFiles_int = sFiles;
        
    %% pre process the LFP signal for burst detection
    for iLFP = 1 : 2
        % define the frequency band defintion
        freqband = [peak_freq(iLFP,isubject)-3 peak_freq(iLFP,isubject)+3];
        % z-score normalise the data
        sFiles_z = bst_process('CallProcess', 'process_baseline_norm', sFiles_int, [], ...
            'baseline',  [], ...
            'method',    'zscore', ...  % Z-score transformation:    x_std = (x - &mu;) / &sigma;
            'overwrite', 0);  
        for iFreq = 1 : size(freqband_name,2) % if you want to use more than one frequency band definiton
            % if you want beta peak bursts but no beta peak --> no need to detect bursts.
            if peak_freq(iLFP,isubject) == 0 && iFreq == 1
                continue
            end
            %% Process: Time-frequency (Morlet wavelets)
            sFiles = bst_process('CallProcess', 'process_timefreq', sFiles_z, [], ...
                'sensortypes', 'SEEG', ...
                'edit',        struct(...
                     'Comment',         'Power,FreqBands', ...
                     'TimeBands',       [], ...
                     'Freqs',           {{'Freq', [num2str(freqband(iFreq,1)) ', ' num2str(freqband(iFreq,2))], 'mean'}}, ...
                     'MorletFc',        1, ...
                     'MorletFwhmTc',    10, ...
                     'ClusterFuncTime', 'none', ...
                     'Measure',         'power', ...
                     'Output',          'all', ...
                     'SaveKernel',      0), ...
                'normalize',   'multiply');  % 1/f compensation: Multiply output values by frequency
            %load the Time-frequency signal
            load ([path_brainsotrm_data sFiles.FileName],'TF');
            %create an empty entry for the preocesed data
            DC_Data_ON.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP}) = [];
            % get the TF data from the channel of interest
            Morlet_data = TF(iLFP,:,1);
            %% temporal smoothing with 0.2s
            % sampling rate
            sRate = 1000;
            Smooth_time = 0.2;%s Smoothing time Window
            Smooth_timepoints = sRate*Smooth_time;
            Smooth_data = smooth(Morlet_data,Smooth_timepoints, 'moving')';
            %% DC-Correction with a time constant of 20 s  
            DC_time = 20;%s DC-correction time constant
            DC_timepoints = sRate*DC_time;
            DC_data = zeros(size(Smooth_data));
            if size(Smooth_data,2) > DC_timepoints
                for iData = 1 : size(Smooth_data,2)
                    if iData <= DC_timepoints/2
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(1:iData+DC_timepoints/2));
                    elseif iData >= size(Smooth_data,2) - DC_timepoints/2
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(iData-DC_timepoints/2:size(Smooth_data,2)));
                    else
                        DC_data(iData) = Smooth_data(iData) - mean(Smooth_data(iData-DC_timepoints/2:iData+DC_timepoints/2));
                    end
                end
            end
            % save pre-processed data for later use
            DC_Data_ON.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP}) = DC_data;
            % determine Trshold based on OFF and ON data
            if isempty(Subj_Thresholds(iFreq,isubject,iLFP)) || Subj_Thresholds(iFreq,isubject,iLFP) == 0
                Subj_Thresholds(iFreq,isubject,iLFP) =  prctile(DC_data, thresh_prc);
            else
                Subj_Thresholds(iFreq,isubject,iLFP) = mean([Subj_Thresholds(iFreq,isubject,iLFP) prctile(DC_data, thresh_prc)]);
            end                
        end
        % Process: Delete selected files
        sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
        'target', 1);  % Delete selected files  
    end
    % Process: Delete selected files
    sFiles = bst_process('CallProcess', 'process_delete', sFiles_z, [], ...
    'target', 1);  % Delete selected files   
end
brainstorm stop 
%% Determine Bursts OFF medication
Med_state = 'OFF';
%get the filenames
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
%% Define the important paths
path_brainstorm_project = ['.../brainstorm_db/STN_COR_' Med_state '/'];
for irun=1:length(subject{isubject,2})
    %load the Time vector for this LF
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/data_block001.mat'], 'Time')
    for iFreq = 1 : size(freqband_name,2)
        for iLFP = 1 : 2
        % create a Data struct to later save the timepoinrts of the burst
        % could be also used directly in brainstorm
        clear DataMat
        DataMat.F.events.label = [];
        DataMat.F.events.color = [];
        DataMat.F.events.epochs = [];
        DataMat.F.events.samples = [];
        DataMat.F.events.times = [];
        DataMat.F.events.reactTimes = [];
        DataMat.F.events.select = [];
        % save the burst amplitude
        Amplitude_cell = {};               
            %define the frequency band
            freqband = [peak_freq(iLFP,isubject)-3 peak_freq(iLFP,isubject)+3];
            % check if data was preprocessed or the re ref was not possible
            try
                DC_data = DC_Data_OFF.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP});
            catch
                continue
            end
            % get the timepoints above threshold
            above_thresh = (DC_data > Subj_Thresholds(iFreq,isubject,iLFP));
            % get the burst which a at least a defined minimal duration
            min_time = 80/2;%ms
            above_min_thresh = zeros(size(above_thresh));
            for n = min_time+1 : size(above_thresh,2)-min_time
                if sum(above_thresh(n-min_time:n+min_time)) >= min_time*2
                    above_min_thresh(n) = 1; 
                end
            end
            above_min_thresh_pos = find(above_min_thresh == 1);
            for n = 1 : size(above_min_thresh_pos,2)
            above_min_thresh(above_min_thresh_pos(n)-min_time:above_min_thresh_pos(n)+min_time) = 1;
            end
            above_min_thresh(1) = 0;
            above_min_thresh(end) = 0;
            % get the start, end, peak timepoint of the burst and als the
            % peak value
            start_beta= [];
            end_beta = [];
            peak_beta = [];
            peak_beta_value = [];
            for n = 1 : size(above_min_thresh,2)-1
                if above_min_thresh(n)-above_min_thresh(n+1) == -1
                    start_beta(end+1) = (n+1);
                end            
            end
            for n = 2 : size(above_min_thresh,2)
                if above_min_thresh(n-1)-above_min_thresh(n) == 1
                    end_beta(end+1) = (n-1);
                end            
            end
            for n = 1 : size(start_beta,2)
                [peak_beta_value(end+1),peak_beta(end+1)] = max(DC_data(start_beta(n):end_beta(n))); 
                peak_beta(n) = peak_beta(n)+start_beta(n);
            end
            % save the beta time points in a struct
            for n = 1 : size(time_point,1)
                 DataMat.F.events(end+1).label = ['Morlet_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' Hemispheres{iLFP} '_' time_point{n}];
                 DataMat.F.events(end).color = [0 0 1];
                 DataMat.F.events(end).epochs = ones(size(peak_beta));
                 if n == 1
                     DataMat.F.events(end).samples = start_beta;
                     DataMat.F.events(end).times = Time(start_beta);
                 elseif n == 2
                     DataMat.F.events(end).samples = peak_beta;
                     DataMat.F.events(end).times = Time(peak_beta);                        
                 elseif n == 3
                     DataMat.F.events(end).samples = end_beta;
                     DataMat.F.events(end).times = Time(end_beta);                        
                 end
                 DataMat.F.events(end).reactTimes = [];
                 DataMat.F.events(end).select = 1;  
            end
            Amplitude_cell{iLFP} = DC_data(peak_beta);
            save([OutputPath subject_ID{isubject}(1:4)  '/Morlet_Bipolar_' Hemispheres{iLFP} '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{isubject} '_run_' num2str(irun) '_' final_file_name '_' Med_state],'DataMat','Amplitude_cell')
        end
    end
end


%% Define the important paths
Med_state = 'ON';
%get the filenames
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
%% Define the important paths
path_brainstorm_project = ['.../brainstorm_db/STN_COR_' Med_state '/'];
for irun=1:length(subject{isubject,2})
    %load the Time vector for this LF
    load([path_brainstorm_project 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/data_block001.mat'], 'Time')
    for iFreq = 1 : size(freqband_name,2)
        for iLFP = 1 : 2
            % create a Data struct to later save the timepoinrts of the burst
            % could be also used directly in brainstorm
            clear DataMat
            DataMat.F.events.label = [];
            DataMat.F.events.color = [];
            DataMat.F.events.epochs = [];
            DataMat.F.events.samples = [];
            DataMat.F.events.times = [];
            DataMat.F.events.reactTimes = [];
            DataMat.F.events.select = [];
            % save the burst amplitude
            Amplitude_cell = {};   
            %define the frequency band
            freqband = [peak_freq(iLFP,isubject)-3 peak_freq(iLFP,isubject)+3];
            % check if data was preprocessed or the re ref was not possible
            try
                DC_data = DC_Data_ON.(subject_ID{isubject}).(Runs{irun}).(freqband_name{iFreq}).(Hemispheres{iLFP});
            catch
                continue
            end
            % get the timepoints above threshold
            above_thresh = (DC_data > Subj_Thresholds(iFreq,isubject,iLFP));
            % get the burst which a at least a defined minimal duration
            min_time = 80/2;%ms
            above_min_thresh = zeros(size(above_thresh));
            for n = min_time+1 : size(above_thresh,2)-min_time
                if sum(above_thresh(n-min_time:n+min_time)) >= min_time*2
                    above_min_thresh(n) = 1; 
                end
            end
            above_min_thresh_pos = find(above_min_thresh == 1);
            for n = 1 : size(above_min_thresh_pos,2)
            above_min_thresh(above_min_thresh_pos(n)-min_time:above_min_thresh_pos(n)+min_time) = 1;
            end
            above_min_thresh(1) = 0;
            above_min_thresh(end) = 0;
            % get the start, end, peak timepoint of the burst and als the
            % peak value
            start_beta= [];
            end_beta = [];
            peak_beta = [];
            peak_beta_value = [];
            for n = 1 : size(above_min_thresh,2)-1
                if above_min_thresh(n)-above_min_thresh(n+1) == -1
                    start_beta(end+1) = (n+1);
                end            
            end
            for n = 2 : size(above_min_thresh,2)
                if above_min_thresh(n-1)-above_min_thresh(n) == 1
                    end_beta(end+1) = (n-1);
                end            
            end
            for n = 1 : size(start_beta,2)
                [peak_beta_value(end+1),peak_beta(end+1)] = max(DC_data(start_beta(n):end_beta(n))); 
                peak_beta(n) = peak_beta(n)+start_beta(n);
            end
            % save the beta time points in a struct
            for n = 1 : size(time_point,1)
                 DataMat.F.events(end+1).label = ['Morlet_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' Hemispheres{iLFP} '_' time_point{n}];
                 DataMat.F.events(end).color = [0 0 1];
                 DataMat.F.events(end).epochs = ones(size(peak_beta));
                 if n == 1
                     DataMat.F.events(end).samples = start_beta;
                     DataMat.F.events(end).times = Time(start_beta);
                 elseif n == 2
                     DataMat.F.events(end).samples = peak_beta;
                     DataMat.F.events(end).times = Time(peak_beta);                        
                 elseif n == 3
                     DataMat.F.events(end).samples = end_beta;
                     DataMat.F.events(end).times = Time(end_beta);                        
                 end
                 DataMat.F.events(end).reactTimes = [];
                 DataMat.F.events(end).select = 1;  
            end
            Amplitude_cell{iLFP} = DC_data(peak_beta);
            save([OutputPath subject_ID{isubject}(1:4)  '/Morlet_Bipolar_' Hemispheres{iLFP} '_' num2str(freqband(iFreq,1)) '_' num2str(freqband(iFreq,2)) '_' subject_ID{isubject} '_run_' num2str(irun) '_' final_file_name '_' Med_state],'DataMat','Amplitude_cell')
        end
    end
end
    
% delete the created files    
try
    delete(data_name_OFF)
catch
end
try
    delete(data_name_ON)
catch
end

  


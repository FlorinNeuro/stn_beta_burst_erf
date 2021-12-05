%% count the number of bursts
%This script is used to count number of burst after bipolar reref. 
%Restrictions can be made on the used, conditions, patients and hemispheres.
%It is also possible to make comparisons between groups.
%Matthias Sure

clear
% Subject ID
subjects = {'S001','S003','S004', ...};
    %name the frequency bands
freqband_name = {'iBP'};% can be more than one but than you have to adap the frequency band defintion in the freqband variable in this script
% if you use individual beta peak frequencies: % beta peak frequencies --> 0 if no freq; shoud be identical as in the burst detction script
peak_freq = [21,15,16,17, ...; ...
             21,23,16,0, ...];
% Hemispheres of the STN
Hemispheres = {'right_STN','left_STN'};
% set which bursts you want to work with         
thresh_prc = 75;
detection_identifier = 'first_z_bipolar';
Data_of_Interest = [detection_identifier '_' mat2str(thresh_prc)];
% Medication states
med_state = {'OFF','ON'};
%path to the burst files
data_path = ['path_to_burst_data'];

for iSubject = 1 : size(subjects,2)
    for iSide = 1 : size(Hemispheres,2)
        for iFreq = 1 : size(freqband_name,2)
            for iMed = 1 : size(med_state,2)
                try
                    load([data_path subjects{iSubject} '\Results_Bursts_' freqband_name{iFreq} '_' Data_of_Interest '_' med_state{iMed} '.mat'])
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time = ... 
                        Results.(Data_of_Interest).MeasurementDuration;
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration = ... 
                        Results.(Data_of_Interest).Duration;
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power = ... 
                        Results.(Data_of_Interest).Power;
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).zScore_Amplitude = ... 
                        Results.(Data_of_Interest).zScore_Amplitude;
                    clear Results                    
                catch
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time = ... 
                        [];
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration = ... 
                        cell(1,2);
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power = ... 
                        cell(1,2);
                    Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).zScore_Amplitude = ...
                        cell(1,2);                        
                end
            end
        end
    end
end
%check for what burst duration time intervalls you want to count the burst
time_intervalls = [0.08 : 0.02 : 6];
time_intervalls = 0.2;
counter = 0; % for how many conditions/subjects burst were detected
counter_l = 0; % for how many conditions/subjects burst were detected only left STN
counter_r = 0; % for how many conditions/subjects burst were detected only right STN
Collect_Duration = []; %all durations
hist_mat = [];% matrix with the number of bursts and the subject number


for iSet = 1
    for iFreq = 1 : size(freqband_name,2)
        for iMed = 1 %: size(med_state,2)
            for iChan = 1 : size(Hemispheres,2)
                    for iSubject = 1 : size(subjects,2)
                        if isempty(Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan})
                            Data_Density(iSubject,iChan,iMed,iFreq,iSet) = NaN;
                            Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = NaN;
                            Data_Power(iSubject,iChan,iMed,iFreq,iSet) = NaN;
%                         elseif iChan == 2 && iSubject == 1
%                             Data_Density(iSubject,iChan,iMed,iFreq,iSet) = NaN;
%                             Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = NaN;
%                             Data_Power(iSubject,iChan,iMed,iFreq,iSet) = NaN;                            
                        else
                            Data_Density(iSubject,iChan,iMed,iFreq,iSet) = size(Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan},2)/(Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time/60);
                            Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = median(Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan});
                            Data_Power(iSubject,iChan,iMed,iFreq,iSet) = mean(Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power{iChan});
                            if iFreq == 1 && iMed == 1 %&& iChan == 2 %&& iSubject ~= 1 % if you only want to count bursts for specific conditions
                                if iChan == 2 && iSubject == 1 %% if you want to skip specific datasets
                                    continue
                                end
                                counter = counter + 1;
                                temp_Duration = Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan};
                                for iTime = 1 : size(time_intervalls,2)
%                                     hist_mat(counter,iTime) = sum(temp_Duration > time_intervalls(iTime) & temp_Duration < time_intervalls(iTime)+0.8);
                                    hist_mat(counter,iTime) = sum(temp_Duration > time_intervalls(iTime));
                                    hist_mat(counter,2) = iSubject;
                                    
%                                     hist_mat(counter,iTime) = sum(temp_Duration > 100.0);
                                end
                                Number_all_bursts(counter,1) = length(temp_Duration);
                                Collect_Duration = [Collect_Duration Data.(Data_of_Interest).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan}];
                                if iChan == 1
                                    counter_r = counter_r + 1;
                                    Number_r_bursts(counter_r,1) = length(temp_Duration);
                                elseif iChan == 2
                                    counter_l = counter_l + 1;
                                    Number_l_bursts(counter_l,1) = length(temp_Duration);
                                end
                            end
                        end                        
                    end
            end
        end
    end
end
hist(Collect_Duration)
bar(time_intervalls,sum(hist_mat));
[mean_max_bursts,max_time_slot] = max(mean(hist_mat));
max_bursts = hist_mat(:,max_time_slot);
relative_max_bursts = max_bursts./Number_all_bursts;
relative_all_bursts = sum(max_bursts)./sum(Number_all_bursts);
[h,p,ci,stats] = ttest2(Number_r_bursts,Number_l_bursts)
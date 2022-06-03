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
        for iFreq = 1 %: size(freqband_name,2)
            for iMed = 1 : size(med_state,2)
                % get the Time, Duration, Power and/or z-Sore Amplitude. If
                % data was not recorded for this patient add empty cells
                try
                    load([data_path subjects{iSubject} '\Results_Bursts_' freqband_name{iFreq} '_' final_file_name '_' med_state{iMed} '.mat'])
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time = ... 
                        Results.(final_file_name).MeasurementDuration;
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration = ... 
                        Results.(final_file_name).Duration;
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power = ... 
                        Results.(final_file_name).Power;
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).zScore_Amplitude = ... 
                        Results.(final_file_name).zScore_Amplitude;
                    clear Results                    
                catch
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time = ... 
                        [];
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration = ... 
                        cell(1,2);
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power = ... 
                        cell(1,2);
                    Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).zScore_Amplitude = ...
                        cell(1,2);                        
                end
            end
        end
    end
end
%investiget different intervalls for the duration
duration_intervalls = [0.1 : 0.1 : 1];
% time_intervalls = 0.2;
% count the number of valid datasets, left, right and both hemispheres
counter = 0;
counter_l = 0;
counter_r = 0;
%collect all burst duration values in one array
Collect_Duration = [];
%create a histogram for the burst duration
hist_mat = [];

%select wich data should be investigated
iFreq = 1;
iMed = 1;
iSet = 1;

for iChan = 1 : size(Hemispheres,2)
    for iSubject = 1 : size(subjects,2)
        if isempty(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan})
            Data_Density(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_short(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_long(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_data(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Power(iSubject,iChan,iMed,iFreq,iSet) = NaN;
        elseif iChan == 2 && iSubject == 1
            Data_Density(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_short(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_long(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Density_data(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = NaN;
            Data_Power(iSubject,iChan,iMed,iFreq,iSet) = NaN;                            
        else
            Data_Density(iSubject,iChan,iMed,iFreq,iSet) = size(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan},2)/(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time/60);
            Data_Density_short(iSubject,iChan,iMed,iFreq,iSet) = size(find(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan} < 0.3),2)/(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time/60);
            Data_Density_long(iSubject,iChan,iMed,iFreq,iSet) = size(find(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan} > 0.5),2)/(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time/60);
            Data_Density_data(iSubject,iChan,iMed,iFreq,iSet) = size(find(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan} > 0.2),2)/(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Recording_Time/60);
            Data_Duration(iSubject,iChan,iMed,iFreq,iSet) = median(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan});
            Data_Power(iSubject,iChan,iMed,iFreq,iSet) = mean(Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Power{iChan});
            if iFreq == 1 && iMed == 1 %&& iChan == 2 %&& iSubject ~= 1
                if iChan == 2 && iSubject == 1
                    continue
                end
                counter = counter + 1;
                temp_Duration = Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan};
                hist_mat(counter,1) = iSubject;
                hist_mat(counter,2) = sum(temp_Duration < duration_intervalls(1));
                for iTime = 1 : size(duration_intervalls,2)-1
                      hist_mat(counter,iTime+2) = sum(temp_Duration > duration_intervalls(iTime) & temp_Duration < duration_intervalls(iTime+1));
                end
                hist_mat(counter,size(duration_intervalls,2)+2) = sum(temp_Duration > duration_intervalls(size(duration_intervalls,2)));
                %count number of all burts
                Number_all_bursts(counter,1) = length(temp_Duration);
                Collect_Duration = [Collect_Duration Data.(final_file_name).(freqband_name{iFreq}).(med_state{iMed}).(subjects{iSubject}).Duration{iChan}];
                if iChan == 1
                    counter_r = counter_r + 1;
                    Number_r_bursts(counter_r,1) = length(temp_Duration);
                    hist_mat_r(counter_r,1) = iSubject;
                    hist_mat_r(counter_r,2) = sum(temp_Duration < duration_intervalls(1));
                    hist_mat_r_perc(counter_r,1) = iSubject;
                    hist_mat_r_perc(counter_r,2) = sum(temp_Duration < duration_intervalls(1))/size(temp_Duration,2);                    
                    for iTime = 1 : size(duration_intervalls,2)-1
                          hist_mat_r(counter_r,iTime+2) = sum(temp_Duration > duration_intervalls(iTime) & temp_Duration < duration_intervalls(iTime+1));
                          hist_mat_r_perc(counter_r,iTime+2) = sum(temp_Duration > duration_intervalls(iTime) & temp_Duration < duration_intervalls(iTime+1))/size(temp_Duration,2);
                    end
                    hist_mat_r(counter_r,size(duration_intervalls,2)+2) = sum(temp_Duration > duration_intervalls(size(duration_intervalls,2))); 
                    hist_mat_r_perc(counter_r,size(duration_intervalls,2)+2) = sum(temp_Duration > duration_intervalls(size(duration_intervalls,2)))/size(temp_Duration,2); 
                elseif iChan == 2
                    counter_l = counter_l + 1;
                    Number_l_bursts(counter_l,1) = length(temp_Duration);
                    hist_mat_l(counter_l,1) = iSubject;
                    hist_mat_l(counter_l,2) = sum(temp_Duration < duration_intervalls(1));
                    hist_mat_l_perc(counter_l,1) = iSubject;
                    hist_mat_l_perc(counter_l,2) = sum(temp_Duration < duration_intervalls(1))/size(temp_Duration,2);
                    for iTime = 1 : size(duration_intervalls,2)-1
                          hist_mat_l(counter_l,iTime+2) = sum(temp_Duration > duration_intervalls(iTime) & temp_Duration < duration_intervalls(iTime+1));
                          hist_mat_l_perc(counter_l,iTime+2) = sum(temp_Duration > duration_intervalls(iTime) & temp_Duration < duration_intervalls(iTime+1))/size(temp_Duration,2);
                    end
                    hist_mat_l(counter_l,size(duration_intervalls,2)+2) = sum(temp_Duration > duration_intervalls(size(duration_intervalls,2)));
                    hist_mat_l_perc(counter_l,size(duration_intervalls,2)+2) = sum(temp_Duration > duration_intervalls(size(duration_intervalls,2)))/size(temp_Duration,2);
                end
            end
        end                        
    end
end
%% Bar plots
time_windows = categorical({'< 0.1','0.1 - 0.2','0.2 - 0.3','0.3 - 0.4','0.4 - 0.5','0.5 - 0.6','0.6 - 0.7','0.7 - 0.8','0.8 - 0.9','0.9 - 1.0','> 1.0'});
time_windows = reordercats(time_windows,{'< 0.1','0.1 - 0.2','0.2 - 0.3','0.3 - 0.4','0.4 - 0.5','0.5 - 0.6','0.6 - 0.7','0.7 - 0.8','0.8 - 0.9','0.9 - 1.0','> 1.0'});
figure
bar(time_windows,sum(hist_mat(:,2:end)));
xlabel('Time windows (s)')
ylabel('Number of bursts')
title('Temporal distribution of all bursts')

figure
subplot(1,2,1)
right = bar(time_windows,sum(hist_mat_r(:,2:end)));
xlabel('Time windows (s)')
ylabel('Number of bursts')
title('Temporal distribution of bursts from the right STN')
ylim([0 2200])
set(gca,'FontWeight','bold')
set(gca,'FontSize',12)
xtips = right(1).XEndPoints;
ytips = right(1).YEndPoints;
labels = string(right(1).YData);
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

% figure
subplot(1,2,2)
left = bar(time_windows,sum(hist_mat_l(:,2:end)));
xlabel('Time windows (s)')
title('Temporal distribution of bursts from the left STN')
ylim([0 2200])
set(gca,'FontWeight','bold')
set(gca,'FontSize',12)
xtips = left(1).XEndPoints;
ytips = left(1).YEndPoints;
labels = string(left(1).YData);
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

figure
both = bar(time_windows,[mean(hist_mat_r(:,2:end));mean(hist_mat_l(:,2:end))]);
xlabel('Time windows (s)')
ylabel('Number of bursts')
title('Mean temporal distribution of all STN beta bursts')
% ylim([0 2200])
set(gca,'FontWeight','bold')
set(gca,'FontSize',12)
xtips = both(1).XEndPoints;
ytips = both(1).YEndPoints;
labels = string(both(1).YData);
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips = both(2).XEndPoints;
ytips = both(2).YEndPoints;
labels = string(both(2).YData);
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

figure
data = [mean(hist_mat_r(:,2:end));mean(hist_mat_l(:,2:end))];
error = [std(hist_mat_r(:,2:end))/sqrt(size(hist_mat_r,1));std(hist_mat_l(:,2:end))/sqrt(size(hist_mat_l,1))];
b = bar(data');

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:ngroups
    x(:,i) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',data,error,'k','linestyle','none');
hold off
legend('right STN','left STN')
xlabel('Time windows (s)')
ylabel('Mean number of bursts per hemisphere')
title('Mean temporal distribution of all STN beta bursts')
% ylim([0 2200])
set(gca,'xticklabel',time_windows)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)


figure
data = [mean(hist_mat_r_perc(:,2:end)*100);mean(hist_mat_l_perc(:,2:end)*100)];
% error = [std(hist_mat_r_perc(:,2:end)*100)/sqrt(size(hist_mat_r_perc,1)*100);std(hist_mat_l_perc(:,2:end)*100)/sqrt(size(hist_mat_l_perc,1)*100)];
error = [std(hist_mat_r_perc(:,2:end)*100);std(hist_mat_l_perc(:,2:end)*100)];

b = bar(data');

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:ngroups
    x(:,i) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',data,error,'k','linestyle','none');
hold off
legend('right STN','left STN')
xlabel('Burst duration window [s]')
ylabel('Amount of bursts per duration window [%]')
title('Relative temporal distribution of all STN beta bursts')
% ylim([0 2200])
set(gca,'xticklabel',time_windows)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)
%% mean burst rate all
nanmean(Data_Density(:))
nanstd(Data_Density(:))
nanmean(Data_Density_data(:))
nanstd(Data_Density_data(:))
nanmean(Data_Density_short(:))
nanstd(Data_Density_short(:))
nanmean(Data_Density_long(:))
nanstd(Data_Density_long(:))
%% mean burst rate right
nanmean(Data_Density(:,1))
nanstd(Data_Density(:,1))
nanmean(Data_Density_data(:,1))
nanstd(Data_Density_data(:,1))
nanmean(Data_Density_short(:,1))
nanstd(Data_Density_short(:,1))
nanmean(Data_Density_long(:,1))
nanstd(Data_Density_long(:,1))
%% mean burst rate left
nanmean(Data_Density(:,2))
nanstd(Data_Density(:,2))
nanmean(Data_Density_data(:,2))
nanstd(Data_Density_data(:,2))
nanmean(Data_Density_short(:,2))
nanstd(Data_Density_short(:,2))
nanmean(Data_Density_long(:,2))
nanstd(Data_Density_long(:,2))

%% percentage of bursts
% over 200 ms
sum(sum(hist_mat(:,4:end)))/sum(sum(hist_mat(:,2:end)));
sum(sum(hist_mat_r(:,4:end)))/sum(sum(hist_mat_r(:,2:end)));
sum(sum(hist_mat_l(:,4:end)))/sum(sum(hist_mat_l(:,2:end)));
% over 500 ms
sum(sum(hist_mat(:,7:end)))/sum(sum(hist_mat(:,2:end)));
sum(sum(hist_mat_r(:,7:end)))/sum(sum(hist_mat_r(:,2:end)));
sum(sum(hist_mat_l(:,7:end)))/sum(sum(hist_mat_l(:,2:end)));
% under 300 ms
sum(sum(hist_mat(:,2:4)))/sum(sum(hist_mat(:,2:end)));
sum(sum(hist_mat_r(:,2:4)))/sum(sum(hist_mat_r(:,2:end)));
sum(sum(hist_mat_l(:,2:4)))/sum(sum(hist_mat_l(:,2:end)));

%% test for differences between hemisphers
%over 200 ms
[h,p,ci,stats] = ttest2(sum(hist_mat_r(:,4:end),2),sum(hist_mat_l(:,4:end),2))
%over 500 ms
[h,p,ci,stats] = ttest2(sum(hist_mat_r(:,7:end),2),sum(hist_mat_l(:,7:end),2))
%under 300 ms
[h,p,ci,stats] = ttest2(sum(hist_mat_r(:,2:4),2),sum(hist_mat_l(:,2:4),2))

[mean_max_bursts,max_time_slot] = max(mean(hist_mat));
max_bursts = hist_mat(:,max_time_slot);
relative_max_bursts = max_bursts./Number_all_bursts;
relative_all_bursts = sum(max_bursts)./sum(Number_all_bursts);
[h,p,ci,stats] = ttest2(Number_r_bursts,Number_l_bursts)
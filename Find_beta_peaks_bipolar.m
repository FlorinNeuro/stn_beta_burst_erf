%% Dieses Skript soll Peaks in den Powert spektren finden
% close all
clear
Med_states = {'OFF','ON'};
chanName='SEEG';%detect Events in the LFP channels
subjects = {'S001','S003','S004','S008','S009','S010','S011','S012','S013','S014','S016','S020','S021','S022','S023','S025','S027','S028','S029','S032','S033','S034','S036','S037','S038','S041','S043'};
data_path = 'E:\brainstorm_db\test\data\Test\';
beta_peak_freq_l = NaN(4,length(subjects));
beta_peak_amp_l = NaN(4,length(subjects));
beta_peak_freq_r = NaN(4,length(subjects));
beta_peak_amp_r = NaN(4,length(subjects));
peak_freqs = [21,15,16,17,24,18,24,20, 0, 0, 0,17,23,27, 0, 0,0, 0,0,0,15,16,15,20,14,15,30; ...
             21,23,16,15,24,21,24,20,23,27,23,17, 0,27,22,20,0,27,0,0,15, 0, 0,18, 0, 0,15];
%%
for iSubj = 1 : size(subjects,2)-1
    try
        cd([data_path subjects{iSubj}])
    catch
        continue
    end
    files = dir;
    n_count = 0;
    for iFile = 1 : size(files,1)
        if contains(files(iFile).name,'multiply')
            n_count = n_count + 1;
            load([data_path subjects{iSubj} '/' files(iFile).name],'TF','Freqs');%load the names
%             Freqs = round(Freqs);
            %get the 1Hz and 45 Hz Frequency 
            Freq_index_1 = find(Freqs == 1);
            Freq_index_45 = find(Freqs == 45);
            Freq_index_55 = find(Freqs == 55);
            Freq_index_95 = find(Freqs == 95);
%             Freqs = Freqs(Freq_index_1:Freq_index_45);
            %rigth LFP
            if sum(TF(1,1,:)) ~= 0
%                 Power = squeeze(TF(1,1,Freq_index_1:Freq_index_45));%/(sum(squeeze(TF(1,1,Freq_index_1:Freq_index_45)))+sum(squeeze(TF(1,1,Freq_index_55:Freq_index_95))));
                Power = squeeze(TF(1,1,:));%Freq_index_1:Freq_index_45));%/(sum(squeeze(TF(1,1,Freq_index_1:Freq_index_45)))+sum(squeeze(TF(1,1,Freq_index_55:Freq_index_95))));
%                 [beta_peak_freq_r(n_count,iSubj),beta_peak_amp_r(n_count,iSubj)] = find_beta_peak(Power);
                [beta_peak_freq_r(n_count,iSubj),beta_peak_amp_r(n_count,iSubj)] = find_beta_peak_freq(Power,Freqs); 
            end
            %left LFO
            if sum(TF(2,1,:)) ~= 0
%                 Power = squeeze(TF(2,1,Freq_index_1:Freq_index_45));%/(sum(squeeze(TF(2,1,Freq_index_1:Freq_index_45)))+sum(squeeze(TF(2,1,Freq_index_55:Freq_index_95))));
                Power = squeeze(TF(2,1,:));%Freq_index_1:Freq_index_45));%/(sum(squeeze(TF(2,1,Freq_index_1:Freq_index_45)))+sum(squeeze(TF(2,1,Freq_index_55:Freq_index_95))));
%                 [beta_peak_freq_l(n_count,iSubj),beta_peak_amp_l(n_count,iSubj)] = find_beta_peak(Power);
                [beta_peak_freq_l(n_count,iSubj),beta_peak_amp_l(n_count,iSubj)] = find_beta_peak_freq(Power,Freqs);
            end          
        end
    end
end

beta_peak_left = [peak_freqs(2,:)' beta_peak_freq_l'];
beta_peak_right = [peak_freqs(1,:)' beta_peak_freq_r'];

function [peak_freq,peak_amp] = find_beta_peak(Power)
%     Prominence = 8.831512942340412e-04;
%     Prominence = 0.0005;
    Prominence = 1e-13;
    min_beta = 12; %Hz
    max_beta = 35; %Hz
    Power_beta = Power(min_beta:max_beta);
    [beta_amp,beta_freq] = max(Power_beta);
    if beta_freq == 1
        if beta_amp < Power(min_beta-1)
            locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
            if isempty(locs)
                peak_freq = NaN;
                peak_amp = NaN; 
            else
                try 
                    [~,temp] = max(Power_beta(locs));
                    peak_freq = locs(temp) + min_beta-1;
                    peak_amp =  Power(locs(temp) + min_beta-1);
                catch
                    peak_freq = NaN;
                    peak_amp = NaN;                    
                end
            end
        else
            peak_freq = beta_freq+min_beta-1;
            peak_amp = beta_amp;
        end
    elseif beta_freq == size(Power_beta,1)
        if beta_amp < Power(max_beta+1)
            locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
            if isempty(locs)
                peak_freq = NaN;
                peak_amp = NaN; 
            else            
                try 
                    [~,temp] = max(Power_beta(locs));
                    peak_freq = locs(temp) + min_beta-1;
                    peak_amp =  Power(locs(temp) + min_beta-1);
                catch
                    peak_freq = NaN;
                    peak_amp = NaN;                    
                end
            end
        end
    else
        locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
        if isempty(locs)
            peak_freq = NaN;
            peak_amp = NaN; 
        else        
            try 
                [~,temp] = max(Power_beta(locs));
                peak_freq = locs(temp) + min_beta-1;
                peak_amp =  Power(locs(temp) + min_beta-1);
            catch
                peak_freq = NaN;
                peak_amp = NaN;                    
            end
        end
    end
end

function [peak_freq,peak_amp] = find_beta_peak_freq(Power,Freqs)
%     Prominence = 8.831512942340412e-04;
%     Prominence = 0.0005;
    Prominence = 1e-13;

%     min_beta = 12; %Hz
%     max_beta = 35; %Hz
    [val,min_beta]=min(abs(Freqs-12));
    [val,max_beta]=min(abs(Freqs-35));
    Power_beta = Power(min_beta:max_beta);
    [beta_amp,beta_freq] = max(Power_beta);
    if beta_freq == 1% max Power at min_beta
        if beta_amp < Power(min_beta-1)%check for local maxima
            locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
            if isempty(locs)
                peak_freq = NaN;
                peak_amp = NaN; 
            else
                try 
                    [~,temp] = max(Power_beta(locs));
                    peak_freq = Freqs(locs(temp) + min_beta-1);
                    peak_amp =  Power(locs(temp) + min_beta-1);
                catch
                    peak_freq = NaN;
                    peak_amp = NaN;                    
                end
            end
        else
            peak_freq = Freqs(beta_freq+min_beta-1);
            peak_amp = beta_amp;
        end
    elseif beta_freq == size(Power_beta,1)% max Power at max_beta
        if beta_amp < Power(max_beta+1)
            locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
            if isempty(locs)
                peak_freq = NaN;
                peak_amp = NaN; 
            else            
                try 
                    [~,temp] = max(Power_beta(locs));
                    peak_freq = Freqs(locs(temp) + min_beta-1);
                    peak_amp =  Power(locs(temp) + min_beta-1);
                catch
                    peak_freq = NaN;
                    peak_amp = NaN;                    
                end
            end
        else%untere Grenze ist der Beta Peak
            peak_freq = Freqs(beta_freq+min_beta-1);
            peak_amp = beta_amp;
        end
    else% max Power in the beta band
        locs = find(islocalmax(Power_beta,'MinProminence',Prominence));
        if isempty(locs)
            peak_freq = NaN;
            peak_amp = NaN; 
        else        
            try 
                [~,temp] = max(Power_beta(locs));
                peak_freq = Freqs(locs(temp) + min_beta-1);
                peak_amp =  Power(locs(temp) + min_beta-1);
            catch
                peak_freq = NaN;
                peak_amp = NaN;                    
            end
        end
    end
end
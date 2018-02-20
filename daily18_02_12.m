% comparing organoid LFP and baby EEG features
load('/Users/rdgao/Documents/data/Muotri/InfantEEGFeatures/preterm_features.mat')
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/aggregate.mat')
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/names.mat')
%% get recording days
wells = 5:12;
%exclude recordings not done after exchanging media, or pharmacology
exclusion = [5 9 12 15 25 26 34 35 36 40];
recs = setdiff(1:length(dates),exclusion);
dayVec = zeros(1,length(recs));
for day = 1:length(recs)
    %parse numerical date
    date = dates(recs(day)).name(5:end);
    disp(date)
    dayVec(day) = datenum(str2num(date(5:6)), str2num(date(1:2)), str2num(date(3:4)));
end
recDays = dayVec-dayVec(1)+1;
dayVec = recDays/7+7; %use weeks instead

%% refind find peaks
for r_idx = 1:length(dates)
    disp(dates(r_idx))
    for well = wells
        pktimes{r_idx}{well} = kernel_findpeak(nws_smo{r_idx}(:,well), 1, 0.8, 200, 0.001);        
    end
end

%% compute feature table for organoid LFP similar to EEG
wells = 5:12;
num_feat = 23;
psd_freq = 0:0.5:500;
event_thresh=0.1;
fs=1000;
% stack all features
stacked_features = zeros(length(recs), length(wells), num_feat);

% loopdiloop
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    disp(rec)
    for well = wells
        % 7: events per hour
        stacked_features(r_idx, well-4, 7) = length(kerTimes{rec}{well})./T{rec}*3600;
        
        % 8-11: event duration
        % might have to redo peak-detection & windowing
        figure
        title(sprintf('Rec %i, Well %i', rec, well))
        pktimes = kernel_findpeak(nws_smo{rec}(:,well), 1, 0.8, 200, 0.001);
        close
%         for evt = 1:size(kernels{rec}{well},2)            
%             pk_time = round(kerTimes{rec}{well}(evt)*1000);
%             
%             % search in entire nws_vector
%             ker_thresh = nws_smo{rec}(pk_time,well)*event_thresh;
%             % look from [-0.5s, 5s] for first threshold crossing
%             if (pk_time-fs/2)<1
%                 evt_start = find(nws_smo{rec}(1:pk_time,well)<ker_thresh, 1, 'last');
%             else
%                 evt_start = find(nws_smo{rec}((pk_time-fs/2):pk_time,well)<ker_thresh, 1, 'last');
%             end
%             
%             evt_end = find(nws_smo{rec}(pk_time:end,well)<ker_thresh, 1, 'first');
% %             if (pk_time+fs*5)>size(nws_smo{rec},1)
% %                 evt_end = find(nws_smo{rec}(pk_time:end,well)<ker_thresh, 1, 'first');
% %             else
% %                 evt_end = find(nws_smo{rec}(pk_time:(pk_time+fs*5),well)<ker_thresh, 1, 'first');
% %             end
%                         
%             % search in snipped kernel
% %             ker_thresh = kernels{rec}{well}(ker_peak_ind,evt)*event_thresh;
% %             evt_start = find(kernels{rec}{well}(1:ker_peak_ind,evt)<ker_thresh, 1, 'last');
% %             evt_end = find(kernels{rec}{well}(ker_peak_ind:end,evt)<ker_thresh, 1, 'first');
%             
%             
%             evt_dur{rec}{well}(evt) = (evt_end-evt_start)/fs;
%         end
        
        
        % 12-15: RMS, 50%, 5%, and 95% inter-event duration
        ied = diff(kerTimes{rec}{well});
        stacked_features(r_idx, well-4, 12) = mean(ied.^2).^0.5;
        stacked_features(r_idx, well-4, 13) = quantile(ied, 0.5);
        stacked_features(r_idx, well-4, 14) = quantile(ied, 0.05);
        stacked_features(r_idx, well-4, 15) = quantile(ied, 0.95);
        
        
        
        % 20-23: relative spectral power
        % delta (0-3), theta(3-8), alpha(8-15), beta(15-30)
        % sum to 1 from 0-30Hz
        psd = PSDw{rec}{well}(1:61,:);
        psd_norm = psd./repmat(sum(psd,1),61,1);
        stacked_features(r_idx, well-4, 20) = mean(sum(psd_norm(1:7,:),1));
        stacked_features(r_idx, well-4, 21) = mean(sum(psd_norm(8:17,:),1));
        stacked_features(r_idx, well-4, 22) = mean(sum(psd_norm(18:31,:),1));
        stacked_features(r_idx, well-4, 23) = mean(sum(psd_norm(32:end,:),1));
                
    end
end


%%
for feat = 1:num_feat
    
    
end
%% plotting
stacked_features = ctc_EMAfeatures;
feat_num = 23;
plot_feat = double(features_fullepoch(:,feat_num));
plot(age(:,2),plot_feat,'k.', 'markersize', 10)
xlabel('Age (weeks)')
ylabel(feature_cols(feat_num))
params = polyfit(age(:,2),plot_feat,1);
hold on
plot(dayVec, stacked_features(:,:,feat_num), '.b')
plot(dayVec, nanmean(stacked_features(:,:,feat_num),2), '.-')
plot([min(age(:,2)) max(age(:,2))], [min(age(:,2)) max(age(:,2))]*params(1)+params(2), 'r', 'linewidth',2)
hold off



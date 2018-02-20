% scratchpad for developing new organoid electrophys workflow
rec = 30;
well = 5;
for evt = 1:length(pktimes{rec}{well}(:,2))
    plot(collect_spikes(nws_smo{rec}(:,well),[],pktimes{rec}{well}(evt,2),[-1000, 3000]));
end

%%
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/aggregate.mat', 'kerTimes')
% compare number of peaks between the two methods
for rec=1:length(pktimes)
    for well=wells
        disp([rec, well])
        disp(size(pktimes{rec}{well},1)),
        disp(length(kerTimes{rec}{well}))
%         if ~isempty(pktimes{rec}{well})
%             extrapk = setdiff(pktimes{rec}{well}(:,2)/1000, kerTimes{rec}{well});
%             if ~isempty(extrapk)                
%                 disp(length(extrapk))
%             end
%         %    disp('---')
%             extrapk = setdiff(kerTimes{rec}{well},pktimes{rec}{well}(:,2)/1000);
%             if ~isempty(extrapk)
%                 disp(length(extrapk))
%             end
%         end
    end
end

%% plot all organoid and EEG features
%% plotting
stacked_features = ctc_EMAfeatures;
all_feats = [7:15 20:23];
figure
for i=1:length(all_feats)
    subplot(3,5,i)
    feat_num = all_feats(i);
    plot_feat = double(features_fullepoch(:,feat_num));
    %plot_feat = double(features_lowact(:,feat_num));
    plot(age(:,2),plot_feat,'k.', 'markersize', 10)
    xlabel('Age (weeks)')
    ylabel(feature_cols(feat_num))
    params = polyfit(age(:,2),plot_feat,1);
    hold on
    plot(dayVec, stacked_features(:,:,feat_num), '.')
    plot(dayVec, nanmean(stacked_features(:,:,feat_num),2), '.-')
    plot([min(age(:,2)) max(age(:,2))], [min(age(:,2)) max(age(:,2))]*params(1)+params(2), 'r', 'linewidth',2)
    hold off
    xlim(dayVec([1, end]))
end
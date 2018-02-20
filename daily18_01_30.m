% screwing around with premature infant EEG features
cd /Users/rdgao/Documents/Data/Muotri/InfantEEGFeatures
load preterm_features.mat
%% plot each baby differently
feat_num = 13;
figure
hold on
plot_feat = double(features_fullepoch(:,feat_num));
for baby = 1:length(unique(age(:,1)))
    inds = age(:,1)==baby;
    plot(age(inds,2), plot_feat(inds),'-o');
end
xlabel('Age (weeks)')
ylabel(feature_cols(feat_num))

%% aggregate plot
%feat_num = 7;
%plot_feat = 3600./double(features_fullepoch(:,feat_num));
feat_num = 13;
plot_feat = double(features_fullepoch(:,feat_num));
%figure
plot(age(:,2),plot_feat,'k.', 'markersize', 10)
xlabel('Age (weeks)')
ylabel(feature_cols(feat_num))
params = polyfit(age(:,2),plot_feat,1);
hold on
plot([min(age(:,2)) max(age(:,2))], [min(age(:,2)) max(age(:,2))]*params(1)+params(2), 'r', 'linewidth',2)
hold off

%% making the baby & organoid plot
figure
subplot(1,3,1)
% run the above cell
keyboard
title('Premature Infant')

subplot(1,3,2)
% go to corticoid_analysis_flat.m and run the first analysis cell (IEI)
keyboard
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
set(gca,'ytick',0:50:150)
ylabel('Inter-Event Interval (s)')
xlim([8 40])
set(gca,'xtick',10:10:40)
text(10,120, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
title('Organoid')

subplot(1,3,3)
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Inter-Event Interval (s)')
set(gca,'xtick',25:5:40)
hold on
plot(age(:,2),plot_feat,'.', 'color', [0 0 0 0.7], 'markersize', 10)
plot([25 40], [25 40]*params(1)+params(2), 'r', 'linewidth', 2)
hold off
xlim([25 40])
ylim([0 50])
title('Organoid & Premature Infant')

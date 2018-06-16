%%
%cd /Users/rdgao/Documents/Data/Lipton/batch3/
F = dir('*.xlsx');

% Import the data and save as mat
for i=1:length(F)
    disp(F(i).name)
    %[~, ~, raw] = xlsread(F(i).name,'Sheet1');
    [~, ~, raw] = xlsread(F(i).name,'Parfocality'); %batch3 55Hz
    raw = raw(2:end,:);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    stringVectors = string(raw(:,2));
    stringVectors(ismissing(stringVectors)) = '';
    raw = raw(:,[1,3]);
    
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    
    % Create output variable
    data = reshape([raw{:}],size(raw));
    
    % Allocate imported array to column variable names
    ImagePlane = data(:,1);
    ROI = categorical(stringVectors(:,1));
    AverageIntensity = data(:,2);
    
    % Clear temporary variables
    clearvars data raw stringVectors R;
    
    %
    cells = unique(ROI);
    data = zeros(length(find(ROI==cells(1))), length(cells));
    for cell = 1:length(cells)
        inds = find(ROI==cells(cell));
        data(:,cell) = AverageIntensity(inds);
    end
    t = (ImagePlane(inds)-1)*(1/55);    
    save([F(i).name(1:end-5) '.mat'], 't','data')    
    clear data cells inds t ROI ImagePlane AverageIntensity
end

%%
cd /Users/rdgao/Documents/data/Lipton/batch1
F = dir('*.mat');
labels = {'PS1 Het','PS1 WT','TL Het','TL WT'};
close all
% do the FT
for i=1:length(F)
    load(F(i).name)
    data_dtd = detrend(data,'constant');
    fs = 10000/(5*60);
    winlen = fs*6*10;
    noverlap = fs*6*8;
    [P, f_axis] = pwelch(data_dtd,winlen,noverlap,winlen,fs);
    [fitparams{i}, fiterr{i}] = robfitPSD(P, f_axis(7:61), f_axis(2));
        
    figure(1)    
    subplot(2,2,i)
    plot(t, data_dtd)
    title(F(i).name(1:end-4))
    
    figure(2)
    subplot(2,2,i)
    loglog(f_axis,P)
    hold on
    loglog(f_axis,median(P,2), 'k', 'linewidth',2)
    hold off
    xlim([0 15])
    title(F(i).name(1:end-4))
    
    figure(3)
    subplot(1,3,1)
    hold on    
    plot(i, fitparams{i}(:,1), '.k')
    plot(i,mean(fitparams{i}(:,1)), 'ok', 'markerfacecolor', 'k')
    errorbar(i,mean(fitparams{i}(:,1)), std(fitparams{i}(:,1)))
    hold off
    xlim([0 5])
    title('Offset')
    subplot(1,3,2)
    hold on    
    plot(i, fitparams{i}(:,2), '.k')
    plot(i,mean(fitparams{i}(:,2)), 'ok', 'markerfacecolor', 'k')
    errorbar(i,mean(fitparams{i}(:,2)), std(fitparams{i}(:,2)))
    hold off
    xlim([0 5])
    title('Slope')
    subplot(1,3,3)
    hold on    
    plot(i, fiterr{i}, '.k')
    plot(i,mean(fiterr{i}), 'ok', 'markerfacecolor', 'k')
    errorbar(i,mean(fiterr{i}), std(fiterr{i}))
    hold off
    xlim([0 5])
    title('Fit Error')
end
figure(1)
xlabel('Time(s)')
ylabel('Intensity')

figure(2)
xlabel('Frequency (Hz)')
ylabel('Power')

figure(3)
for i=1:3
    subplot(1,3,i)
    xticks(1:4)
    xticklabels(labels)
end
%% batch 2 analysis with drugs
cd /Users/rdgao/Documents/Data/Lipton/batch2
F = dir('*.mat');
close all

% define drug windows
win1 = [1001 4001; 4001 6001; 8001 10001];
win2 = [1 3001; 3001 5001; 8001 10001];
win3 = [2001 5001; 5001 7001; 8001 10001];

fs = 10000/(5*60);
winlen = fs*6*1.2;
noverlap = fs*6;
fit_range = [2, 10];
labels = {'Pre', 'Drug', 'Wash'};
param_labels={'Offset', 'Slope', 'Fit Error'};
fit_inds = 2:14;
diff_plots = {[1,2],[3,2],[1,3]};
PSDs = {};
PFITs = {};
ad_ind = zeros(length(F),1);
for f=1:length(F)    
    disp(F(f).name);
    load(F(f).name);  
    
    if contains(F(f).name,'3000')
        drug_win = win2;
    elseif contains(F(f).name,'5000')
        drug_win = win3;
    else
        drug_win = win1;
    end        
    
    % plot raw data
    figure(1)
    plot(t,data)
    hold on
    plot(t(drug_win(2,1))*[1 1], ylim,'k--')
    plot(t(drug_win(2,2))*[1 1], ylim,'k--')
    hold off
    xlabel('Time (s)')
    ylabel('Intensity')
    title(F(f).name(1:end-4))
    %nice_figure(gcf, ['NS_figures/raw_' F(f).name(1:end-4)], [12, 3.5])
    close
    
    % plot PSD & PSD ratios
    figure(2)
    clear P_rec P_fit
    for seg = 1:3
        % change window length
        %inds = drug_win(seg):(drug_win(seg+1)-1);
        inds = drug_win(seg,1):(drug_win(seg,2)-1);
        data_dtd = detrend(data(inds,:),'constant');
        [P, f_axis] = pwelch(data_dtd,winlen,noverlap,winlen,fs);
        P_rec(:,:,seg) = P;
        f_fit = f_axis(f_axis>=fit_range(1) & f_axis<=fit_range(2));
        [fitparams{f}(seg,:,:), fiterr{f}(seg,:)] = robfitPSD(P, f_fit, f_axis(2));
        P_fit(:,:,seg) = (10.^repmat(fitparams{f}(seg,:,1),length(f_axis),1)).*(f_axis.^fitparams{f}(seg,:,2));
        
        subplot(2,3,seg)
        loglog(f_axis, P, 'color', [0 0 0 0.2])
        hold on
        loglog(f_axis, median(P,2), 'r', 'linewidth',2)
        hold off
        xlim(f_axis([1, end]))
        title(labels(seg))        
    end
    PSDs{f} = P_rec;
    PFITs{f} = P_fit;
    YL = ylim;
    for seg=1:3
        subplot(2,3,seg)
        ylim(YL)
    end
    
    for seg=1:3
       P_ratio = P_rec(:,:,diff_plots{seg}(2))./P_rec(:,:,diff_plots{seg}(1));
       subplot(2,3,seg+3)
       loglog(f_axis, P_ratio,'color', [0 0 0 0.2])
       hold on
       loglog(f_axis, median(P_ratio,2), 'r', 'linewidth',2)
       loglog(f_axis([2 end]), [1 1], 'b--', 'linewidth',1)
       hold off
       xlim(f_axis([1, end]))
       title([labels{diff_plots{seg}(2)} ':' labels{diff_plots{seg}(1)}])       
       ylim([0.01, 100])
    end
    %nice_figure(gcf, ['NS_figures/PSDs_' F(f).name(1:end-4)], [12, 7])
    close
    
    % plot slope-fitting params
    figure(3)
    ind=1;
    for param = {fitparams{f}(:,:,1), fitparams{f}(:,:,2), fiterr{f}}
        subplot(1,3,ind)
        plot(squeeze(param{1}), 'o-', 'color', [0 0 0 0.2])
        hold on
        plot(median(param{1},2), 'ro-', 'linewidth',2)
        hold off
        xlim([0.5, 3.5])
        xticks(1:3)
        xticklabels(labels)        
        title(param_labels(ind))
        ind=ind+1;
    end        
    %nice_figure(gcf, ['NS_figures/fits_' F(f).name(1:end-4)], [12, 3.5])
    close
    
    figure(4)
    % plot PSD change due to pharmacology
    P_ratio = P_rec(:,:,2)./P_rec(:,:,1);
        
    if contains(F(f).name,'WT')
        ad_ind(f)=0;
        C = 'k';
    else
        ad_ind(f)=1;
        C='r';
    end
    
    subplot(1,2,1)
    PSDm = mean(P_rec(:,:,1),2);       
    loglog(f_axis,mean(PSDs{f}(:,:,1)./PFITs{f}(:,:,1),2),C, 'linewidth',1)
    hold on
    xlim(f_axis([1, end]))
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title('1/f-Removed PSD')
    
%     subplot(1,3,2)
%     loglog(f_axis,PSDm./min(PSDm(1:find(f_axis>=10, 1, 'first'))),C, 'linewidth',1)
%     hold on
%     xlim(f_axis([1, end]))
%     title('Pre-Drug Period')
    
    subplot(1,2,2)
    loglog(f_axis, mean(P_ratio,2), C, 'linewidth',1)
    hold on
    xlim(f_axis([1, end]))   
    ylim([0.01, 100])
end
figure(4)
subplot(1,2,2)
loglog(f_axis([2 end]), [1 1], 'b--', 'linewidth',1)
title([labels{2} ':' labels{1}])
ylabel('Power Ratio')
hold off
%% plot fit power
figure
for f = 1:length(PSDs)
    p_pre = PSDs{f}(:,:,1)./PFITs{f}(:,:,1);
    p_post = PSDs{f}(:,:,2)./PFITs{f}(:,:,2);
    loglog(f_axis, mean(p_post./p_pre,2))    
    hold on
end


%% fooofed PSD power
power_range = [3,7];
power_inds = find(f_axis>=power_range(1) & f_axis<=power_range(2));
theta_power = zeros(length(PSDs),1);
for f = 1:length(PSDs)
    theta_power(f) = mean(sum(log10(PSDs{f}(power_inds,:,1)./PFITs{f}(power_inds,:,1)),1));
end
figure
boxplot(theta_power, ad_ind)
hold on
for ind=1:2
plot(ind, theta_power(ad_ind==ind-1), 'ok')
end
xticklabels({'WT', 'AD'})
ylabel('Log10 Power')
title('Mean Log10 Oscillatory Theta (3-7Hz) Power')
[h,p]=ttest2(theta_power(ad_ind==0), theta_power(ad_ind==1), 'tail', 'right')
cohensD = (mean(theta_power(ad_ind==0))-mean(theta_power(ad_ind==1)))/std(theta_power)
nice_figure(gcf, 'comp_theta', [4,4])
%% WT vs. AD low-freq power
power_range = [0.1, 1];
power_inds = find(f_axis>=power_range(1) & f_axis<=power_range(2));
low_freq_power = zeros(length(PSDs),1);
for f = 1:length(PSDs)
%     P_adj = PSDs{f}(:,:,1)./repmat(PSDs{f}(find(f_axis>=10,1),:,1),length(f_axis),1);
%     low_freq_power(f) = mean(sum(log10(P_adj(power_inds,:)),2),1);
    low_freq_power(f) = mean(sum(log10(PSDs{f}(power_inds,:,1)./PFITs{f}(power_inds,:,1)),1));
end
figure
boxplot(low_freq_power, ad_ind)
hold on
for ind=1:2
    plot(ind, low_freq_power(ad_ind==ind-1), 'ok')
end
ylabel('Log10 Power')
xticklabels({'WT', 'AD'})
title({'Mean Low Frequency (0.2-1Hz) Power'})
nice_figure(gcf, 'comp_lowfreqWTAD', [4,4])
[h,p]=ttest2(low_freq_power(ad_ind==0), low_freq_power(ad_ind==1))


%% drug treatment low_freq power change
power_range = [0.2, 1];
power_inds = find(f_axis>=power_range(1) & f_axis<=power_range(2));
low_freq_pr = zeros(length(PSDs),1);
for f = 1:length(PSDs)
    low_freq_pr(f) = mean(sum(log10(PSDs{f}(power_inds,:,2)./PSDs{f}(power_inds,:,1)),1));
end
figure
boxplot(low_freq_pr, ad_ind)
hold on
for ind=1:2
    plot(ind, low_freq_pr(ad_ind==ind-1), 'ok')
end
ylabel('Power Ratio')
xticklabels({'WT', 'AD'})
title({'Mean Low Frequency (0.2-1Hz) Power','Change after NitroSynapsin'})
nice_figure(gcf, 'comp_lowfreq', [4,4])
[h,p]=ttest2(low_freq_pr(ad==0), low_freq_pr(ad==1))

%% batch 3 (May 26)
cd('/Users/rdgao/Documents/data/Lipton/batch3/Calcium imaging_50 Hz')
% 50Hz exps
F = dir('*.mat');
labels = {'TL Het','TL WT'};
close all
% do the FT
for i=1:length(F)
    load(F(i).name)
    if contains(F(i).name,'WT')
        ad_ind(i)=0;
        C = 'k';
    else
        ad_ind(i)=1;
        C='r';
    end
    data_dtd = detrend(data,'constant');
    fs = 50;
    winlen = fs*6*3;
    noverlap = fs*4*3;
    thresh_freq = 0.8;
    [P, f_axis] = pwelch(data_dtd,winlen,noverlap,winlen,fs);
    [fitparams{i}, fiterr{i}] = robfitPSD(P, f_axis(2:find(f_axis>=thresh_freq,1)), f_axis(2));
        
    figure(1)    
    subplot(2,7,i)
    plot(t, data_dtd)
    title(F(i).name(1:end-4))
    xlim([0, t(end)])
    
    figure(2)
    subplot(2,7,i)
    loglog(f_axis,P)
    hold on
    loglog(f_axis,median(P,2), 'k', 'linewidth',2)
    hold off
    xlim([0 15])
    title(F(i).name(1:end-4))
    
    figure(3)
    subplot(1,3,1)
    hold on    
    plot(i, fitparams{i}(:,1), '.', 'color', C)
    plot(i,mean(fitparams{i}(:,1)), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fitparams{i}(:,1)), std(fitparams{i}(:,1)), 'color', C)
    hold off
    %xlim([0 5])
    title('Offset')
    subplot(1,3,2)
    hold on    
    plot(i, fitparams{i}(:,2), '.', 'color', C)
    plot(i,mean(fitparams{i}(:,2)), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fitparams{i}(:,2)), std(fitparams{i}(:,2)), 'color',C)
    hold off
    %xlim([0 5])
    title('Slope')
    subplot(1,3,3)
    hold on    
    plot(i, fiterr{i}, '.', 'color', C)
    plot(i,mean(fiterr{i}), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fiterr{i}), std(fiterr{i}), 'color', C)
    hold off
    %xlim([0 5])
    title('Fit Error')
end
figure(1)
xlabel('Time(s)')
ylabel('Intensity')

figure(2)
xlabel('Frequency (Hz)')
ylabel('Power')

figure(3)
%% 55Hz exps
cd('/Users/rdgao/Documents/data/Lipton/batch3/Calcium imaging_55 Hz')
F = dir('*.mat');
labels = {'TL Het','TL WT'};
close all
% do the FT
for i=1:length(F)
    load(F(i).name)
    if contains(F(i).name,'WT')
        ad_ind(i)=0;
        C = 'k';
    else
        ad_ind(i)=1;
        C='r';
    end
    data_dtd = detrend(data,'constant');
    fs = 55;
%     winlen = fs*6*2;
%     noverlap = fs*4*2;
    [P, f_axis] = pwelch(data_dtd,winlen,noverlap,winlen,fs);
    [fitparams{i}, fiterr{i}] = robfitPSD(P, f_axis(2:find(f_axis>=thresh_freq,1)), f_axis(2));
        
    figure(1)    
    subplot(2,6,i)
    plot(t, data_dtd)
    xlim([0,t(end)])
    title(F(i).name(1:end-4))
    
    figure(2)
    subplot(2,6,i)
    loglog(f_axis,P)
    hold on
    loglog(f_axis,median(P,2), 'k', 'linewidth',2)
    hold off
    xlim([0 15])
    title(F(i).name(1:end-4))
    
    figure(3)
    subplot(1,3,1)
    hold on    
    plot(i, fitparams{i}(:,1), '.', 'color', C)
    plot(i,mean(fitparams{i}(:,1)), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fitparams{i}(:,1)), std(fitparams{i}(:,1)), 'color', C)
    hold off
    %xlim([0 5])
    title('Offset')
    subplot(1,3,2)
    hold on    
    plot(i, fitparams{i}(:,2), '.', 'color', C)
    plot(i,mean(fitparams{i}(:,2)), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fitparams{i}(:,2)), std(fitparams{i}(:,2)), 'color',C)
    hold off
    %xlim([0 5])
    title('Slope')
    subplot(1,3,3)
    hold on    
    plot(i, fiterr{i}, '.', 'color', C)
    plot(i,mean(fiterr{i}), 'ok', 'markerfacecolor', C)
    errorbar(i,mean(fiterr{i}), std(fiterr{i}), 'color', C)
    hold off
    %xlim([0 5])
    title('Fit Error')
end
figure(1)
xlabel('Time(s)')
ylabel('Intensity')

figure(2)
xlabel('Frequency (Hz)')
ylabel('Power')

figure(3)
% for i=1:3
%     subplot(1,3,i)
%     xticks(1:4)
%     xticklabels(labels)
% end
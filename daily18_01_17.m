%%
cd /Users/rgao/Documents/Data/Lipton
F = dir('*.xlsx');

% Import the data and save as mat
for i=1:length(F)
    disp(F(i).name)
    [~, ~, raw] = xlsread(F(i).name,'Sheet1');
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
    t = (ImagePlane(inds)-1)*0.03;
    save([F(i).name(1:end-5) '.mat'], 't','data')    
    clear data cells inds t ROI ImagePlane AverageIntensity
end

%%
cd /Users/rgao/Documents/Data/Lipton
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
%%



%% 
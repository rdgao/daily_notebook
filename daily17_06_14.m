% convert infant .dat eeg files to .mat files
cd /Users/rgao/Documents/Data/Muotri/InfantEEG
F = dir('*.dat');
%%
for ff=1:length(F)
    disp(F(ff).name)
    save_file = [F(ff).name(1:end-3) 'mat'];
    if ~exist(save_file)
        [A Delim Nhead] = importdata(F(ff).name, ' ', 14);
        labels = A.textdata{9};
        rate = A.textdata{5};
        data = A.data;
        save([F(ff).name(1:end-3) 'mat'], 'data', 'rate', 'labels');
    else
        disp('skipped')
    end
end

%%
F = dir('*.mat');
file_num = 6;
load(F(file_num).name)
% filter data
data_filt = eegfilt(data',fs,0,45)';
plot((1:length(data_filt))/fs,data_filt(:,1:5))
%%
clear SPG
fs = 167;
for chan = 1:14
    [SPG(:,:,chan), Faxis, T] = spectrogram(data(:,chan),fs*2,45,fs*2,fs);
end
%% plot spectrogram
figure
imagesc(T, Faxis(1:100), log10(abs(SPG(1:100,20:end))))
%% plot PSD 
figure
loglog(Faxis,squeeze(median(abs(SPG),2)))


%% plotting and labeling eeg12466_17089.mat
labels = importdata('eeg12466_17089.tcl');
evt_t = labels.data(2:32,1);
evt_names = labels.textdata(37:end,2);
load('eeg12466_17089.mat')
fs = 167;
data = detrend(data,'constant');
data_filt = eegfilt(data',fs,0,45)';
%%
figure
plot((1:length(data_filt))/fs,data_filt(:,5))
%plot((1:length(data_filt))/fs,data(:,1))
hold on
plot([evt_t evt_t]/fs, [-200 300], 'k--')
hold off
for i=1:length(evt_t)
    text(evt_t(i)/fs, -200, evt_names(i))
end
xlim([1790 1850])
xlabel('Time (s)')

%%
fs = 167;
burst_per = 1790*fs:1850*fs;
burst_inds = 1790:1850;
fs = 167;
SPG = [];
for chan = 1:14
    %[SPG(:,:,chan), Faxis, T] = spectrogram(data(burst_per,chan),fs*2,fs,fs*2,fs);
    [SPG(:,:,chan), Faxis, T] = spectrogram(data(:,chan),fs*2,fs,fs*2,fs);
end
figure
imagesc(T(burst_inds), Faxis(1:21), (abs(SPG(1:21,burst_inds,1))))
xlabel('Time (s)')
ylabel('Freq (Hz)')

%%
%figure
loglog(Faxis,mean(abs(SPG(:,burst_inds,1)),2))
    
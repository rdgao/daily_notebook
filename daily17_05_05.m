%% re-doing spike detection
cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017
load well_5.mat
%%
fs = 12500;
chan = 10;
f_o = [20, 50, 100, 200];
for i=1:4
    b = fir1(f_o(i),[300 3000]/(fs/2));
    filtered(:,i) = filtfilt(b,1,MEA(:,chan));    
end
figure 
plot(t, filtered, 'linewidth',1)

%% and causal vs. non-causal
figure
chan = 2;
f_o = 50;
b = fir1(f_o,[300 3000]/(fs/2));
f = filter(b,1,MEA(:,chan));
ff = filtfilt(b,1,MEA(:,chan));
plot(ff, 'linewidth', 1)
hold on
plot(f(f_o/2+1:end), 'linewidth', 1)
hold off

%% see the spikes
[spks, spkcnt] = spike_detect_abs(filtered, fs, 5.5);
%%
for i=1:4
    figure
    plot(filtered(:,i), 'linewidth',1)
    hold on
    plot(spks{i},filtered(spks{i}), 'or')
    hold off
    title(f_o(i))
end

%% compare with butterworth filter
chan=10;
f = butterpass(MEA(:,chan),fs,[300 3000],3);
figure
plot(f, 'linewidth',1)
hold on
plot(filtered(:,4),'linewidth',1)
hold off
figure
for i=1:4
    subplot(2,2,i)
    plot(f,filtered(:,i), '.')
    line([-1 1]*2e-4, [-1 1]*2e-4, 'color', 'r')
    xlabel('butter')
    ylabel('fir')
    title(f_o(i))
end
%%
[spksbw, spkcntbw] = spike_detect_abs(f, fs, 5.5);
figure
plot(f, 'linewidth',1)
hold on
plot(spksbw{1}, f(spksbw{1}), 'o')
hold off

%% look at spatial neighbors
fbw = butterpass(MEA,fs,[300 3000],3);
%%
chans = [10, 9, 11, 2, 18];
figure
hold on
for i=chans
    plot(fbw(:,i), 'linewidth', 1)
end
hold off

%% compare axion spikes vs mine
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoid spikes/'
load CTC031017_Axion.mat
%%
[spksbw, spkcntbw] = spike_detect_abs(fbw(:,10), fs, 4);
figure
plot(fbw(:,10), 'linewidth',1)
hold on
plot(spikes{5,10},fbw(spikes{5,10},10), 'o')
plot(spksbw{1},fbw(spksbw{1},10), 'o')
hold off
legend({'data', 'axion', 'rg: bw'})
%% look at ISI
[spksbw, spkcntbw] = spike_detect_abs(fbw, fs, 4);
%%
chan = 1;
figure
hold on
for chan = 1:64
    plot(spksbw{chan}(1:end-1)/fs, diff(spksbw{chan})/fs, '.')
end
%% population spike
spks = vertcat(spksbw{:})/fs;
bins = 0:0.0004:2;
figure
loglog(bins, (hist(diff(sort(spks)), bins)), 'o-')
xlabel('population ISI')
ylabel('count')
%%
figure
dt = 0.0004*10;
bins = dt:dt:2;
for chan=1:64
    loglog(bins,hist(diff(spksbw{chan})/fs,bins), 'linewidth',1, 'color', [0 0 0 0.1])
    hold on
end
hold off
xlabel('channel ISI')
ylabel('count')

%% 
bsp = squeeze(binarize_spikes(t,fs,spksbw,1000));
nws = sum(bsp);
win_len = 50;
smo_win = gausswin(win_len)/sum(gausswin(win_len));
nws_smo = conv(nws, smo_win, 'same');
plot(nws_smo)

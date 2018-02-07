% load data
load /Users/rgao/Documents/Data/Carmena/020608B/LFP65_128.mat
%%
% bandpass oscillations
B = fir1(625,[8 13]/(fs/2));
alpha = filtfilt(B,1,LFP);

B = fir1(333,[15 35]/(fs/2));
beta = filtfilt(B,1,LFP);

B = fir1(143,[35 70]/(fs/2));
gamma = filtfilt(B,1,LFP);

%%
chan = 1;
plot((1:length(LFP))/1000,LFP(:,chan), 'k')
hold on
plot((1:length(LFP))/1000,alpha(:,chan))
plot((1:length(LFP))/1000,beta(:,chan))
plot((1:length(LFP))/1000,gamma(:,chan))
hold off

%%
smo_len = 1000;
figure
beta_amp = medfilt1(abs(hilbert(beta(1:end-1,chan))),smo_len);
beta_freq = medfilt1(diff(unwrap(angle(hilbert(beta(:,chan)))))*fs/(2*pi),smo_len);
gamma_amp = medfilt1(abs(hilbert(gamma(1:end-1,chan))),smo_len);
gamma_freq = medfilt1(diff(unwrap(angle(hilbert(gamma(:,chan)))))*fs/(2*pi), smo_len);
plot(beta_freq,gamma_freq, '.')
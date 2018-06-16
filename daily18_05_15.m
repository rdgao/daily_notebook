% loading data from binary
file='~/Documents/data/CRCNS/pfc2/EE.049/EE.049.eeg';
[data, OrigIndex]= LoadBinary(file, 1:104);
lfp=data;
%% 
% load lfp, spikes & behavior
load('EE.049_Behavior.mat');
load('EEG049_LFP.mat')
%%
% compute PSD
[Pxx,F] = pwelch(lfp',1250,1250/2,1250,1250);
%%
% bin spikes
%spikes = {};
for cell = 1:length(spikeph(:,1))
    spikes{cell} = spiket(spikeind==cell);    
end
% 
bsp = squeeze(binarize_spikes(1790,20000,spikes,1000));
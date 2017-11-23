%% looking at rett organoids
cd('/Users/rgao/Documents/Data/Muotri/Rett_CTC/rtt_021017')
load LFP_Sp.mat

%%
fs = 12500;
nws = spikes2pop(spikes, t_ds, fs, fs_ds);
nws_smo = tsmovavg(nws', 'e', 100)';
%% GO TO 17_05_15 notebook

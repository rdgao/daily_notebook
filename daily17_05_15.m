% analyzing in detail WTF went wrong with the control data
cd('/Volumes/My Passport for Mac/Dish/CDKL5_CTC/CDKL5_CTC_022417/')
load LFP_Sp.mat
%%
fs=12500;
bsp = binarize_spikes(ceil(t_ds(end)), fs,spikes,fs_ds);
nws = squeeze(sum(bsp,2))'; 
nws = nws(1:length(t_ds),:);
nws_smo = tsmovavg(nws', 'e', 25);
plot_tight(nws_smo', [3 4], [], [0 3])
%% raster and population spiking
well = 9;
figure
spike_raster(t_s,spikes,well,1)
hold on
plot(t_ds, tsmovavg(nws(:,well)','e', 25)*10, 'r');
hold off

%% well oscillations 
well = 9;
data = LFP{well};
delta = butterpass(data, fs_ds, [0.2 4], 3);
MUA = butterpass(data, fs_ds, [100 300], 3);
data = butterpass(data, fs_ds, [0.2 55], 3);
%%
cmap = [colormap('parula'); flipud(colormap('parula'))]; %make circular
figure
imagesc(t_ds, [], angle(hilbert(delta(:,:)))')
colormap(cmap)
colorbar
hold on
plot(t_ds, tsmovavg(nws(:,well)','e', 25)*10, 'color',[0 0 0 0.5]);
hold off
set(gca, 'ydir', 'normal')

%% delta & mua for matt
well = 11;
chan = 31;
figure
plot(t_ds, data(:,chan), 'linewidth', 1)
hold on
plot(t_ds, abs(hilbert(MUA(:,chan))), 'g', 'linewidth', 1)
plot(t_ds, nws_smo(well,:)/5e5, 'k', 'linewidth', 1)
plot(t_ds, delta(:,chan),'linewidth', 2)
hold off

%% spikes and oscillation
well = 11;
chan = 31;
data = LFP{well}(:,chan);
delta = butterpass(data, fs_ds, [0.5 4], 3);
figure
spike_raster(t_s,spikes,well,1)
hold on
plot(t_ds, delta*2e6+30,'linewidth', 1)
hold off
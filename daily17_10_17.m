% read in klusta sorted spikes in matlab and compare with my spikes
% raw data file: CTC_031017
fs=12500;

%% klusta prb: separate shank
filename = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/spike_sorting/klusta_separate_shank/binary_well_5.kwik';
cells = {};
for i=1:64
    spks{i} = hdf5read(filename, sprintf('/channel_groups/%i/spikes/time_samples',i-1));
    clus{i} = hdf5read(filename, sprintf('/channel_groups/%i/spikes/clusters/main',i-1));
    for nc=1:length(unique(clus{i}))
        cells = [cells double(spks{i}(clus{i}==nc))];
    end
end
% binarize spikes
bsp_klu = squeeze(binarize_spikes(248, fs, cells, 1000));

%% klusta prb: unconnected electrode, single shank
filename = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/spike_sorting/klusta_unconnected_electrodes/binary_well_5.kwik';
spks = hdf5read('binary_well_5.kwik', '/channel_groups/0/spikes/time_samples');
clus = hdf5read('binary_well_5.kwik', '/channel_groups/0/spikes/clusters/main');
cells = {};
for nc=1:length(unique(clus))
   cells{nc} = double(spks(clus==nc));
end
% binarize spikes
bsp_klu = squeeze(binarize_spikes(248, fs, cells, 1000));

%% get my spikes
myfile = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017/LFP_Sp.mat';
load(myfile, 'spikes');
bsp_rg = squeeze(binarize_spikes(248, fs, spikes(5,:), 1000));

%% compare
bsp_klu = bsp_klu(3:end,:);
t = 0:0.001:248-0.001;
figure
plot(t, sum(bsp_klu))
hold on
plot(t, sum(bsp_rg))
hold off
legend('Klusta', 'RG')

figure
plot(t, tsmovavg(sum(bsp_klu),'w',gausswin(100)'))
hold on
plot(t, tsmovavg(sum(bsp_rg),'w',gausswin(100)'))
hold off
legend('Klusta', 'RG')
%% correlation of pop. spiking vector
figure
plot(tsmovavg(sum(bsp_klu),'w',gausswin(100)'), tsmovavg(sum(bsp_rg),'w',gausswin(100)'), '.');
hold on
plot(sum(bsp_klu), sum(bsp_rg),'.')
hold off
xlabel('klusta')
ylabel('RG')
title(corr(sum(bsp_klu)', sum(bsp_rg)'))

%% cross correlation of cells
cross_prod = bsp_rg*bsp_klu';
[val, ind] = max(cross_prod);
[val ind] = sort(ind);
cross_prod = cross_prod(:, ind);
figure
imagesc(cross_prod./repmat(max(cross_prod),[size(bsp_rg,1),1]))
xlabel('klusta cells')
ylabel('RG cells')
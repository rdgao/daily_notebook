% with reference to daily19_04_17
% comparing axion spikes with mine
file = AxisFile('34-018_CTC_031017_re-recording(000).spk');
S = file.DataSets.LoadData('B1');
spks = [S{2,1,1,1}.Start]*12500;

spikes = load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017/LFP_Sp.mat', 'spikes');
%%
bspAx = bin_spikes(spks, [12500, 242], 1000);
bspAxAdj = bin_spikes(spks+11, [12500, 242], 1000);
bsp = bin_spikes(spikes.spikes{5,1}, [12500, 242], 1000);
%axion spike times have to be adjusted by the "Start" time from the Axion file

figure
hold on
plot(bsp)
plot(bspAxAdj*0.9)
hold off

%% re-interpret axion spikes and compare with mine
AX = load('CTC031017_Axion.mat');
RG = load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017/LFP_Sp.mat', 'spikes');

%%
spksAX = [AX.spikes{5,:}];
spksRG = vertcat(RG.spikes{5,:})';
%%
bspAx = bin_spikes(spksAX, [12500, 250], 1000);
bspRG = bin_spikes(spksRG, [12500, 250], 1000);
figure
plot(bspRG)
hold on
plot(bspAx)
hold off
legend({'RG','AX'})
%%
fs=12500;
bspAx = binarize_spikes(0:fs*250,12500,AX.spikes,1000);
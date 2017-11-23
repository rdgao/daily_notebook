% dish PAC figures
cd('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416')
load LFP_Sp.mat
%%
fs = 1000;
well = 12;
chan = 45; %30/45
delta = eegfilt(LFP{well}(:,chan)',fs,0,4)';
delta_ph = angle(hilbert(delta));
gamma = eegfilt(LFP{well}(:,chan)',fs,200,300)';
gamma_amp = abs(hilbert(gamma));
gamma_amp = eegfilt(abs(hilbert(gamma))',fs,0,4)';
%%
figure
plot(t_ds, delta, 'color', [0 0 0 0.5], 'linewidth',1)
hold on
plot(t_ds, gamma_amp, 'color', [1 0 0 0.5], 'linewidth',1)
hold off
figure
plot(delta,gamma_amp, '.')
plot(angle(hilbert(delta)),gamma_amp, '.')

%%
load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/aggregate.mat', 'kerTimes')
load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/aggregate.mat', 'LFP_ker')
%%
agg_date = 18;
pktimes = kerTimes{agg_date}{well}*fs;
dl_epc = collect_spikes(delta, [], pktimes, [-500 2500]);
dlph_epc = collect_spikes(angle(hilbert(delta)), [], pktimes, [-500 2500]);
gm_epc = collect_spikes(gamma, [], pktimes, [-500 2500]);
gmap_epc = collect_spikes(gamma_amp, [], pktimes, [-500 2500]);

figure
colors = get(gca,'colororder');
plot((-500:2500)/fs, LFP_ker{agg_date}{well}(:,1,chan), 'color', [colors(1,:) 0.5], 'linewidth',1);
hold onPAC
plot((-500:2500)/fs,dl_epc(:,1), 'color', 'k', 'linewidth',1)
plot((-500:2500)/fs,1e6*gmap_epc(:,1).^2, 'color', colors(2,:), 'linewidth',1)
hold off
legend({'Raw LFP', 'Delta', 'Gamma Amp.'}, 'location', 'southeast')
xlabel('Time (s)')
ylabel('Voltage / Power')
xlim([-0.5 2.5])
set(gca, 'xtick', [0 2])
set(gca, 'ytick', 0)

%% bin gamma power
bin_edge = linspace(-pi,pi,20);
[count, bins] = histc(angle(hilbert(delta)), bin_edge);
gm_binned = {};
for b = 1:length(bin_edge)-1
    gm_binned{b} = gamma_amp(bins==b);
end
N = cellfun(@length,gm_binned);
M = cellfun(@mean,gm_binned);
S = cellfun(@std,gm_binned)/sqrt(N);
figure
plot_filled(bin_edge(1:end-1), M,S,'k')
xlabel('Delta Phase')
ylabel('Gamma Power')

%%
[count, bins] = histc(reshape(dlph_epc,1,[]), bin_edge);
gm_lin = reshape(gmap_epc,1,[]);
gm_binned = {};
for b = 1:length(bin_edge)-1
    gm_binned{b} = gm_lin(bins==b);
end
N = cellfun(@length,gm_binned);
M = cellfun(@mean,gm_binned);
S = cellfun(@std,gm_binned)/sqrt(N);

figure
plot_filled(bin_edge(1:end-1), M,S,'k')
xlabel('Delta Phase')
ylabel('Gamma Power')

%% compute PAC during and outside of events
load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/aggregate.mat', 'kerTimes')

agg_date = 18;
fs = 1000;
well = 12;
evtimes = kerTimes{agg_date}{well}*fs;
win = -500:2500;
num_bins = 21;
bin_edge = linspace(-pi,pi,num_bins);

% get event and non-event indices
ev_inds = reshape((repmat(evtimes,1,length(win))+repmat(win,length(evtimes),1))',1,[]);
nev_inds = setdiff(1:length(LFP{well}(:,1)),ev_inds);

% filtering
chan = 5;
for chan=1:64
delta = eegfilt(LFP{well}(:,chan)',fs,0,4)';
delta_ph = angle(hilbert(delta));
gamma = eegfilt(LFP{well}(:,chan)',fs,200,400)';
gamma_amp = abs(hilbert(gamma)).^2;

[count, ev_bins] = histc(delta_ph(ev_inds), bin_edge); % event delta phase
[count, nev_bins] = histc(delta_ph(nev_inds), bin_edge); % non-event delta phase
gamma_ev = gamma_amp(ev_inds);
gamma_nev = gamma_amp(nev_inds);

gm_ev = {};
gm_nev = {};
for b = 1:length(bin_edge)-1    
    gm_ev{b} = gamma_ev(ev_bins==b);
    gm_nev{b} = gamma_nev(nev_bins==b);   
end


%figure
subplot(8,8,chan)
% plot(bin_edge(1:end-1),cellfun(@mean,gm_ev)/sum(cellfun(@mean,gm_ev)))
% hold on
% plot(bin_edge(1:end-1),cellfun(@mean,gm_nev)/sum(cellfun(@mean,gm_nev)))
% hold off

plot(bin_edge(1:end-1),zscore(cellfun(@mean,gm_ev)))
hold on
plot(bin_edge(1:end-1),zscore(cellfun(@mean,gm_nev)))
hold off

% stairs(linspace(-pi,pi,num_bins-1),cellfun(@mean,gm_ev)/sum(cellfun(@mean,gm_ev)), 'linewidth', 2)
% hold on
% stairs(linspace(-pi,pi,num_bins-1),cellfun(@mean,gm_nev)/sum(cellfun(@mean,gm_nev)), 'linewidth', 2)
% hold off

xlim([-pi pi])
set(gca, 'xtick', [-pi,0, pi])
set(gca, 'xticklabel', {'-pi','0', 'pi'})
end

xlabel('Delta Phase')
ylabel('Normalized Gamma Amp.')
legend({'Event', 'Non-Event'})
title(sprintf('Well %i, Chan %i', well, chan))
%% PAC analysis
well = 12;
agg_date = 18;
fs = 1000;
evtimes = kerTimes{agg_date}{well}*fs;
win = -500:2500;
num_bins = 21;
bin_edge = linspace(-pi,pi,num_bins);
GM_EV = {};
GM_NEV = {};

% get event and non-event indices
ev_inds = reshape((repmat(evtimes,1,length(win))+repmat(win,length(evtimes),1))',1,[]);
nev_inds = setdiff(1:length(LFP{well}(:,1)),ev_inds);

% filtering
for chan=1:64
    delta = eegfilt(LFP{well}(:,chan)',fs,0,4)';
    delta_ph = angle(hilbert(delta));
    gamma = eegfilt(LFP{well}(:,chan)',fs,200,400)';
    gamma_amp = abs(hilbert(gamma)).^2;
    
    [count, ev_bins] = histc(delta_ph(ev_inds), bin_edge); % event delta phase
    [count, nev_bins] = histc(delta_ph(nev_inds), bin_edge); % non-event delta phase
    gamma_ev = gamma_amp(ev_inds);
    gamma_nev = gamma_amp(nev_inds);
    
    gm_ev = {};
    gm_nev = {};
    for b = 1:length(bin_edge)-1
        gm_ev{b} = gamma_ev(ev_bins==b);
        gm_nev{b} = gamma_nev(nev_bins==b);
    end
    GM_EV{chan} = gm_ev;
    GM_NEV{chan} = gm_nev;
end
%%
max_ph = zeros(64,2);
MI = zeros(64,2);
N = num_bins-1;
figure
hold on
for chan = 1:64
    [y, ind] = max(cellfun(@mean,GM_EV{chan}));
    %gmap_ev(:,chan) =  circshift(zscore(cellfun(@mean,GM_EV{chan})),-ind+1,2); %zscore
    gmap_ev(:,chan) =  circshift(cellfun(@mean,GM_EV{chan})/sum(cellfun(@mean,GM_EV{chan})),-ind+1,2); %normed (pdf for MI)
    MI(chan,1) = log2(N)+sum(gmap_ev(:,chan).*log2(gmap_ev(:,chan))); %MI calculation
    max_ph(chan,1) = ind;
    plot(gmap_ev(:,chan), 'color', [0 0 0 0.2]);
    
    [y, ind] = max(cellfun(@mean,GM_NEV{chan}));
    %gmap_nev(:,chan) =  circshift(zscore(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
    gmap_nev(:,chan) =  circshift(cellfun(@mean,GM_NEV{chan})/sum(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
    MI(chan,2) = log2(N)+sum(gmap_nev(:,chan).*log2(gmap_nev(:,chan)));
    max_ph(chan,2) = ind;
    plot(gmap_nev(:,chan), 'color', [1 0 0 0.2]);    
       
end
hold off
%%
figure
plot_filled(linspace(-pi,pi,num_bins-1), mean(gmap_nev,2), std(gmap_nev,1,2), colors(2,:))
hold on
plot_filled(linspace(-pi,pi,num_bins-1), mean(gmap_ev,2), std(gmap_ev,1,2), 'k')
hold off
set(gca, 'xtick', [-pi 0 pi])
set(gca, 'xticklabel', {'-\pi', '0', '\pi'})
xlabel('Delta Phase')
ylabel('Normalized Gamma Amp.')
axes_h = gca;
legend(axes_h.Children([1 3]), {'Event' 'Non-Event'}, 'location', 'best')
nice_figure(gcf, [fig_folder '4_PAC_MI'],[3 3])
%%
figure
boxplot(MI, 'colors', [0 0 0; colors(2,:)], 'outliersize', 1)
%plot(1,MI(:,1), 'k.')
ylabel('Modulation Index')
set(gca,'xticklabel',{'Event' 'Non-Event'})
nice_figure(gcf, [fig_folder '4_MI'],[2 3])
%%
figure
h = rose(bin_edge(max_ph(:,1)));
h.Color = 'k'
hold on
h = rose(bin_edge(max_ph(:,2)));
hold off
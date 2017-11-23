% dev for figures for rett 
cd '/Users/rgao/Documents/Data/Muotri/Rett_CTC/rtt_021017'
fig_folder =  '~/Dropbox/Research/Reports/Muotri/CorticoidFigs/';
load LFP_Sp.mat
%%
fs = 12500;
fs_ds=1000;
smo_len = 100;
smo_mask = gausswin(smo_len)./sum(gausswin(smo_len));
bsp = binarize_spikes(ceil(t_ds(end)), fs,spikes,fs_ds);
nws = squeeze(sum(bsp,2))'; nws = nws(1:length(t_ds),:);
for well = 1:12
    nws_smo(:,well) = conv(nws(:,well), smo_mask,'same');
end
plot_tight(nws_smo, [3 4], [], [0 3])
%% example figure
figure
colors = get(gca,'colororder');
subplot(1,2,1)
plot(t_ds, nws_smo(:,10), 'k', 'linewidth', 1)
xlim([35 40])
set(gca, 'xtick', [35 40])
set(gca, 'xtickLabel', {'0', '5'})
xlabel('Time (s)')
ylim([0 1])
set(gca, 'ytick', [0 1])
ylabel('Population Spiking (spk/s)')
box off
legend('Control')

subplot(1,2,2)
plot(t_ds, nws_smo(:,3), 'color', colors(2,:), 'linewidth', 1)
xlim([35 40])
set(gca, 'xtick', [35 40])
set(gca, 'xtickLabel', {'0', '5'})
ylim([0 1])
set(gca, 'ytick', [0 1])
box off
legend('RTT')
%nice_figure(gcf, [fig_folder '5_rettosc'],[6 3])

%% collecting kernels and quantification
ker_win = [-500,2500];
smo_len = 100;
min_event_spike = 0.2; % peak is at least 1spk tall to be considered event

wells = 1:12;
for well = 1:12
    [PK,IND] =findpeaks(nws_smo(:,well), 'minpeakheight', max(max(nws_smo(:,well))*0.75,min_event_spike), 'minpeakdistance',fs_ds);
    kernels{well} = collect_spikes(nws_smo(:,well),[],IND,ker_win);
    %record spike stamps
    kerTimes{well} = IND./fs_ds;
end

%% plotting all kernels
for well = 1:12
    subplot(3,4,well)
    if ~isempty(kernels{well})
        plot((ker_win(1):ker_win(end))/1000, kernels{well}, 'k', 'linewidth',1)
    end    
    set(gca, 'xtick', '')
    xlim([-0.5 2.5])
end
subplot(3,4,9)
set(gca, 'xtick', [0 2])
xlabel('Time (s)')
ylabel('Population Spikes')
nice_figure(gcf, [fig_folder '5supp_rettkers'],[8 6])

%% bar plot of number of events
nker = cellfun(@length, kerTimes);
inds = {[1 2 5 6 9 10] [3 4 7 8 11 12]};
figure
hold on
bar(1, mean(nker(inds{1})), 'w', 'edgecolor', 'k')
line([1 1], mean(nker(inds{1}))+[-std(nker(inds{1})) std(nker(inds{1}))], 'color', 'k')
bar(2, mean(nker(inds{2})))
scatter(1+randn([1,6])/10, nker(inds{1}), 15, 'k', 'fill')
scatter(2+randn([1,6])/10, nker(inds{2}), 15, colors(2,:), 'fill')
hold off
xlim([0.5 2.5])
ylabel('Network Events')
set(gca, 'ytick',[0 7])
set(gca, 'xtickLabel',{'Control' , 'RTT'})
box off
[h p] = ttest2(nker(inds{1}),nker(inds{2}))
[p h] = ranksum(nker(inds{1}),nker(inds{2}))
%nice_figure(gcf, [fig_folder '5_rettnker'],[2 3])

%% population firing rate rate
thresh = 3*6; % 1 spike for 10 second threshold
pop_fr = (sum(spike_cnt.*(spike_cnt>thresh),2)./sum(spike_cnt>thresh,2))/t_ds(end);
%pop_fr = sum(spike_cnt,2)/64/t_ds(end);

inds = {[1 2 5 6 9 10] [3 4 7 8 11 12]};
figure
hold on
bar(1, mean(pop_fr(inds{1})), 'w', 'edgecolor', 'k')
scatter(1+randn([1,6])/10, pop_fr(inds{1}), 15, 'k', 'fill')
line([1 1], mean(pop_fr(inds{1}))+[-std(pop_fr(inds{1})) std(pop_fr(inds{1}))], 'color', 'k')

bar(2, mean(pop_fr(inds{2})), 'w', 'edgecolor', colors(2,:))
scatter(2+randn([1,6])/10, pop_fr(inds{2}), 15, colors(2,:), 'fill')
line([2 2], mean(pop_fr(inds{2}))+[-std(pop_fr(inds{2})) std(pop_fr(inds{2}))], 'color', colors(2,:))
hold off

xlim([0.5 2.5])
ylabel('Mean Channel Firing Rate (Hz)')
%set(gca, 'ytick',[0 3])
set(gca, 'xtickLabel',{'Control' , 'RTT'})
box off
[h p] = ttest2(pop_fr(inds{1}),pop_fr(inds{2}))
[p h] = ranksum(pop_fr(inds{1}),pop_fr(inds{2}))
%nice_figure(gcf, [fig_folder '5_rettpop_fr'],[2 3])


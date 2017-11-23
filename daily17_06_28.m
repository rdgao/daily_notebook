% analysis and plotting for the pharmacology experiments
cd('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids')
load aggregate.mat
load names.mat
%%
pharm_days = [];
for n = 1:length(dates)
    if strfind(dates(n).name,'Drugs')        
        disp(dates(n-1).name)
        disp(dates(n).name)
        disp(n)
        pharm_days = [pharm_days n];
    end
end
drug_labels = {{'CT' 'CT' 'AP5+CNQX' 'Bicuculline'}, ...
    {'CT' 'Baclofen' 'AP5+CNQX' 'Bicuculline'},...
    {'CT' 'Muscimol' 'AP5' 'CNQX'}};

%% plot kernels
figure
comps = 3;
for i=1:8
    subplot(2,4,i)
    plot((-500:2500)/1000, kernels{pharm_days(comps)-1}{i+4}, 'color', [0 0 0 0.5], 'linewidth',1)   
    hold on
    try
        plot((-500:2500)/1000, kernels{pharm_days(comps)}{i+4}, 'color', [1 0 0 0.5],'linewidth',1)
    end
    hold off
    title(drug_labels{comps}{mod(i-1,4)+1})
    xlim([-500 2500]/1000)
end

%% plot mean/std kernel
figure
colors = get(gca,'colororder');
comps = 2;
for i=1:8
    subplot(2,4,i)
    plot_filled((-500:2500)/1000, mean(kernels{pharm_days(comps)-1}{i+4},2), std(kernels{pharm_days(comps)-1}{i+4},1,2), 'k')
    hold on
    try
        plot_filled((-500:2500)/1000, mean(kernels{pharm_days(comps)}{i+4},2), std(kernels{pharm_days(comps)}{i+4},1,2), colors(2,:))
    end
    hold off
    title(drug_labels{comps}{mod(i-1,4)+1})
    xlim([-500 2500]/1000)
end
axes_h = gca;
legend(axes_h.Children([1 3]), {'Pre' 'Post'})
subplot(2,4,5)
xlabel('Time (s)')
ylabel('Population Spikes')

%% ----- figure S7
comp_sel = [2 3 2 2];
cond_sel = [1 6 2 3];
to_plot = {kernels{pharm_days(comps)-1}{i+4}};
figure
for i=1:4
    subplot(1,4,i)
    comps = comp_sel(i);
    j = cond_sel(i);
    plot_filled((-500:2500)/1000, mean(kernels{pharm_days(comps)-1}{j+4},2), std(kernels{pharm_days(comps)-1}{j+4},1,2), 'k')
    hold on
    try
        plot_filled((-500:2500)/1000, mean(kernels{pharm_days(comps)}{j+4},2), std(kernels{pharm_days(comps)}{j+4},1,2), colors(2,:))
    end
    hold off
    title(drug_labels{comps}{mod(j-1,4)+1})
    xlim([-500 2500]/1000)
end
subplot(1,4,1)
axes_h = gca;
legend(axes_h.Children([3 1]), {'Pre' 'Post'})
xlabel('Time (s)')
ylabel('Population Spikes')
fig_folder =  '~/Dropbox/Research/Reports/Muotri/CorticoidFigs/';
%%
nice_figure(gcf, [fig_folder 'S7_pharm'],[12 3.5]);
%% plot example LFP kernels
figure
colors = get(gca,'colororder');
comps = 1;
chan = [14 28 14 14 14 14 14 30];
for i=1:8
    subplot(2,4,i)
    plot_filled((-500:2500)/1000, mean(LFP_ker_smo{pharm_days(comps)-1}{i+4}(:,:,chan(i)),2), std(LFP_ker_smo{pharm_days(comps)-1}{i+4}(:,:,chan(i)),1,2), colors(1,:))
    hold on
    try
        plot_filled((-500:2500)/1000, mean(LFP_ker_smo{pharm_days(comps)}{i+4}(:,:,chan(i)),2), std(LFP_ker_smo{pharm_days(comps)}{i+4}(:,:,chan(i)),1,2), colors(2,:))
    end
    hold off
    title(drug_labels{comps}{mod(i-1,4)+1})
    xlim([-500 2500]/1000)
end
axes_h = gca;
legend(axes_h.Children([3 1]), {'CT' 'Drug'})

%% plot LFP kernel
well = 6;
figure
for chan = 1:64
    subplot(8,8,chan)
%     plot((-500:2500)/1000,LFP_ker_smo{pharm_days(comps)-1}{well}(:,:,chan), 'linewidth',1, 'color', [0 0 0 0.2]);
%     hold on
%     plot((-500:2500)/1000,LFP_ker_smo{pharm_days(comps)}{well}(:,:,chan), 'linewidth',1, 'color', [1 0 0 0.2]);
%     hold off
    plot_filled((-500:2500)/1000, mean(LFP_ker_smo{pharm_days(comps)-1}{well}(:,:,chan),2), std(LFP_ker_smo{pharm_days(comps)-1}{well}(:,:,chan),1,2), colors(1,:))
    hold on
    plot_filled((-500:2500)/1000, mean(LFP_ker_smo{pharm_days(comps)}{well}(:,:,chan),2), std(LFP_ker_smo{pharm_days(comps)}{well}(:,:,chan),1,2), colors(2,:))
    hold off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([-500 2500]/1000)
end

%% plot PSDs
comps = 1;
well = 12;
fs = 1000;
figure
for chan = 1:64
    subplot(8,8,chan)
    loglog(0:0.5:fs/2, PSDw{pharm_days(comps)-1}{well}(:,chan), 'k')        
    hold on
    loglog(0:0.5:fs/2, PSDw{pharm_days(comps)}{well}(:,chan), 'r')   
    hold off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([0 400])
end
%% plot meanPSD
fs = 1000;
comps = 1;
well = 12;
figure
loglog(0:0.5:fs/2, mean(PSDw{pharm_days(comps)-1}{well},2))
hold on
loglog(0:0.5:fs/2, mean(PSDw{pharm_days(comps)}{well},2))
hold off

%% plot mean spiking change
comps = 1;
FR = [sum(nws_smo{pharm_days(comps)-1}(:,5:end));sum(nws_smo{pharm_days(comps)}(:,5:end))];
figure
colors = get(gca,'colororder');
h1 = plot(FR(:,[1 2 5 6]), 'ok-');
hold on
h2 = plot(FR(:,[3 7]), 'o-', 'color', colors(2,:));
h3 = plot(FR(:,[4 8]), 'o-', 'color', colors(1,:));
hold off
xlim([0.5 2.5])
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'Before' 'After'})
ylabel('Total Spike Count')
legend([h1(1) h2(1) h3(1)], {'Ctrl', 'CNQX+AP5', 'Bicu'}, 'Location', 'Best')


% visualizing spike components from w/ vs. w/o bicuculline (well 12)
% also need to grab the network spikes that actually got captured
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/Spikes/components/';
load nwaydecomp_32root.mat
%%
n_comps = length(nwaycomp.comp);
tr_per_rec = 40;
n_trials = length(nwaycomp.comp{1}{3});
%%
figure
hold on
for i=1:n_comps
    plot(1+(1:n_trials)/tr_per_rec, nwaycomp.comp{i}{3}./max(nwaycomp.comp{i}{3}), 'linewidth',1)
end
hold off
legend(num2str((1:n_comps)'))
xlabel('Recording Number')
ylabel('Component Weight (Normed)')
%% looking at component congruence
figure
subplot(1,3,1)
imagesc(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrall)
colorbar()
title('Congr All')

subplot(1,3,2)
imagesc(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrglobmin)
colorbar()
title('Congr Global Min')

subplot(1,3,3)
plot(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrcumul(:,:,2))
title('Congr Culm')
%% visualizing single components
dish_dim = [8 8];
comp = 3;
[max_amp max_chan] = max(nwaycomp.comp{comp}{1});

figure

% trial profile
subplot(2,1,1)
tr_label = ((0:n_trials-1)/tr_per_rec)+1;
%plot(((1:n_trials)/tr_per_rec)+1, nwaycomp.comp{comp}{3}.^2, 'linewidth',1)
bar(tr_label, nwaycomp.comp{comp}{3}.^2, 'facecolor', 'k')
%plot(((1:n_trials)/tr_per_rec)+1, nwaycomp.comp{comp}{3}./max(nwaycomp.comp{comp}{3}), 'linewidth',1)
xlabel('Recording Number (Days)')
xlim([min(tr_label), max(tr_label)])

% amplitude profile
subplot(2,3,4)
imagesc(reshape(nwaycomp.comp{comp}{1}, dish_dim)')
title('Amplitude Profile')
colorbar
axis off

% frequency profile
subplot(2,3,5)
plot(nwaycomp.freq, nwaycomp.comp{comp}{2})
xlabel('Frequency (Hz)')
ylabel('Component Strength')
title('Frequency Profile')

% phase/time profile
subplot(2,3,6)
[temp max_freq] = max(nwaycomp.comp{comp}{2}); % get max freq
%imagesc(reshape(nwaycomp.comp{comp}{4}(:,max_freq), dish_dim)) %spacefsp
imagesc(reshape(nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan)', dish_dim)) %spacetime
colorbar
caxis([-0.002 0.002])
axis off
title('Time Profile')

[~, amp_sorted] = sort(nwaycomp.comp{comp}{1},'descend');
plot_mea_grid([],[],nwaycomp.comp{comp}{1}, (nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan))*1e3, [], amp_sorted(1:10))



%% grab the time indices again to see what the network spike looks like per trial
load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/aggregate.mat', 'T', 'nws_smo');
%%
well = 12;
fs = 12500;
fs_ds=1000;
% partition parameters
num_part = 40;
part_len = 10; %seconds

recs = [25, 26];
nws_stack = {};
for i=1:length(recs)
    len = T{recs(i)}*fs;
    part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
    tr_info = [part_inds' part_inds'+(part_len-0.5)*fs-1 zeros(num_part,1)];
    SI = round(tr_info(:,1)/fs*fs_ds)+1;
    for j=1:size(tr_info,1)
        nws_stack = [nws_stack nws_smo{recs(i)}(SI(j):(SI(j)+(part_len-0.1)*fs_ds)-1,well)];
    end    
end
nws_stack = cell2mat(nws_stack);
%%
comp = 3;
thresh = 0.15;
figure
plot(nws_stack(:,nwaycomp.comp{comp}{3}./max(nwaycomp.comp{comp}{3})>thresh))
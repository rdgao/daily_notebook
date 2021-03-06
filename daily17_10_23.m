%% I use this code quite a lot so probably should turn them into functions

% looking at the 10component solution
%load '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/LFP/all_lfp_FC/comps/nwaydecomp_set10comps.mat'
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/LFP/good_lfp_FC/'
load '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/LFP/good_lfp_FC/nwaydecomp_set10comp.mat';
load trial_indices
%%
%cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/LFP/good_lfp_FC/'
n_comps = length(nwaycomp.comp);
tr_per_rec = 24; %24;% 16;
n_trials = length(nwaycomp.comp{1}{3});
if exist('select_chans.mat')
    load select_chans.mat
    tr_label = 1+((find(tr_include)-1)/tr_per_rec);
else
    tr_label = 1+(1:n_trials)/tr_per_rec;
end
%% looking at component congruence
comp_names = {'SpatialAmp' 'Freq' 'Trial' 'Time' 'magicD'};
figure
subplot(1,3,1)
imagesc(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrall)
colorbar()
title('Congr All')
xticks(1:5)
xticklabels(comp_names)

subplot(1,3,2)
imagesc(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrglobmin)
colorbar()
title('Congr Global Min')
xticks(1:5)
xticklabels(comp_names)

subplot(1,3,3)
plot(nwaycomp.splitrelstat.allrandomstatfull{n_comps}.congrcumul(:,:,1)) % last index is feature
title('Congr Culm')
legend(num2str((1:n_comps)'))

%% plot all trial profiles and correlation
tr_profs = zeros(length(nwaycomp.comp{1}{3}), n_comps);
profs_agg = {[],[],[],[]};
for i=1:n_comps
    tr_profs(:,i) = nwaycomp.comp{i}{3}./max(nwaycomp.comp{i}{3});
    for j=1:4
       profs_agg{j}(:,i) = nwaycomp.comp{i}{j};
    end
end
figure
subplot(2,2,1)
plot(tr_label, tr_profs, '.-')
legend(num2str((1:n_comps)'))
xlabel('Recording Number')
ylabel('Component Weight (Normed)')
xlim(tr_label([1 end]))

subplot(2,2,2)
imagesc(tr_label, 1:n_comps, tr_profs')
colorbar

for j=1:4
    subplot(2,4,j+4)
    CC = corr(profs_agg{j});
    CC(find(eye(n_comps)))=0;
    imagesc(CC)
    colorbar
    title(comp_names{j})
    axis off
end

%% visualizing single components
%close all
dish_dim = [8 8];
comp = 2;
%max_chan = 9;
[max_amp max_chan] = max(nwaycomp.comp{comp}{1});

%close all
figure

% trial profile
subplot(2,1,1)
bar(tr_label, nwaycomp.comp{comp}{3}./max(nwaycomp.comp{comp}{3}), 'facecolor', 'k')
xlim([1 43])
xlabel('Recording Day')
title(sprintf('Component #%i, Trial Profile', comp))

% amplitude profile
subplot(2,3,4)
imagesc(reshape(nwaycomp.comp{comp}{1}, dish_dim)')
axis off
title('Amplitude Profile')
colorbar

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
%imagesc(reshape(nwaycomp.comp{comp}{4}, dish_dim)') %spacetime
imagesc(reshape(nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan), dish_dim)') %spacetime
colorbar
axis off
title('Time Profile')
%caxis([-0.1 0.1])
[~, amp_sorted] = sort(nwaycomp.comp{comp}{1},'descend');
%plot_mea_grid([],[],nwaycomp.comp{comp}{1}, nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan), [], amp_sorted(1:20))
%%
% plot profiles of all components
thresh_by_amp = 1; % 1 = dont plot channels that dont meet threshold
amp_thresh = 0.25; % threshold as a proportion of strongest channel
grid_setting = [2,5];
for comp=1:n_comps
    if thresh_by_amp
        mask = nan(prod(dish_dim),1);
        [maskv, maskind] = sort(nwaycomp.comp{comp}{1},'descend');
        mask(maskind(maskv>maskv(1)*amp_thresh))=1;
    else
        mask = ones(prod(dish_dim),1);
    end
    
    figure(1)
    subplot(grid_setting(1),grid_setting(2),comp)
    %pcolor(reshape(nwaycomp.comp{comp}{1}.*mask, dish_dim)')    
    imagesc(reshape(nwaycomp.comp{comp}{1}.*mask, dish_dim)')
    %C = colormap;
    %colormap([[0 0 0]; C])
    axis square
    axis off
    
    figure(2)
    subplot(grid_setting(1),grid_setting(2),comp)
    plot(nwaycomp.freq, nwaycomp.comp{comp}{2})
    YL = ylim;
    ylim([0 YL(2)])
    
    figure(3)
    subplot(grid_setting(1),grid_setting(2),comp)
    plot(tr_label, tr_profs(:,comp), '.k')
    xlim(tr_label([1 end]))

    figure(4)
    subplot(grid_setting(1),grid_setting(2),comp)
    [temp max_freq] = max(nwaycomp.comp{comp}{2}); % get max freq
    imagesc(1000*reshape((nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan)).*mask, dish_dim)') %spacetime
    axis square
    axis off
    %C = colormap;
    %colormap([[0 0 0]; C])
    colorbar
end
temp = [1 2 3 4];
for f=1:4
    figure(f)
    subplot(grid_setting(1),grid_setting(2),1)
    title(comp_names{temp(f)})
end
%%
% go into a specific recording to see the time series for a channel
cd '/Volumes/My Passport for Mac/Dish/CTC'
F = dir('CTC_*');
%reorder folders because 2017 recordings are in front
for f=1:length(F)
    if strcmp(F(f).name,'CTC_073116')
        offset = f;
    end
end
F = [F(offset:end);F(1:offset-1)];
date = 25;
disp(F(date).name)
load([F(date).name '/LFP_Sp.mat']);
%% plot trialed data for strongest channel at that comp
well=5;
chan = max_chan; %max_chan;
%chan = 19;

lp_freq = 20;
B = fir1(11./lp_freq*fs_ds, lp_freq/(fs_ds/2));
bp_freq = [0.1 20];
BP = fir1(5./bp_freq(1)*fs_ds, lp_freq/(fs_ds/2));

%data_filt = filter(B,1,LFP{well}(:, chan));
data_filt = filter(B,1,LFP{well}(:, chan));

tr_data = zeros((tr_all{date}(1,2)-tr_all{date}(1,1)+1),tr_per_rec);
for tr=1:tr_per_rec
    %tr_data(:,tr) = LFP{well}(tr_all{date}(tr,1):tr_all{date}(tr,2), chan);    
    tr_data(:,tr) = data_filt(tr_all{date}(tr,1):tr_all{date}(tr,2));
end
plot_tight(tr_data, [4,tr_per_rec/4], (1:size(tr_data,1))/1000, [min(min(tr_data)) max(max(tr_data))])

%% plot a trial filtered
tr = 15;
top_N = 5;
[~, amp_sorted] = sort(nwaycomp.comp{comp}{1},'descend');
figure
plot((1:size(tr_data,1))/fs_ds, filter(B,1,LFP{well}(tr_all{date}(tr,1):tr_all{date}(tr,2), amp_sorted(1:top_N))));
legend(num2str(amp_sorted(1:top_N)))


%% components look like garbage because some recordings are broken.
% go through data to check data quality visually and by range/standard deviation
%% CHECK all corticoid data
% compute fourier coefficients for all recordings
fs = 1000;

% handle file pathing
cd '/Volumes/My Passport for Mac/Dish/CTC/';
F = dir('CTC_*');
%reorder folders because 2017 recordings are in front
for f=1:length(F)
    if strcmp(F(f).name,'CTC_073116')
        offset = f;
    end
end
F = [F(offset:end);F(1:offset-1)];

%% traversing and compute
num_well = 8;
num_chan = 64;

range_all = zeros(length(F), num_well, num_chan);
std_all = zeros(length(F), num_well, num_chan);

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    load LFP_Sp.mat LFP spikes spike_cnt t_s t_ds
    
    %figure
    for well=5:5+num_well-1
        disp(sprintf('Well %i', well));
        
        %skip well if there's no LFP
        if isempty(LFP{well})            
            disp(sprintf('Well %i skipped, no data.', well))
            continue
        end
              
        data = LFP{well}';
        [numchan, len] = size(data);
                
        %subplot(num_well,1,well-4)
        %plot(t_ds, data')
        %title(sprintf('Well %i', well))
        
        range_all(f,well-4,:) = range(data,2);
        std_all(f,well-4,:) = std(data,1,2);
    end
    %nice_figure(gcf, ['../LFP_ts/' F(f).name], [12, 24])    
    %close
    
    cd ..
end
save data_qual.mat range_all std_all

%% visually examine saved plots (in LFP_ts) and enter indexing for well 5
quality = zeros(42,1);
for i=1:42
    quality(i) = input([F(i).name ': ']);
end
%% now grab only the trials from recordings with good data
% quality_5 is the vector for well 5, which was named "quality" above
cd '~/Documents/Data/Muotri/Pri_Corticoids/SPACE/LFP/all_lfp_FC/'
fc_all = load('fourfull.mat');
%% 
good_recs = find(quality_5==1);
tr_include = zeros(1,size(fourier,3));
for rec = 1:length(good_recs)
    tr_include((good_recs(rec)-1)*16+(1:16))=1;
end
disp(all(unique(ceil(find(tr_include)/16))==good_recs')); %should be 1
fourier = fc_all.fourier(:,:,find(tr_include),:); % grab good indices

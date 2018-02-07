% grabbing all the spiking data from harddrive and saving it to disk so i
% don't have to lug this thing around all the time to analyze the dish
% spike data

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
%%
%% traversing and compute
T = zeros(1,length(F));
spiketrains = cell(1,length(F));
spikeshapes = cell(1,length(F));
spikecounts = cell(1,length(F));
names = cell(1,length(F));

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)       
    load LFP_Sp.mat spikes spike_cnt t_s %spike_shape
    T(f) = t_s(end);
    names{f} = F(f).name(5:end);
    spiketrains{f} = spikes;
    %spikeshapes{f} = spike_shape;
    spikecounts{f} = spike_cnt;
    cd ..
end
%%
save spikes_aggregate.mat T names spiketrains spikecounts
%save spikes_aggregate.mat T names spiketrains spikeshapes spikecounts -v7.3
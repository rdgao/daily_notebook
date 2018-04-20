%%% analyzing Spencer's mouse tracking data
% import data from CSV
% import the whole sheet as a variable, pad empty with NANs, then process each matrix to grab
% the actual data


data_raw = {WT_M, WT_F, MT_M, MT_F};
data_proc = {};
data_labels = {'WT_M', 'WT_F', 'MT_M', 'MT_F'};
%% convert to array and save out
for d = 1:length(data_raw)
    data = data_raw{d};
    start_cols = 1:4:size(data,2)
    dataset = {};
    for col = 1:length(start_cols)
        pos = data(:,start_cols(col):start_cols(col)+2);
        pos = pos(~isnan(pos(:,1)),:);
        dataset{col} = pos;
    end
    data_proc{d} = dataset;
end
save mouse_position.mat data_proc data_labels
%% analyze position
load mouse_position.mat
inner = [3,6];
for ds = 1:length(data_proc)
    for ms = 1:length(data_proc{ds})
        coors = data_proc{ds}{ms}(:,2:3);
        inner_x = coors(:,1)>=inner(1) & coors(:,1)<=inner(2);
        inner_y = coors(:,2)>=inner(1) & coors(:,2)<=inner(2);
        IC = find(inner_x.*inner_y);
        OC = setdiff(1:length(coors),IC);
        
        durations = [diff(data_proc{ds}{ms}(:,1)); 0];
        
        ms_times{ds}(ms,1) = sum(durations(IC));
        ms_times{ds}(ms,2) = sum(durations(OC));                               
    end
end
%% plotting
concat_times = cell2mat(ms_times');
x = concat_times(:,1)./sum(concat_times,2);
figure
boxplot(x, [1*ones(1,9),2*ones(1,9),3*ones(1,10),4*ones(1,10)])

hold on
for i=1:4
plot(i, ms_times{i}(:,1)./sum(ms_times{i},2), 'ro')
end
hold off
xticklabels(data_labels)
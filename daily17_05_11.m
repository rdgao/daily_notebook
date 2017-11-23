cd('/Volumes/My Passport for Mac/Dish/CDKL5_CTC')
%% get recording dates from folder string
F = dir('*CTC*');
F = [F(18:end); F(1:17)];
for i=1:length(F)
    dates{i} = F(i).name(end-5:end);
    dayVec(i) = datenum(str2num(dates{i}(5:6)), str2num(dates{i}(1:2)), str2num(dates{i}(3:4)));
end
dayVec = dayVec-dayVec(1)+1;
%dayVec = recDays/7+7; %use weeks instead

%% analyzing aggregate CDKL5 data
wells = 1:12;

for well = wells
    for rec = 1:length(kernels)
        n_events = size(kernels{rec}{well},2);
        
        
        
        % mean and std of inter event latency
        ker_latM(rec,well) = mean(diff(kerTimes{rec}{well}));
        ker_latS(rec,well) = std(diff(kerTimes{rec}{well}));
        
        % event count
        ker_cnt(rec, well) = n_events;
        
        % mean self correlation
        if n_events<2
            %only 1 event, don't compute
            self_cor(rec,well) = NaN;
            self_cor_LFP(rec,well,:) = NaN;
        else
            %more than 1 event, compute pair-wise correlation and average
            self_cor(rec,well) = sum(sum(triu(corr(kernels{rec}{well}),1)))/nchoosek(n_events,2);
            
            % LFP kernel correlation
            for chan = 1:64
                self_cor_LFP(rec,well, chan) = sum(sum(triu(corr(LFP_ker_smo{rec}{well}(:,:,chan)),1)))/nchoosek(n_events,2);                
            end
        end
        
        % total spiking under the curve
        ker_totalM(rec, well) = mean(sum(kernels{rec}{well}));
        
        % initial peak amplitude
        peakAmp = max(kernels{rec}{well});
        ker_peakM(rec, well) = mean(peakAmp);
        ker_peakS(rec, well) = std(peakAmp);                      
    end    
end
%% cycling through days
conds = {[1 2 5 6 9 10], [3 4 7 8 11 12]};
wells = 1:12;
figure
t = (-500:2500)/1000;
for day = 1:length(kernels)
    for well = wells
        subplot(3,4,well)
        plot(kernels{day}{well}, 'color', [0 0 0 0.3], 'linewidth',1)
        %xlim(t([1 end]))
        
%         loglog(0:0.5:500, PSDw{day}{well}, 'color', [0 0 0 0.1], 'linewidth', 1);
%         xlim([0 100])
    end
    subplot(3,4,1)
    title(dates{day})
    pause
end

%% overlap all days
figure
t = (-500:2500)/1000;
for day = 1:length(kernels)
    for well = wells
        subplot(3,4,well)
        hold on
        plot(t, mean(kernels{day}{well},2), 'linewidth',1)
        xlim(t([1 end]))
        ylim([0 5])
    end   
end
hold off
%legend(dates)

%% kernel aggregate
nsubs = zeros(length(kernels),12);
peakAmp = zeros(length(kernels),12);
FWHM = zeros(length(kernels),12);
total_spikes = zeros(length(kernels),12);
for day = 1:length(kernels)
    disp(dates{day})
    for well = wells
        [nsubs_, peakAmp_,FWHM_, subpks_, subpksT_] = kernel_pkfinding(kernels{day}{well});
        nsubs(day, well) = mean(nsubs_);
        peakAmp(day, well) = mean(peakAmp_);
        FWHM(day, well) = mean(FWHM_);
        total_spikes(day, well) = mean(sum(kernels{day}{well},1));
    end
end

%% mean kernel period
figure
plot(dayVec,ker_latM(:,conds{1}), '-o', 'color', [0 0 0 0.4])
hold on
plot(dayVec,ker_latM(:,conds{2}), '-o', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('Event Latency')

%% mean % of subpeaks
figure
plot(dayVec,nsubs(:,conds{1}), '-o', 'color', [0 0 0 0.4])
hold on
plot(dayVec,nsubs(:,conds{2}), '-o', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('# of subpeaks')

%% mean first peak height
figure
plot(dayVec,peakAmp(:,conds{1}), '-o', 'color', [0 0 0 0.4])
hold on
plot(dayVec,peakAmp(:,conds{2}), '-o', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('Peak 1 Amplitude')

%% mean first peak width
figure
plot(dayVec,FWHM(:,conds{1}), '-o', 'color', [0 0 0 0.4])
hold on
plot(dayVec,FWHM(:,conds{2}), '-o', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('Peak 1 Width')

%% mean total spikes
figure
plot(dayVec,total_spikes(:,conds{1}), '-o', 'color', [0 0 0 0.4])
hold on
plot(dayVec,total_spikes(:,conds{2}), '-o', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('Total spikes in event')
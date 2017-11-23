% looking at drug induced differences in slope in organoids LFP
% first get slope matrix from corticoid_flat_analysis

%pairs of drug recordings
drug_idx = {[25 26], [34 35], [39 40]};
drug_labels = {{'ctrl','ctrl','-GLUT','-GABA'},...
                {'ctrl','+GABAB','-GLUT','-GABA'},...
                {'ctrl','-GABAA','-NMDA','-AMPA'}};
            
slp_label = {'1-10 Hz' '10-30 Hz' '30-80 Hz'};

figure
for slp_range = 1:3
    for drug_cond = 1:3
        disp(dates(drug_idx{drug_cond}(2)).name)
        subplot(3,3,(slp_range-1)*3+drug_cond)
        slp_diff = squeeze(slopes(drug_idx{drug_cond}(1),:,:,slp_range)-slopes(drug_idx{drug_cond}(2),:,:,slp_range));        
        boxplot(slp_diff(:,:)')
        set(gca,'xtickLabel', drug_labels{drug_cond})
        line([0 9], [0 0], 'color', 'k')        
        [h pv] = ttest(slp_diff');
        yl=ylim();
        hold on
        plot(find(pv<0.05), yl(2)*ones(1,length(find(pv<0.05))), 'ok', 'markersize', 15)
        plot(find(pv<0.01), yl(2)*ones(1,length(find(pv<0.01))), 'or')
        hold off
        title(slp_label{slp_range})
    end
end

%%
% effect of pharm on oscillation 
figure
for drug_cond = 1:3
    subplot(1,3,drug_cond)    
    disp(dates(drug_idx{drug_cond}(2)).name)
    slp_diff=squeeze(osc_power_foof(drug_idx{drug_cond}(2),:,:)-osc_power_foof(drug_idx{drug_cond}(1),:,:));
    boxplot(slp_diff(:,:)')
    set(gca,'xtickLabel', drug_labels{drug_cond})
    line([0 9], [0 0], 'color', 'k')
    [h pv] = ttest(slp_diff');
    yl=ylim();
    hold on
    plot(find(pv<0.05), yl(2)*ones(1,length(find(pv<0.05))), 'ok', 'markersize', 15)
    plot(find(pv<0.01), yl(2)*ones(1,length(find(pv<0.01))), 'or')
    hold off    
end

%% progression of slope over time
for slp_range=1:3
    subplot(3,1,slp_range)
    plot(dayVec, median(slopes(recs,:,:,slp_range),3),'.-')
    title(slp_label{slp_range})
end
xlabel('Weeks')
ylabel('Slope')
% PCA on firing rate
fs = 12500;
fs_ds = 1000;
bsp = binarize_spikes(ceil(t_ds(end)), fs,spikes,fs_ds);
%%
well = 5;
FR = (tsmovavg(squeeze(bsp(well,:,:)),'e',50));
%%
figure
plot(sum(FR,1))
%%
T0 = 144000;
T = T0:(T0+4000-1);
[coef, sc, lat, tsq, exp] = pca(FR(:,T));
figure
bar(exp)
xlabel('# of PC')
ylabel('variance explained')
xlim([1,10])
figure
plot(coef(:,1:2))
legend({'PC1' 'PC2'})
figure
scatter(coef(:,1),coef(:,2),1, T-T(1))
xlabel('PC1')
ylabel('PC2')
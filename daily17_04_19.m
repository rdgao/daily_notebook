%% trying different filter length for spiking filtering, 
cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017
load well_5.mat
%%
fs = 12500;
chan = 10;
f_o = [20, 50, 100, 200];
for i=1:4
    b = fir1(f_o(i),[300 3000]/(fs/2));
    filtered(:,i) = filtfilt(b,1,MEA(:,chan));    
end
figure 
plot(t, filtered, 'linewidth',1)

%% and causal vs. non-causal
figure
chan = 2;
f_o = 50;
b = fir1(f_o,[300 3000]/(fs/2));
f = filter(b,1,MEA(:,chan));
ff = filtfilt(b,1,MEA(:,chan));
plot(ff, 'linewidth', 1)
hold on
plot(f(f_o/2+1:end), 'linewidth', 1)
hold off


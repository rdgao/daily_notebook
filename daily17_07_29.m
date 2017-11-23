% playing with 4k recording
load('/Users/rgao/Documents/Data/Muotri/4k/3Braintest.mat')
%%
V = whos;
V = V(2:4097);
%%
Vrng = zeros(4096,2);
for i=1:4096
    Vrng(i,1) = min(eval(V(i).name));
    Vrng(i,2) = max(eval(V(i).name));
end
%%
chan = 1029;
fs = 1.7856e4;
[F A T] = spectrogram(double(eval(V(chan).name)),fs*2,fs,fs*2,fs);
loglog(A, mean((abs(F)),2))
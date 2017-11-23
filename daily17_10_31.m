% making movie for organoid
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416';
load LFP_Sp.mat spikes t_s
fs = 12500;
fs_ds = 1000;
well=5;
bsp = squeeze(binarize_spikes(t_s(end), fs, spikes(well,:), fs_ds));
%%
% smooth
conv_win = gausswin(101);
numchan = size(bsp,1);
bsp_smo = zeros(size(bsp));
for chan = 1:numchan
    bsp_smo(chan,:) = conv(bsp(chan,:),conv_win, 'same');
end
%%
bsp_sq = reshape(bsp_smo, [8,8,size(bsp_smo,2)])./max(max(bsp_smo));
%bsp_sq = round(reshape(bsp_smo, [8,8,size(bsp_smo,2)]));
%%
start_ind = 108000;
end_ind = start_ind+5*fs_ds;
figure
colormap(gray)
pause
for i=start_ind:4:end_ind    
    imagesc(imgaussfilt(bsp_sq(:,:,i),0.5))
    caxis([0 1])
    title(sprintf('%.3f',(i-start_ind)/fs_ds))
    pause(0.01)
end
 

%% WTF?
start_ind = 109000;
end_ind = start_ind+5*fs_ds;
video = VideoWriter('dish', 'MPEG-4'); %create the video object
video.FrameRate = 30;
open(video); %open the file for writing
for ii=1:(end_ind-start_ind) %where N is the number of images  
  writeVideo(video,bsp_sq(:,:,start_ind+ii-1)); %write the image to file
end
close(video); %close the file
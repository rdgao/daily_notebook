% write datafile as binary for spike sorting in klusta
%datafolder = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_031017/';
datafolder = '/Volumes/My Passport for Mac/Dish/CTC/CTC_121716/';
datafile = 'well_12.mat';
load([datafolder datafile])

%% save data to binary .dat format
volt_conv = 5.51e-8;
output_file = ['binary_' datafile(1:end-4) '.dat'];
fileID = fopen(output_file,'w');
flat_MEA = reshape(MEA',1,[])/volt_conv;

disp('Writing..')
count = fwrite(fileID, int16(flat_MEA), 'int16')
fclose(fileID);
disp('Done.')

%% read in data file to make sure it saved correctly
input_file = ['binary_well_5.dat'];
%input_file = 'hybrid_10sec.dat';
fileID = fopen(input_file,'r');
x = fread(fileID, inf, 'int16');
fclose(fileID);
disp('Done.')

%% take klusta cells and save out spike times in mat cell format
kwik_file = '';


%% making probe geometry
clc
for i=1:64
    x = mod(i-1,8)*50;
    y = floor((i-1)/8)*50;
    disp(sprintf('%i: (%i, %i),',i-1,x,y))
end
%% making adjacency graph
clc

%% making disconnected probe .prb file
clc
for i=0:63
    x = mod(i,8)*50;
    y = floor((i)/8)*50;
    disp(sprintf('%i: {''channels'': [%i], ''graph'':[], ''geometry'':{%i:(%i,%i)} },', i, i,i, x, y));
end



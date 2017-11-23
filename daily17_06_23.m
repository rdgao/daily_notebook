% analysis and plotting for the pharmacology experiments
%cd('/Volumes/My Passport for Mac/Dish/CTC')
cd('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids')
load names.mat
for n = 1:length(dates)
    if strfind(dates(n).name,'Drugs')
        disp(n)
        disp(dates(n).name)
    end
end
%%
plot_tight(nws_smo{25}, [3 4], [], [0 25])
plot_tight(nws_smo{26}, [3 4], [], [0 25])
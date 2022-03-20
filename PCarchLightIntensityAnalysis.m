function [allInt] = PCarchLightIntensityAnalysis()
if ismac
    cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/M39GADCre_04062020')
else
    cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M39GADCre_04062020');
end
%% params
clear;clc
normFR = 1;
lightInt = [0 18 25 30];%mW/mm^2
%% data
dataSet(1).exps = 4;dataSet(1).clu = 1:3;
dataSet(2).exps = 5;dataSet(2).clu = 1:2;
dataSet(3).exps = 7;dataSet(3).clu = 1:4;
dataSet(4).exps = 7;dataSet(4).clu = 1;
%%
allInt = [];
for i = 1:size(dataSet,2)
    for j = 1:size(dataSet(i).clu,2)
        if ismac
            [exp] = readExpDescFile_mac(pwd); fileindx = dataSet(i).exps;
        else
            [exp] = readExpDescFile(pwd); fileindx = dataSet(i).exps;
        end
        lp = initLightParam(exp,fileindx,1); ...
            lp = readLpData(lp,[-10 10],[1 1 1],dataSet(i).clu(j));
        [~,orig] = odorLightIntensitiesAnalysis(lp,'inhalationonset');
        close;
        if normFR==1
%             allInt = [allInt;(orig-min(orig))./max(orig-min(orig))];
              allInt = [allInt;orig./orig(1)];
        else
            allInt = [allInt;orig];
        end
    end
end
figure;
plot(lightInt,allInt','-o')
hold on;
errorbar(lightInt,mean(allInt),std(allInt)./sqrt(size(allInt,1)),'-ok', 'linewidth',3)
ylabel('Firing rate (Spikes/Sec)','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% xlim([-5 35])
hold off;
function [allInt] = OBchr2LightIntensityAnalysis(normFR)
cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M42GADCre_24062020');
%% experiment with no specific script
%% params
lightInt = [0 18 25 30];%mW/mm^2
%% Experiment list

dataSet(1).exps = 4;dataSet(1).clu = 1;
dataSet(2).exps = 6;dataSet(2).clu = 1;
dataSet(3).exps = 7;dataSet(3).clu = 1;
dataSet(4).exps = 8;dataSet(4).clu = 2;
dataSet(5).exps = 9;dataSet(5).clu = 1;
dataSet(6).exps = 10;dataSet(6).clu = [1, 2];
dataSet(7).exps = 11;dataSet(7).clu = 1;
dataSet(8).exps = 12;dataSet(8).clu = [1,2];
%%
allInt = [];
for i = 1:size(dataSet,2)
    for j = 1:size(dataSet(i).clu,2)
        [exp] = readExpDescFile(pwd); fileindx = dataSet(i).exps;
        lp = initLightParam(exp,fileindx,1); ...
            lp = readLpData(lp,[-10 10],[1 1 1],dataSet(i).clu(j));
        [~,orig] = odorLightIntensitiesAnalysis(lp,'inhalationonset');
        close;
        if normFR==1
%             allInt = [allInt;(orig-min(orig))./max(orig-min(orig))];
              allInt = [allInt;orig./orig(1)];
        else
%             allInt = [allInt;orig];
              allInt = [allInt;(orig-orig(1))./orig(1)];
        end
    end
end
figure;
plot(lightInt,allInt','-o')
hold on;
errorbar(lightInt,mean(allInt),std(allInt)./sqrt(size(allInt,1)),'-ok', 'linewidth',3)
ylabel('Normalized firing rate (Spikes/Sec)','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xlim([-5 35])
hold off;
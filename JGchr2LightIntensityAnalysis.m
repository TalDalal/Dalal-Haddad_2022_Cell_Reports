function [allInt] = JGchr2LightIntensityAnalysis(normFR)
%% PC recordings, while activating JG
%% experiment with no specific script
%% params
lightInt = [0 0.3 .7 1];%mW/mm^2

%% Experiment list

dataSet(1).exps = 11;dataSet(1).clu = 1;dataSet(1).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\ChR2_PG11022020';
dataSet(2).exps = 13;dataSet(2).clu = [1,2,3];dataSet(2).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\ChR2_PG11022020';
dataSet(3).exps = 15;dataSet(3).clu = [1,2];dataSet(3).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\ChR2_PG11022020';
dataSet(4).exps = 16;dataSet(4).clu = [1,2,3];dataSet(4).dir = ...
   'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\ChR2_PG11022020';

allInt = [];
for i = 1:size(dataSet,2)
    for j = 1:size(dataSet(i).clu,2)
        cd(dataSet(i).dir)
        [exp] = readExpDescFile(pwd); fileindx = dataSet(i).exps;
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
ylabel('Normalized firing rate (Spikes/Sec)','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xlim([-.1 1.1])
hold off;
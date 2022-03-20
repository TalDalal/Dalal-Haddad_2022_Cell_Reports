function [allInt] = PCchr2LightIntensityAnalysis(normFR)
cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M41GADCre_10062020');
%% experiment with no specific script
%% params
exps = [1 2 3];
clusters = 1:2;
lightInt = [0 18 25 30];%mW/mm^2
bgLight = [];bgOdor = [];
fgLight = [];fgOdor = [];
%%
for i = 1:length(clusters)
    for j = 1:length(exps)
        [exp] = readExpDescFile(pwd); fileindx = exps(j);
        lp = initLightParam(exp,fileindx,1); lp = readLpData(lp,[-10 10],[1 1 1],clusters(i));
        res = odorLightStimOrResp(lp,'inhalationOnset');
        close;
        %odor+light
        bgLight(i,j) = mean(res(1).BgFg(:,1)); fgLight(i,j) = mean(res(1).BgFg(:,2));
        %odor
        bgOdor(i,j) = mean(res(2).BgFg(:,1)); fgOdor(i,j) = mean(res(2).BgFg(:,2));
    end
end
odorMat = mean(mean(fgOdor-bgOdor));
clus1 = [odorMat fgLight(2,3) fgLight(1,2) fgLight(1,1)];
clus2 = [odorMat fgLight(1,3) fgLight(2,2) fgLight(2,1)];
if normFR==1
%     clus1 = (clus1-min(clus1))./max((clus1-min(clus1)));
%     clus2 = (clus2-min(clus2))./max((clus2-min(clus2)));
        clus1 = clus1./clus1(1);
        clus2 = clus2./clus2(1);
end
allInt = [clus1;clus2];

%% specific script

dataSet(1).exps = 8;dataSet(1).clu = [1,2,4];dataSet(1).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M41GADCre_10062020';
dataSet(2).exps = 19;dataSet(2).clu = 1;dataSet(2).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M41GADCre_10062020';
dataSet(3).exps = 22;dataSet(3).clu = 1;dataSet(3).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M41GADCre_10062020';
dataSet(4).exps = 1;dataSet(4).clu = 2;dataSet(4).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M42GADCre_24062020';
dataSet(5).exps = 2;dataSet(5).clu = 1;dataSet(5).dir = ...
    'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M42GADCre_24062020';


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
xlim([-5 35])
hold off;
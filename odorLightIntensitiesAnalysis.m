function [deltaFR, origEvFR] = odorLightIntensitiesAnalysis(lp,alighnType)
%Description:
% alighnType: OdorOnset or InhalationOnset
%% Extract the data
clrVal = unique(lp.clrVec);
odor = lp.TTLEvents(:,find(lp.clrVec==0)); % baseline: clr==0, only odor

odorLight = [];
for ii = 2:length(clrVal)
    odorLight = [odorLight ; lp.TTLEvents(:,find(lp.clrVec==clrVal(ii)))];
    
end
titArry = {'clr:0','clr:80','clr:150','clr:255'};
EventTime = [odor;odorLight];% create an array of all conditions.


%% plot the PSTH and raster plot for each condition.
if strcmpi(alighnType,'inhalationonset')
    lp.sniffTimes = [lp.inhaleTimes;lp.exhaleTimes];% combine it to 2d matrix.
end
c = 0;Ylim = [];allRes = [];
figure;
for i = 1:2:length(clrVal)*2%run over all conditions.
    c = c + 1;
    subplot(2,3,c)
    
    [res] = labPlotOdor(EventTime(i:i+1,:), lp.t, lp.sniffTimes,  alighnType);% plot the PSTH and raster plot for each condition.
    allRes(c).R = res.R;allRes(c).t = res.t;allRes(c).spikesPerSniff = res.brFrAll;allRes(c).beforeAndAfter = res.odorBeforeAndAfterFr;
    
    % take the ylim from the first figure and apply it to the rest.
    a = gca;
    Ylim = [Ylim get(a,'YLim')];
    title([res.title titArry(c)])
    ylim(minmax(Ylim))
    
end
suptitle([num2str(lp.fileIndx) '_' num2str(lp.mClusterIndx) ', ' alighnType])
%% superposition
if strcmpi(alighnType,'odoronset')
    
    evOdorOnly = mean(allRes(1).beforeAndAfter(:,2));%-mean(allRes(1).beforeAndAfter(:,1)));%evoked fr for clr==0.
    evFRodorLight = [];
    for i = 2:length(clrVal)
%         fg = allRes(i).beforeAndAfter(:,2);
%         bg = allRes(i).beforeAndAfter(:,1);
        evFRodorLight(end+1) = mean(allRes(i).beforeAndAfter(:,2));%-evOdorOnly;%mean(allRes(i).beforeAndAfter(:,1)));%other light int.
    end
    for jj = 1:4
        origEvFR(jj) = mean(allRes(jj).beforeAndAfter(:,2));
    end
    deltaFR = [evOdorOnly evFRodorLight];
    
    subplot(2,3,5)
    plot(clrVal,(deltaFR),'-ok','lineWidth',2)
    xlabel('clr')
    ylabel('evoked fr')
    box off;
elseif strcmpi(alighnType,'inhalationonset')
    bgSniffs = 2:4;
    fgSniffs = 5:8;
    evOdorOnly = mean(mean(allRes(1).spikesPerSniff(:,fgSniffs)));%-mean(mean(allRes(1).spikesPerSniff(:,bgSniffs))));%evoked fr for clr==0
    evFRodorLight = [];
    for i = 2:length(clrVal)
        evFRodorLight(end+1) = mean(mean(allRes(i).spikesPerSniff(:,fgSniffs)))-evOdorOnly;%mean(mean(allRes(i).spikesPerSniff(:,bgSniffs))));%other light int.
    end
    deltaFR = [evOdorOnly evFRodorLight];
    
    for i = 1:size(allRes,2)
        origEvFR(i) = mean(mean(allRes(i).spikesPerSniff(:,fgSniffs)));%;-mean(mean(allRes(i).spikesPerSniff(:,bgSniffs))));
    end
    subplot(2,3,5)
    plot(clrVal,(origEvFR),'-ok','lineWidth',2)
    xlabel('clr')
    ylabel('evoked fr')
    box off;
end

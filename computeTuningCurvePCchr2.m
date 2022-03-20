function [] = computeTuningCurvePCchr2(virusType)
%%% Description
% This function counts for each neuron the number of odor it responded to
% (p<0.05 corrected by the number of odors).
%% load
if strcmpi(virusType, 'chr2')
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_27Neurons.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_27Neurons.mat')
        end
%     if ismac
%         load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_panelOdors.mat')
%     else
%         load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_panelOdors.mat')
%     end
    data = dataStructPCX_ChR2;
    titStr = 'ChR2';
elseif strcmpi(virusType, 'arch')
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataPCarch_6OdorsAllNeurons.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataPCarch_6OdorsAllNeurons.mat')
    end
    data = dataPCarch_6OdorsAllNeurons;
    titStr = 'ArchT';
end
%%
%%%%%%%%%%%%%%%%% sparseness analysis %%%%%%%%%%%%%%%
jitterFactor = 0.1;
%%% build the tuning curve for each neuron
numOdors = 6;c = 0;alphaVal = 0.05/numOdors;
lightResp = [];odorResp = [];
for i = 1:size(data,2)/numOdors
    for j = 1:numOdors
        c = c+1;
        lightResp(i,j) = sum(data(c).odorOnset(1).pVal<=alphaVal);
        odorResp(i,j) = sum(data(c).odorOnset(2).pVal<=alphaVal);
    end
end



%%%%%%%%%%%%%%%%% per neuron %%%%%%%%%%%%%%
neuronsLightResp = sum(lightResp');
neuronsOdorResp = sum(odorResp');

%%% plot by pair

figure;
for i = 1:length(neuronsLightResp)
    plot(1:2,[neuronsLightResp(i) neuronsOdorResp(i)]+jitterFactor*rand(1),'o-','color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none')
    hold on;
end
errorbar([0.75 2.25],[mean(neuronsLightResp) mean(neuronsOdorResp)],[std(neuronsLightResp) std(neuronsLightResp)]./sqrt(length(neuronsLightResp)) ,'ok','MarkerFaceColor','k','MarkerEdgeColor','none')
xlim([0.25 2.75])
box off;
ylabel('# Effective odors','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
[pRank] = signrank(neuronsLightResp, neuronsOdorResp);
disp(['ttets for # effective odors: p = ' num2str(pRank)])
ylim([0 6.2])
title(titStr);

% Plot by delta histogram
responseDelta = neuronsLightResp-neuronsOdorResp;
figure;
meanDelta = mean(responseDelta);
seDelta = std(responseDelta)./sqrt(length(responseDelta));
%for jitter
a = 0.7;
b = 1.3;
r = (b-a).*rand(1,length(responseDelta)) + a;

barwitherr(seDelta, meanDelta)
hold on;
scatter(r.*ones(1, length(responseDelta)), responseDelta, 'k');
set(gca,'linewidth',0.25)
ylabel('\Delta number of eefective odors')
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);


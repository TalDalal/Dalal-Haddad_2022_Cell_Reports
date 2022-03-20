clear;clc
if ismac
    cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2/Neuronexus_firstExp/dataWithSniffs')
else
    cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2\Neuronexus_firstExp\dataWithSniffs')
end
load('exp14_TT1_1.mat')
%%
odorTrials = find(expLFP_SPIKES(1).stimuliTyps==0);
odorLightTrials = find(expLFP_SPIKES(1).stimuliTyps==1);

channel2use = 1;
trialNum = 2;
dataXLIM = [-.5 2];
Fs = expLFP_SPIKES(1).decimateFs;
 currLFP = expLFP_SPIKES(channel2use).LFP(odorTrials(trialNum), :);
 currSpikes = expLFP_SPIKES(1).spikes{odorTrials(trialNum)};
%%

figure;
subplot(3,1,1)
plot(expLFP_SPIKES(1).t, currLFP);
xlim(dataXLIM);
box off;

subplot(3,1,2);
currSniffSignalT = expLFP_SPIKES(1).sniffsData.tSniffSignalByTrial(trialNum).odor;
currSniffSignal = expLFP_SPIKES(1).sniffsData.sniffSignalByTrial(trialNum).odor;
plot(currSniffSignalT, currSniffSignal);
xlim(dataXLIM);
box off;

subplot(3,1,3);
hold all;
for spikeind = 1:length(currSpikes)
    plot([currSpikes(spikeind) currSpikes(spikeind)], [0 1], 'k-')
end
% scatter(currSpikes, ones(1,length(currSpikes)))
xlim(dataXLIM);
box off;hold off;
ylim([-1 2])

figure;
cwt(currLFP, Fs, 'VoicesPerOctave', 32, 'NumOctaves', 10)
xlim(dataXLIM+expLFP_SPIKES(1).t(end));
caxis([300 800]);
box off;
set(gcf, 'Renderer', 'painters')
colormap(othercolor('BuOr_10'))


function [sAllNeuronsOL,sAllNeuronsOdor, numOfNeurons] = decodingAnalysisPCchr2(virusType)
%% Decoding analysis
%%% Description
% This function uses a leave-obe-out decoder to classify the odor identity
% in a test trial based on the average response for each odor, using
% euclidian distance.
%% params
sCalcTime = [0 1.5];% stim duration
odorVec  = [1:4 6 7];
decoderReps = 100; % how many times to run the decoder.
maxTrialNum = 5;
allMean = [];allSEs = [];
bootstrapFlag = 0;
%% load data
%%%%%%%%%%%%%%%%% PC - GC chr2 %%%%%%%%%%%%%%%%%%%%%%%
if strcmpi ('chr2', virusType)

    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_panelOdors.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_panelOdors.mat')
    end

    % combine datasets
    data = dataStructPCX_ChR2;
    %%%%%%%%%%%%%%%%% PC - GC ArchT %%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi ('arch', virusType)
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataPCarch_6OdorsAllNeurons.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataPCarch_6OdorsAllNeurons.mat')
    end
    data = dataPCarch_6OdorsAllNeurons;
    data(7:18) = [];% no valid respiration signal.
end
labels = [];dbOL = [];dbOdor = [];
%% odor+light database
dbOL = [];
for i = 1:size(data,2)% run over all cell-odor pairs
    if ~isempty(data(i).InhOnset)
        fg =[];bg = [];
        for j = 1:size(data(i).InhOnset(1).rasterPlot,2)
            % run over trials.
            spikes = data(i).InhOnset(1).rasterPlot(j).times;% run over trials.

            bg(j) = length(find(spikes >= -sCalcTime(2) ...
                & spikes < sCalcTime(1)))/sCalcTime(2);% baseline

            fg(j) = length(find(spikes >= sCalcTime(1) ...
                & spikes < sCalcTime(2)))/sCalcTime(2); % forground

        end
        labels(end+1) = data(i).odorNum;
        % save the evoked response.
        dbOL(end+1).features = fg-mean(bg);
    end
end
%% odor only database
fg = [];bg = [];dbOdor = [];
for i = 1:size(data,2)% run over all cell-odor pairs
    if ~isempty(data(i).InhOnset)
        fg = [];bg = [];
        for j = 1:size(data(i).InhOnset(2).rasterPlot,2)
            spikes = data(i).InhOnset(2).rasterPlot(j).times;% run over trials.

            bg(j) = length(find(spikes >= -sCalcTime(2) ...
                & spikes < sCalcTime(1)))/sCalcTime(2);% baseline

            fg(j) = length(find(spikes >= sCalcTime(1) ...
                & spikes < sCalcTime(2)))/sCalcTime(2); % forground
        end
        %         labels(end+1) = data(i).odorNum;
        % save the evoked response.
        dbOdor(end+1).features = fg-mean(bg);
    end
end


%% decoder
%% params
sucsVecOL = [];sucsVecOdor = [];sAllNeuronsOL = [];sAllNeuronsOdor = [];
odors = [1 2 3 4 6 7];
jumpVal = 4;
numOfNeurons = [1:jumpVal:size(dbOdor,2)/length(odors)-1 size(dbOdor,2)/length(odors)]; % number of neurons in dataset.
%% run the decoder
for k = numOfNeurons
    for j = 1:decoderReps % run the decoder n times to compute the mean+-sem.
%         if k ~= size(dbOdor,2)/length(odors)
            if bootstrapFlag == 0
                p = randperm((size(dbOdor,2)/length(odors))-1,k); % rand k neurons
            else%bootstrap
                p = randsample(1:(size(dbOdor,2)/length(odors))-1, k, true);
            end
            loc = [];
            %         p = p*length(odors)+1;%%% note the the +1%% convert p from neuron number to cell-odor number (6 odors).
            p = (p*length(odors))-length(odors)+1;
            for l = 1:length(p)
                loc = [loc p(l):p(l)+length(odors)-1]; % concatenate all locations
            end
%         else % for maximal number of neurons.
%             loc = 1:size(dbOdor,2);
%             sOL = generalLeav1OutDecoder_Tal(labels(loc),dbOL(loc),odors,maxTrialNum);
%             sOdor = generalLeav1OutDecoder_Tal(labels(loc),dbOdor(loc),odors,maxTrialNum);
%             sucsVecOL(:,1) = ones(decoderReps,1).*sOL;
%             sucsVecOdor(:,1) = ones(decoderReps,1).*sOdor;
%             break;
%         end
        sOL = generalLeav1OutDecoder_Tal(labels(loc),dbOL(loc),odors,maxTrialNum);
        sOdor = generalLeav1OutDecoder_Tal(labels(loc),dbOdor(loc),odors,maxTrialNum);
        sucsVecOL(j,1) = sOL;
        sucsVecOdor(j,1) = sOdor;
    end
    % save all success rate for each condition, following a decoding using
    % k neurons
    sAllNeuronsOL(:,end+1) = sucsVecOL;
    sAllNeuronsOdor(:,end+1) = sucsVecOdor;
    sucsVecOL = [];
    sucsVecOdor = [];
end

%% plot

figure;
meanOL = mean(sAllNeuronsOL);steOL = std(sAllNeuronsOL);
meanOdor = mean(sAllNeuronsOdor);steOdor = std(sAllNeuronsOdor);

errorbar_patch(numOfNeurons,meanOL,steOL,'b')%% odor+light
hold on;
errorbar_patch(numOfNeurons,meanOdor,steOdor,'k')%% odor
chance = 1/length(odors);
plot([0 numOfNeurons(end)],[chance chance],'k--')
box off;
xlabel('# Neurons', 'FontSize',14)
ylabel('Classification success rate','FontSize',14);
box off;
set(gca,'linewidth',.25)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gcf,'Renderer', 'painter');
%%
% %%%% Decoder CV
% M = mean([(steOL./meanOL)',(steOdor./meanOdor)']);
% ste = std([(steOL./meanOL)',(steOdor./meanOdor)'])./sqrt(length(meanOdor));
% [h,p] = ttest(steOL./meanOL,steOdor./meanOdor)
% figure;
% barwitherr(ste,M)
% box off;
% ylabel('Decoder CV','fontSize',14)
% set(gca,'fontSize',14)
% box off;
% set(gca,'linewidth',1)
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);

% dbPC = dbOdor;
% matOdors = [];c = 0;
% numOdors = 6;numTrials = 5;neuron = [];
% for i = 1:size(dbPC,2)/numOdors
%     for j = 1:numOdors
%         c = c+1;
%        neuron = [neuron; dbPC(c).features(1:numTrials)'];
%     end
%     matOdors(:,i) = neuron;
%     neuron = [];
% end
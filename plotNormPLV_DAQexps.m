clear;clc
% compMethod = 1;%0 = peak detection, 1 = hilbert transform
mouseType = 'ArchT';%'Tbet_NpHR'
%% load
if strcmpi(mouseType,'ArchT')
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data/Awake_data_updated/vsShiftPredictorAnalysis/ArchT')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data\Awake_data_updated\vsShiftPredictorAnalysis\ArchT')
    end
%     load('Anesth_40_70Hz_hilbert_PLV_SFC.mat')
    load('Anesth_40_70Hz_PhaseCorrected_3.mat')%Used for the analyses
%     before the review of cell reports
%     load('Anesth_40_70Hz_PhaseCorrected_ArchT_correctedPowerRatio')
elseif strcmpi(mouseType,'Tbet_NpHR')
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data/Awake_data_updated/vsShiftPredictorAnalysis/Tbet_NpHR')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data\Awake_data_updated\vsShiftPredictorAnalysis\Tbet_NpHR')
    end
%     load('Anesth_40_70Hz_hilbert_PLV_SFC_Tbet.mat')
      load('Anesth_40_70Hz_PhaseCorrected_NpHR_3.mat')
end

%% extract the data
odorPLVCorrected = [];odorLightPLVCorrected = [];
odorFR = [];odorLightFR = [];
shuffledData = [];odorBiasedPLV = [];odorLightBiasedPLV = [];
firingTH = 1;%Hz
minTrialNum = 5;
spectrumOdor = [];spectrumOdorLight = [];
spectrumOdorRelative = [];spectrumOdorLightRelative = [];
relativeGammaPowerOdor = [];relativeGammaPowerOdorLight = [];
%-pValue indicating odor response
odorResponse = nan(1,size(allPLV,2));
odorLightResponse = odorResponse;
FRstructMTC = [];

for i = 1:size(allPLV, 2)
    
    lowSpikeTrialsOdor = length(find(allPLV(i).data.FR.odor < firingTH));
    lowSpikeTrialsOdorLight = length(find(allPLV(i).data.FR.odorLight < firingTH));
    
    if lowSpikeTrialsOdor > minTrialNum | lowSpikeTrialsOdorLight > minTrialNum
        
        fprintf(['Neuron # ' num2str(i) ' was not included in the analysis \n'])
    else
        
        %------ analyze FR -----------
        odorFR(i,:) = ones(1,length(allPLV(i).data.normPLV.odor)).*mean(allPLV(i).data.FR.odor);
        odorLightFR(i,:) = ones(1,length(allPLV(i).data.normPLV.odorLight)).*mean(allPLV(i).data.FR.odorLight);
        
        odorResponse(i) = allPLV(i).data.odorResponsePvalue.odor;
        odorLightResponse(i) = allPLV(i).data.odorResponsePvalue.odorLight;
        
        %------------- For odor response variability---------------
        FRstructMTC(i).odorOnset(1).BgFg = [allPLV(i).data.baselineFR.odorLight' allPLV(i).data.FR.odorLight'];
        FRstructMTC(i).odorOnset(2).BgFg = [allPLV(i).data.baselineFR.odor' allPLV(i).data.FR.odor'];
        %----------------- PLV ---------------------------------------
        %------ analyze normalize PLV -----------
        odorPLVCorrected(i) = allPLV(i).data.normPLV.odor;
        odorLightPLVCorrected(i) = allPLV(i).data.normPLV.odorLight;
        
        %------ analyze biased PLV -----------
        odorBiasedPLV(i) = allPLV(i).data.realPLV.odor;
        odorLightBiasedPLV(i) = allPLV(i).data.realPLV.odorLight;
        
        %------ shuffled PLV -----------
        shuffledDataOdorPLV(i) = (allPLV(i).data.realPLV.odor-allPLV(i).data.normPLV.odor);
        shuffledDataOdorLightPLV(i) = (allPLV(i).data.realPLV.odorLight-allPLV(i).data.normPLV.odorLight);
        
        %---------------------------- SFC ----------------------------
        %------ analyze normalize SFC -----------
        odorSFCCorrected(i,:) = allPLV(i).data.normSFC.odor;
        odorLightSFCCorrected(i,:) = allPLV(i).data.normSFC.odorLight;
        
        %------ analyze biased PLV -----------
        odorBiasedSFC(i,:) = allPLV(i).data.realSFC.odor;
        odorLightBiasedSFC(i,:) = allPLV(i).data.realSFC.odorLight;
        
        %------ shuffled PLV -----------
        shuffledDataOdorSFC(i,:) = (allPLV(i).data.realSFC.odor-allPLV(i).data.normSFC.odor);
        shuffledDataOdorLightSFC(i,:) = (allPLV(i).data.realSFC.odorLight-allPLV(i).data.normSFC.odorLight);
        
        %------------- Phases -------------------
        odorMeanPhase(i) = allPLV(i).data.meanPhase.odor;
        odorLightMeanPhase(i) = allPLV(i).data.meanPhase.odorLight;
        
        % spike phase dist
        allPhaseDistOdor(i,:) = allPLV(i).data.phaseDist.odor;
        allPhaseDistOdorLight(i,:) = allPLV(i).data.phaseDist.odorLight;
        
        %--------------- PSD --------------------------
        spectrumOdor(end+1, :) = 10*log10(allPLV(i).data.PSD.FGspectrum.odor);
        spectrumOdorLight(end+1, :) = 10*log10(allPLV(i).data.PSD.FGspectrum.odorLight);
        
        spectrumOdorRelative(end+1, :) = allPLV(i).data.PSD.FGspectrum.odor./allPLV(1).data.PSD.BGspectrum.odor;
        spectrumOdorLightRelative(end+1, :) = allPLV(i).data.PSD.FGspectrum.odorLight./allPLV(1).data.PSD.BGspectrum.odorLight;
        
        relativeGammaPowerOdor(end+1) = allPLV(i).data.PSD.relativePower.odor;
        relativeGammaPowerOdorLight(end+1) = allPLV(i).data.PSD.relativePower.odorLight;

    end

end
if strcmpi(mouseType,'ArchT')
    % Save the FR structure for CV analysis
    save('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureArchTMTC_for_CV_analysis', 'FRstructMTC');
    
elseif strcmpi(mouseType,'Tbet_NpHR')
    save('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureNpHRMTC_for_CV_analysis', 'FRstructMTC');
end
%% plot

figure;
suptitle(mouseType)
%%%---------------scatter uncorrected and corrected PLV--------------
subplot(3,4,1)
scatter((odorFR'), (odorBiasedPLV'), 'm')
hold on;
scatter((odorFR'), (odorPLVCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('PLV');
legend('Uncorrected PLV', 'Corrected PLV')
title('Odor')

subplot(3,4,2)
scatter((odorLightFR'), (odorLightBiasedPLV'), 'm')
hold on;
scatter((odorLightFR'), (odorLightPLVCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('PLV');
legend('Uncorrected PLV', 'Corrected PLV')
title('Odor + Light')

%%%---------------scatter uncorrected and corrected SFC--------------
subplot(3,4,3)
scatter((odorFR'), (odorBiasedSFC'), 'm')
hold on;
scatter((odorFR'), (odorSFCCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('SFC');
legend('Uncorrected SFC', 'Corrected SFC')
title('Odor')

subplot(3,4,4)
scatter((odorLightFR'), (odorLightBiasedSFC'), 'm')
hold on;
scatter((odorLightFR'), (odorLightSFCCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('SFC');
legend('Uncorrected SFC', 'Corrected SFC')
title('Odor + Light')

hold off;
%%
%%%---------------mean+SEM PLV---------------------
subplot(3,4,5)
meanSEMplv = [mean(shuffledDataOdorPLV(:)) mean(odorBiasedPLV(:)) mean(shuffledDataOdorLightPLV(:)) mean(odorLightBiasedPLV(:)) ;
              std(shuffledDataOdorPLV(:))./sqrt(length(shuffledDataOdorPLV(:))) std(odorBiasedPLV(:))./sqrt(length(odorBiasedPLV(:))) ... 
              std(shuffledDataOdorLightPLV(:))./sqrt(length(shuffledDataOdorLightPLV(:))) std(odorLightBiasedPLV(:))./sqrt(length(odorLightBiasedPLV(:)))];
barwitherr(meanSEMplv(2,:),meanSEMplv(1,:))

[pPLV] = signrank(odorPLVCorrected(:), odorLightPLVCorrected(:));

box off;
ylabel('Mean PLV');
xticklabels({'Shuffled', 'Odor', 'Shuffled' '+Light'});xtickangle(45)
title(['Biased and shuff. Pcorrected = ' num2str(pPLV)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)

%%%---------------mean+SEM SFC---------------------
subplot(3,4,6)
meanSEMsfc = [mean(shuffledDataOdorSFC(:)) mean(odorBiasedSFC(:)) mean(shuffledDataOdorLightSFC(:)) mean(odorLightBiasedSFC(:)) ;
              std(shuffledDataOdorSFC(:))./sqrt(length(shuffledDataOdorSFC(:))) std(odorBiasedSFC(:))./sqrt(length(odorBiasedSFC(:))) ... 
              std(shuffledDataOdorLightSFC(:))./sqrt(length(shuffledDataOdorLightSFC(:))) std(odorLightBiasedSFC(:))./sqrt(length(odorLightBiasedSFC(:)))];
barwitherr(meanSEMsfc(2,:),meanSEMsfc(1,:))

[pSFC] = signrank(odorSFCCorrected(:), odorLightSFCCorrected(:));

box off;
ylabel('Mean SFC');
xticklabels({'Shuffled', 'Odor', 'Shuffled' '+Light'});xtickangle(45)
title(['Biased and shuff. Pcorrected = ' num2str(pSFC)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)

% ------Plot the shuffle corrected paired data with connecting lines--------

%PLV(VS)
subplot(3,4,11)
deltaCorrectedPLVvalues = odorLightPLVCorrected - odorPLVCorrected;%Odor+light - odor
% for random scatter
a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(deltaCorrectedPLVvalues)) + a;

barwitherr(std(deltaCorrectedPLVvalues)./sqrt(length(deltaCorrectedPLVvalues)) ,mean(deltaCorrectedPLVvalues))
hold on;
% plot([odorPLVCorrected' odorLightPLVCorrected']', 'Color', [.5 .5 .5])
scatter(r.*ones(1,length(deltaCorrectedPLVvalues)), deltaCorrectedPLVvalues,10, 'filled', 'MarkerFaceColor', [.5 .5 .5])

box off;
ylabel('\Delta VS (shuffle-corrected');
xticklabels({'VS'});xtickangle(45)
title(['Shuffle corrected values P = ' num2str(pPLV)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)

%SFC
subplot(3,4,12)

deltaCorrectedSFCvalues = odorLightSFCCorrected - odorSFCCorrected;%Odor+light - odor
% for random scatter
a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(deltaCorrectedSFCvalues)) + a;

barwitherr(std(deltaCorrectedSFCvalues)./sqrt(length(deltaCorrectedSFCvalues)) ,mean(deltaCorrectedSFCvalues))
hold on;
% plot([odorPLVCorrected' odorLightPLVCorrected']', 'Color', [.5 .5 .5])
scatter(r.*ones(1,length(deltaCorrectedSFCvalues)), deltaCorrectedSFCvalues,10, 'filled', 'MarkerFaceColor', [.5 .5 .5])

box off;
ylabel('\Delta SFC (shuffle-corrected');
xticklabels({'SFC'});xtickangle(45)
title(['Shuffle corrected values P = ' num2str(pSFC)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
%% SIGNIFICANT RESPONSE + THRESHOLD CROSSING ANALYSIS
%(Only neurons that responded to odor and that after shuft predictor
%correction had a synchronization value > 'minPLvalue'.

% Quantify the precetn change
%Rate
percentChangeSYNCHsfc = (mean(odorLightSFCCorrected)-mean(odorSFCCorrected))./mean(odorSFCCorrected);
percentChangeSYNCHplv = (mean(odorLightPLVCorrected)-mean(odorPLVCorrected))./mean(odorPLVCorrected);
%binomial dist. std
percentChhangeSYNCHsfcSTD = sqrt(abs(percentChangeSYNCHsfc)*(1-abs(percentChangeSYNCHsfc))/length(odorSFCCorrected));
percentChangeSYNCHplvSTD = sqrt(abs(percentChangeSYNCHplv)*(1-abs(percentChangeSYNCHplv))/length(odorPLVCorrected));

subplot(3,4,7)
barwitherr([percentChhangeSYNCHsfcSTD percentChangeSYNCHplvSTD].*100, [percentChangeSYNCHsfc percentChangeSYNCHplv].*100)
box off;
ylabel('% synchronized change +ArchT');
xticklabels({'SFC', 'VS'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('% synch. change')

%------- PLV and SFC for significant odor responses ----
%--- Neurons that responded significantly to odor (P < 0.05, paired t-test)
significtOdorResponseLocs = find(odorResponse < 0.05);
%-PLV
odorResponsePLVodor = odorPLVCorrected(significtOdorResponseLocs);% only P < 0.05
odorResponsePLVodorLight = odorLightPLVCorrected(significtOdorResponseLocs);
%-SFC
odorResponseSFCodor = odorSFCCorrected(significtOdorResponseLocs);
odorResponseSFCodorLight = odorLightSFCCorrected(significtOdorResponseLocs);

%%%%------ threshold for removeing non-locked neurons follwoing shift
%%%%predictor correction (removed ~zeros)
minPLvalue = 0.02;
%%%%%%%%%%%%----------- PLV----------------------
odorPLVCorrected = odorPLVCorrected(find(odorPLVCorrected >= minPLvalue));
odorLightPLVCorrected = odorLightPLVCorrected(find(odorLightPLVCorrected >= minPLvalue));
%%%%%%%%%%%%----------- SFC ----------------------
odorSFCCorrected = odorSFCCorrected(find(odorSFCCorrected >= minPLvalue));
odorLightSFCCorrected = odorLightSFCCorrected(find(odorLightSFCCorrected >= minPLvalue));
 



% Quantify the rate of phase-locked neurons
subplot(3,4,8)
%View neuron-LFP pairs with significant odor response, and phase-locking above threshold.
%%%%%%%%%%%%----------- PLV----------------------
odorResponsePLVodor = odorResponsePLVodor(find(odorResponsePLVodor >= minPLvalue));
odorResponsePLVodorLight = odorResponsePLVodorLight(find(odorResponsePLVodorLight >= minPLvalue));
%%%%%%%%%%%%----------- SFC ----------------------
odorResponseSFCodor = odorResponseSFCodor(find(odorResponseSFCodor >= minPLvalue));
odorResponseSFCodorLight = odorResponseSFCodorLight(find(odorResponseSFCodorLight >= minPLvalue));
%%%%%%%%%%%%----------- PLV----------------------
% rate out of total # neurons
remainingNeuronsOdorTotal = length(odorResponsePLVodor)./length(odorBiasedPLV);
remainingNeuronsOdorLightTotal = length(odorResponsePLVodorLight)./length(odorLightBiasedPLV);

%%%%%%%%%%%%----------- SFC----------------------
% rate out of total # neurons
remainingNeuronsSFCOdorTotal = length(odorResponseSFCodor)./length(odorBiasedSFC);
remainingNeuronsSFCOdorLightTotal = length(odorResponseSFCodorLight)./length(odorLightBiasedSFC);

bar([ remainingNeuronsSFCOdorTotal remainingNeuronsSFCOdorLightTotal remainingNeuronsOdorTotal remainingNeuronsOdorLightTotal]'.*100)
box off;
ylabel('% synchronized  neuron-LFP pairs');
xticklabels({'Odor SFC', '+ChR2 SFC', 'Odor VS', '+ChR2 VS'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('% synch. neurons')
% legend('Out of responding cells', 'Out of all cells')

subplot(3,4,9)
barwitherr([std(odorSFCCorrected)./sqrt(length(odorSFCCorrected)) std(odorLightSFCCorrected)./sqrt(length(odorLightSFCCorrected))], ...
    [mean(odorSFCCorrected) mean(odorLightSFCCorrected)])
box off;
ylabel('SFC corrected');
xticklabels({'Odor', '+ArchT'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold on;
% plot([odorSFCCorrected odorLightSFCCorrected]', 'k')

subplot(3,4,10)
barwitherr([std(odorPLVCorrected)./sqrt(length(odorPLVCorrected)) std(odorLightPLVCorrected)./sqrt(length(odorLightPLVCorrected))], ...
    [mean(odorPLVCorrected) mean(odorLightPLVCorrected)])
box off;
ylabel('VS corrected');
xticklabels({'Odor', '+ArchT'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold on;
% plot([odorPLVCorrected' odorLightPLVCorrected']', 'k')
%%
% if strcmpi(mouseType,'Tbet_NpHR')
% %------------------SYNCH vs FIRING------------------------
% %%%%%% QUANTIFING RATIOS %%%%%%
% % VS ratio corrected
% vsRatioCorrected = odorLightPLVCorrected./odorPLVCorrected;
% % VS ratio biased
% vsRatioBiased = (odorLightBiasedPLV./odorBiasedPLV);
% 
% % SFC ratio corrected
% sfcRatioCorrected = odorLightSFCCorrected./odorSFCCorrected;
% % SFC ratio biased
% sfcRatioBiased = odorLightBiasedSFC./odorBiasedSFC;
% 
% % change of firing rate
% percentChangeFR = odorLightFR./odorFR;
% 
% % bins according to the redusction of MTC
% 
%     percentChangeBins = [0.81 0.61];
%     
%     %%
%     meanSynchLevel.PLVbiased(1,:) =  [100 0];
%     meanSynchLevel.PLVcorrected(1,:) =  [100 0];
%     meanSynchLevel.SFCbiased(1,:) =  [100 0];
%     meanSynchLevel.SFCcorrected(1,:) =  [100 0];
%     firingWindow = 0.2;% the window size above and below the mean value of reduction
%     
%     for i = 1:length(percentChangeBins)
%         %%%% PLV Biased %%%
%         tmpPercent = vsRatioBiased(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%         tmpPercent = tmpPercent(~isinf(tmpPercent));tmpPercent = tmpPercent(~isnan(tmpPercent));
%         meanSynchLevel.PLVbiased(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%         tmpPercent = [];
%         %%%% PLV CORRECTED %%%
%         tmpPercent = vsRatioCorrected(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%         tmpPercent = tmpPercent(~isinf(tmpPercent));tmpPercent = tmpPercent(~isnan(tmpPercent));
%         meanSynchLevel.PLVcorrected(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%         tmpPercent = [];
%         %%%% SFC Biased %%%
%         tmpPercent = sfcRatioBiased(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%         tmpPercent = tmpPercent(~isinf(tmpPercent));tmpPercent = tmpPercent(~isnan(tmpPercent));
%         meanSynchLevel.SFCbiased(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%         tmpPercent = [];
%         %%%% SFC CORRECTED %%%
%         tmpPercent = [];
%         tmpPercent = sfcRatioCorrected(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%         tmpPercent = tmpPercent(~isinf(tmpPercent));tmpPercent = tmpPercent(~isnan(tmpPercent));
%         meanSynchLevel.SFCcorrected(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%         tmpPercent = [];
%         
%     end
%     
%     %%%% Biased %%%
%     subplot(3,4,9)
%     plotBins = [1 percentChangeBins];
%     hold all;
%     errorbar(plotBins, meanSynchLevel.PLVbiased(:,1), meanSynchLevel.PLVbiased(:,2), 'r-o')
%     errorbar(plotBins, meanSynchLevel.SFCbiased(:,1), meanSynchLevel.SFCbiased(:,2), 'b-o')
%     plot([min(plotBins) 1], [100 100], 'k--')
%     set ( gca, 'xdir', 'reverse' )
%     xtickangle(45)
%     ylabel('% Synchronization ratio')
%     xlabel('% Firing rate change')
%     box off;
%     set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
%     set(gca,'linewidth',0.25)
%     title('Biased')
%     legend('VS', 'SFC')
%     %%%% corrected
%     subplot(3,4,10)
%     errorbar(plotBins, meanSynchLevel.PLVcorrected(:,1), meanSynchLevel.PLVcorrected(:,2), 'r-o')
%     hold all;
%     errorbar(plotBins, meanSynchLevel.SFCcorrected(:,1), meanSynchLevel.SFCcorrected(:,2), 'b-o')
%     plot([min(plotBins) 1], [100 100], 'k--')
%     set ( gca, 'xdir', 'reverse' )
%     xtickangle(45)
%     hold off;
%     
%     
%     ylabel('% Synchronization ratio')
%     xlabel('% Firing rate change')
%     box off;
%     set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
%     set(gca,'linewidth',0.25)
%     title('Corrected')
%     legend('VS', 'SFC')
%     
% end
%% ---------------- additional plots -------
%% ---- Phase histogram for significant odor responses ----
% phaseCycleBins = 0:0.1667*pi:2*pi;
phaseCycleBins = 0:0.2*pi:2*pi;%-pi:0.2*pi:pi;
figure;
subplot(1,2,1)
h = histogram(wrapTo2Pi(odorMeanPhase),phaseCycleBins, 'Normalization', 'probability','FaceColor','k','FaceAlpha',.3);
h.DisplayStyle = 'stairs';h.LineWidth = 1.5;h.EdgeColor = 'k';
hold all;
h = histogram(wrapTo2Pi(odorLightMeanPhase),phaseCycleBins, 'Normalization', 'probability','FaceColor','b','FaceAlpha',.3);
h.DisplayStyle = 'stairs';h.LineWidth = 1.5;h.EdgeColor = 'b';
title('Significant phase-locked')
box off;
xlabel('Mean phase preference (radians)');
ylabel('Rate')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
scatter(pi, 0)
hold off;


%% average spike phase dist. (from each neuron-LFP pair, take the normalized spike phase distribution, and average across all significant neuron-LFP pairs)

spikePhaseBins = deg2rad(15:30:360);%-165:30:180;
subplot(1,2,2)
plot(spikePhaseBins , circshift(mean(allPhaseDistOdor),length(spikePhaseBins)/2),'-ks', 'MarkerFaceColor',[.5 .5 .5])
hold all;
plot(spikePhaseBins, circshift(mean(allPhaseDistOdorLight),length(spikePhaseBins)/2),'-bs', 'MarkerFaceColor','b')
title('Spike-dist phase-locked')
% plot the mean phase
spikeDistMeanPhaseOdor = circ_mean(spikePhaseBins', circshift(mean(allPhaseDistOdor),length(spikePhaseBins)/2)');
spikeDistMeanPhaseOdorLight = circ_mean(spikePhaseBins', circshift(mean(allPhaseDistOdorLight),length(spikePhaseBins)/2)');
plot([spikeDistMeanPhaseOdor spikeDistMeanPhaseOdor], [0.087 0.09], 'k')
plot([spikeDistMeanPhaseOdorLight spikeDistMeanPhaseOdorLight], [0.087 0.09], 'b')


xlabel('Spikes phase preference (radians)');
ylabel('Rate')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
xlim([0 2*pi])
hold off;

%% PSD analysis


%Each experiment has its own power spectrum - take only one recording from
%each experiment.
if strcmpi(mouseType,'ArchT')
    psdExperiments = [1 3 5 7 9 11 12 13 15];
elseif strcmpi(mouseType,'Tbet_NpHR')
    psdExperiments = [1:6 8:10 11 13 15];
end
    
uniquePSDstruct = allPLV(psdExperiments);
fPSD = uniquePSDstruct(1).data.PSD.Freq;%Hz

spectrumOdor = [];spectrumOdorLight = [];
spectrumOdorRelative = [];spectrumOdorLightRelative = [];
relativeGammaPowerOdor = [];relativeGammaPowerOdorLight = [];
for psdInd = 1:size(uniquePSDstruct, 2)
    spectrumOdor(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.FGspectrum.odor);
    spectrumOdorLight(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.FGspectrum.odorLight);
    
    spectrumOdorRelative(end+1, :) = uniquePSDstruct(psdInd).data.PSD.FGspectrum.odor./uniquePSDstruct(psdInd).data.PSD.BGspectrum.odor;
    spectrumOdorLightRelative(end+1, :) = uniquePSDstruct(psdInd).data.PSD.FGspectrum.odorLight./uniquePSDstruct(psdInd).data.PSD.BGspectrum.odorLight;
    
    relativeGammaPowerOdor(end+1,1) = uniquePSDstruct(psdInd).data.PSD.relativePower.odor;
    relativeGammaPowerOdorLight(end+1,1) = uniquePSDstruct(psdInd).data.PSD.relativePower.odorLight;

end
%------- PLOT -------------
figure;
suptitle('PSD')
% plot the power in dB
subplot(2,2,1)
errorbar_patch(fPSD, mean(spectrumOdor), std(spectrumOdor)./sqrt(size(spectrumOdor,1)), 'k')
hold all;
errorbar_patch(fPSD, mean(spectrumOdorLight), std(spectrumOdorLight)./sqrt(size(spectrumOdorLight,1)), 'b')
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
title('abs power')
xlim([0 120])
set(gcf, 'Renderer', 'painters')


% plot the power during the trial divided by the bg power
subplot(2,2,2)
errorbar_patch(fPSD, mean(spectrumOdorRelative), std(spectrumOdorRelative)./sqrt(size(spectrumOdorRelative,1)), 'k')
hold all;
errorbar_patch(fPSD, mean(spectrumOdorLightRelative), std(spectrumOdorLightRelative)./sqrt(size(spectrumOdorLightRelative,1)), 'b')
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('relative power')
hold off;
xlim([0 150])

% Bar of the relative power compared with baseline
subplot(2,2,3)
barwitherr(std([relativeGammaPowerOdor relativeGammaPowerOdorLight])./sqrt(length(relativeGammaPowerOdor)), ...
    mean([relativeGammaPowerOdor relativeGammaPowerOdorLight]))
hold all;
plot(1:2, [relativeGammaPowerOdor relativeGammaPowerOdorLight], 'Color', [0.5 .5 .5])
plot(1:2, [1 1], 'k--')
box off;
ylabel('Relative power');
xticklabels({'odor','+Light'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
[~,pPSD] = ttest(relativeGammaPowerOdor, relativeGammaPowerOdorLight);
title(['P = ' num2str(pPSD)])
    
    
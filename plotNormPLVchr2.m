clear;clc
compMethod = 1;%0 = peak detection, 1 = hilbert transform
refractoryPeriodPercentTH = 0.05;
numChannels = 32;

if ismac
    cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data/Awake_data_updated/vsShiftPredictorAnalysis/ChR2')
else
    cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data\Awake_data_updated\vsShiftPredictorAnalysis\ChR2')
end
%% load
%^^^^^^^^^^^^^^^ AWAKE ^^^^^^^^^^^^^^^^^^^
if compMethod==0%0 = peak detection
    load('Awake_40_70Hz_peakDetectionMethod_PLV.mat')
elseif compMethod==1%1 = hilbert transform
%     load('Awake_40_70Hz_PhaseCorrected_SFCcurve_10_120Hz_50msWindow.mat')
      load('Awake_40_70Hz_PhaseCorrected_SFCcurve_10_120Hz_50msWindow_fullFreqs.mat')
end
awake = allPLV;
% %%^^clear neurons awake with >XX% in ref period of 2ms^^
load('allAwakeExps_NEURONEXUS_CHR2.mat')
awakeRefPeriod=cell2mat([allAwakeExps(2:4).refPeriod]);%exps 2-4
awakeBadCells = find(awakeRefPeriod > refractoryPeriodPercentTH);
awake(awakeBadCells) = [];

%^^^^^^^^^^^^^^^ ANESTHETIZED ^^^^^^^^^^^^^^^^^^^
if compMethod==0%0 = peak detection
    load('Anesth_40_70Hz_peakDetectionMethod_PLV.mat')
elseif compMethod==1%1 = hilbert transform
%     load('Anesth_40_70Hz_PhaseCorrected_SFCcurve_10_120Hz_50msWindow.mat')
      load('Anesth_40_70Hz_PhaseCorrected_SFCcurve_10_120Hz_50msWindow_fullFreqs.mat')
end
anesth = allPLV;
% %%%clear neurons anesth with >XX% in ref period of 2ms
load('allAnesthExps_NEURONEXUS_CHR2.mat')
anesthRefPeriod=cell2mat([allAnesthExps.refPeriod]);
anesthBadCells = find(anesthRefPeriod > refractoryPeriodPercentTH);
anesth(anesthBadCells) = [];

clear allPLV;
allPLV = [anesth awake];
stateIndex = [zeros(1,length(anesth)) ones(1,length(awake))];
% allPLV = [awake];
% stateIndex = [ones(1,length(awake))];
% allPLV = [anesth];
% stateIndex = [ones(1,length(anesth))];
%% remove bad neurons and plot stats about the data
badNeuronsState = [];
firingTH = 1;%Hz
minTrialNum = 5;% minimal number of trials with XX spikes/sec.
for ind = 1:size(allPLV, 2)
    
    lowSpikeTrialsOdor = length(find(allPLV(ind).data.FR.odor < firingTH));
    lowSpikeTrialsOdorLight = length(find(allPLV(ind).data.FR.odorLight < firingTH));
    
    if lowSpikeTrialsOdor > minTrialNum | lowSpikeTrialsOdorLight > minTrialNum
        
        fprintf(['Neuron # ' num2str(ind) ' was not included in the analysis \n'])
        badNeuronsState(end+1,:) = [stateIndex(ind) ind];
    end
end
finalAnesthNeurons = length(anesth) - length(find(badNeuronsState(:,1)==0));
finalAwakeNeurons = length(awake) - length(find(badNeuronsState(:,1)==1));

% remove those neurons
allPLV(badNeuronsState(:,2)) = [];
stateIndex(badNeuronsState(:,2)) = [];
%% # of neurons in both anesth. \ awake states.
figure;
subplot(1,2,1)
bar([finalAnesthNeurons finalAwakeNeurons]);
box off;
ylabel('# neurons in analysis');
xticklabels({'Anesth', 'Awake'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title(['Total ' num2str(finalAnesthNeurons+finalAwakeNeurons) ' neurons'])
suptitle('Stats over the data')
%% Init and params
%-PLV
odorPLVCorrected = nan(size(allPLV,2) , size(allPLV(1).data.normPLV.odor,2));
odorLightPLVCorrected = odorPLVCorrected;

odorBiasedPLV= odorPLVCorrected;
odorLightBiasedPLV= odorPLVCorrected;

shuffledDataOdorPLV= odorPLVCorrected;
shuffledDataOdorLightPLV= odorPLVCorrected;

%-phases
odorMeanPhase = odorPLVCorrected;
odorLightMeanPhase = odorPLVCorrected;
odorPcircTest = odorPLVCorrected;
odorLightPcircTest = odorPLVCorrected;

allSpikePhaseOdor = cell(size(allPLV,2) , size(allPLV(1).data.normPLV.odor,2));
allSpikePhaseOdorLight = cell(size(allPLV,2) , size(allPLV(1).data.normPLV.odor,2));

allPhaseDistOdor = allSpikePhaseOdor;
allPhaseDistOdorLight = allSpikePhaseOdorLight;

%-SFC=
odorSFCCorrected= odorPLVCorrected;
odorLightSFCCorrected = odorPLVCorrected;

odorBiasedSFC = odorPLVCorrected;
odorLightBiasedSFC = odorPLVCorrected;
shuffledDataOdorSFC = odorPLVCorrected;
shuffledDataOdorLightSFC = odorPLVCorrected;

%SFC curves
% Frequency vector for SFC curve
if isfield(allPLV(1).data, 'SFCcurvesFrequencyVector')
    SFCcurvesFrequencyVector = allPLV(1).data.SFCcurvesFrequencyVector;
end
SFCcurveOdor = {};SFCcurveOdorLight = {};
SFCcurveOdorSHUFFLED = {};SFCcurveOdorLightSHUFFLED = {};

%-FR (Hz)
odorFR = odorPLVCorrected;odorLightFR = odorPLVCorrected;
odorFRBaseline = odorPLVCorrected;odorLightFRBaseline = odorPLVCorrected;

%-pValue indicating odor response
odorResponse = nan(1,size(allPLV,2));
odorLightResponse = odorResponse;
%FR structure
FRstructMTC = [];
deltaFRMTC = [];


sameElectrodePLVodor = [];
otherElectrodesPLVodor = [];


%% extract the data
for i = 1:size(allPLV, 2)

        
        %------ analyze FR -----------
        odorFR(i,:) = ones(1,length(allPLV(i).data.normPLV.odor)).*mean(allPLV(i).data.FR.odor);
        odorLightFR(i,:) = ones(1,length(allPLV(i).data.normPLV.odorLight)).*mean(allPLV(i).data.FR.odorLight);

        odorFRBaseline(i,:) = ones(1,length(allPLV(i).data.normPLV.odor)).*mean(allPLV(i).data.baselineFR.odor);
        odorLightFRBaseline(i,:) = ones(1,length(allPLV(i).data.normPLV.odorLight)).*mean(allPLV(i).data.baselineFR.odorLight);
        
        odorResponse(i) = allPLV(i).data.odorResponsePvalue.odor;
        odorLightResponse(i) = allPLV(i).data.odorResponsePvalue.odorLight;
        
        %------------- For odor response variability---------------
        FRstructMTC(i).odorOnset(1).BgFg = [allPLV(i).data.baselineFR.odorLight' allPLV(i).data.FR.odorLight'];
        FRstructMTC(i).odorOnset(2).BgFg = [allPLV(i).data.baselineFR.odor' allPLV(i).data.FR.odor'];
        
        %------ DELTA FR ----------
        deltaFRMTC(i) = mean(allPLV(i).data.FR.odorLight) - mean(allPLV(i).data.FR.odor);
        %----------------- PLV ---------------------------------------
        %------ analyze normalize PLV -----------
        odorPLVCorrected(i,:) = allPLV(i).data.normPLV.odor;
        odorLightPLVCorrected(i,:) = allPLV(i).data.normPLV.odorLight;
        
        if ~strcmpi(allPLV(i).recordingElectrodeIndx, 'nan')
            currChannel = str2double(allPLV(i).recordingElectrodeIndx);
            sameElectrodePLVodor(end+1) = allPLV(i).data.normPLV.odor(currChannel);
            otherElectrodesPLVodor(end+1) = mean(setdiff(allPLV(i).data.normPLV.odor, sameElectrodePLVodor));
        end
        
        %------ analyze biased PLV -----------
        odorBiasedPLV(i,:) = allPLV(i).data.realPLV.odor;
        odorLightBiasedPLV(i,:) = allPLV(i).data.realPLV.odorLight;
        
        %------ shuffled PLV -----------
        shuffledDataOdorPLV(i,:) = (allPLV(i).data.realPLV.odor-allPLV(i).data.normPLV.odor);
        shuffledDataOdorLightPLV(i,:) = (allPLV(i).data.realPLV.odorLight-allPLV(i).data.normPLV.odorLight);
        
        %------ phases ---------
        odorMeanPhase(i,:) = allPLV(i).data.meanPhase.odor;
        odorLightMeanPhase(i,:) = allPLV(i).data.meanPhase.odorLight;
        
        
        for spikeInd = 1:size(allPLV(i).data.allSpikePhases.odor,2)
            % Spike dist
            allSpikePhaseOdor{i,spikeInd} = allPLV(i).data.allSpikePhases.odor{spikeInd};
            allSpikePhaseOdorLight{i,spikeInd} = allPLV(i).data.allSpikePhases.odorLight{spikeInd};
            % spike phase dist
            allPhaseDistOdor{i, spikeInd} = allPLV(i).data.phaseDist.odor{spikeInd};
            allPhaseDistOdorLight{i, spikeInd} = allPLV(i).data.phaseDist.odorLight{spikeInd};
        end
        
        %------ Circular test ---------
        odorPcircTest(i,:) = allPLV(i).data.pRayleighTest.odor;
        odorLightPcircTest(i,:) = allPLV(i).data.pRayleighTest.odorLight;
        
        %---------------------------- SFC ----------------------------
        %------ analyze normalize SFC -----------
        odorSFCCorrected(i,:) = allPLV(i).data.normSFC.odor;
        odorLightSFCCorrected(i,:) = allPLV(i).data.normSFC.odorLight;
        
        %------ analyze biased SFC -----------
        odorBiasedSFC(i,:) = allPLV(i).data.realSFC.odor;
        odorLightBiasedSFC(i,:) = allPLV(i).data.realSFC.odorLight;
        
        %------ shuffled SFC -----------
        shuffledDataOdorSFC(i,:) = (allPLV(i).data.realSFC.odor-allPLV(i).data.normSFC.odor);
        shuffledDataOdorLightSFC(i,:) = (allPLV(i).data.realSFC.odorLight-allPLV(i).data.normSFC.odorLight);
        
        %------ SFC curves -----------
%         SFCcurveOdor = [SFCcurveOdor ; allPLV(i).data.realSFCcurve.Odor];
%         SFCcurveOdorLight = [SFCcurveOdorLight ; allPLV(i).data.realSFCcurve.OdorLight];
          SFCcurveOdor{i} = allPLV(i).data.realSFCcurve.Odor';
          SFCcurveOdorLight{i} = allPLV(i).data.realSFCcurve.OdorLight';
          
          %Shuffled curves
          SFCcurveOdorSHUFFLED{i} =  allPLV(i).data.shuffledSFCcurve.Odor';
          SFCcurveOdorLightSHUFFLED{i} =  allPLV(i).data.shuffledSFCcurve.OdorLight';

end
% Save the FR structure for CV analysis
% if ismac
%     save('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureChR2MTC_for_CV_analysis', 'FRstructMTC', 'stateIndex')
% else
%     save('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureChR2MTC_for_CV_analysis', 'FRstructMTC', 'stateIndex');
% end

% Save the deltaFR for MTC firing rate analysis
% save('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\MTCdeltaDRchr2', 'deltaFRMTC');
%% same vs other electrodes
subplot(1,2,2)
bar(mean([sameElectrodePLVodor' otherElectrodesPLVodor']))
hold on;
plot(1:2, [sameElectrodePLVodor' otherElectrodesPLVodor']', 'Color', [0.5 0.5 0.5])
[~,pSameVsOtherElec] = ttest(sameElectrodePLVodor, otherElectrodesPLVodor);

box off;
ylabel('Corrected PLV value');
xticklabels({'Same electrode', 'Other electrodes'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title(['P = ' num2str(pSameVsOtherElec)])
%% plots

figure;
%%%---------------scatter uncorrected and corrected PLV--------------
subplot(3,4,1)
scatter(mean(odorFR'), mean(odorBiasedPLV'), 'm')
hold on;
scatter(mean(odorFR'), mean(odorPLVCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('PLV');
legend('Uncorrected PLV', 'Corrected PLV')
title('Odor')

subplot(3,4,2)
scatter(mean(odorLightFR'), mean(odorLightBiasedPLV'), 'm')
hold on;
scatter(mean(odorLightFR'), mean(odorLightPLVCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('PLV');
legend('Uncorrected PLV', 'Corrected PLV')
title('Odor + ChR2')

%%%---------------scatter uncorrected and corrected SFC--------------
subplot(3,4,3)
scatter(mean(odorFR'), mean(odorBiasedSFC'), 'm')
hold on;
scatter(mean(odorFR'), mean(odorSFCCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('SFC');
legend('Uncorrected SFC', 'Corrected SFC')
title('Odor')

subplot(3,4,4)
scatter(mean(odorLightFR'), mean(odorLightBiasedSFC'), 'm')
hold on;
scatter(mean(odorLightFR'), mean(odorLightSFCCorrected'), 'g')
box off;
xlabel('Firing rate (Hz)');
ylabel('SFC');
legend('Uncorrected SFC', 'Corrected SFC')
title('Odor + ChR2')

hold off;
%% ALL DATA MEAN+-SEM
%%%---------------mean+SEM PLV---------------------
subplot(3,4,5)
meanSEMplv = [mean(shuffledDataOdorPLV(:)) mean(odorBiasedPLV(:)) mean(shuffledDataOdorLightPLV(:)) mean(odorLightBiasedPLV(:)) ;
              std(shuffledDataOdorPLV(:))./sqrt(length(shuffledDataOdorPLV(:))) std(odorBiasedPLV(:))./sqrt(length(odorBiasedPLV(:))) ... 
              std(shuffledDataOdorLightPLV(:))./sqrt(length(shuffledDataOdorLightPLV(:))) std(odorLightBiasedPLV(:))./sqrt(length(odorLightBiasedPLV(:)))];
barwitherr(meanSEMplv(2,:),meanSEMplv(1,:))

[~,pPLV] = ttest(odorPLVCorrected(:), odorLightPLVCorrected(:));

box off;
ylabel('Mean PLV');
xticklabels({'Shuffled', 'Odor', 'Shuffled' '+ChR2'});xtickangle(45)
title(['All data: corrected: P = ' num2str(pPLV)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)

%%%---------------mean+SEM SFC---------------------
subplot(3,4,6)
meanSEMsfc = [mean(shuffledDataOdorSFC(:)) mean(odorBiasedSFC(:)) mean(shuffledDataOdorLightSFC(:)) mean(odorLightBiasedSFC(:)) ;
              std(shuffledDataOdorSFC(:))./sqrt(length(shuffledDataOdorSFC(:))) std(odorBiasedSFC(:))./sqrt(length(odorBiasedSFC(:))) ... 
              std(shuffledDataOdorLightSFC(:))./sqrt(length(shuffledDataOdorLightSFC(:))) std(odorLightBiasedSFC(:))./sqrt(length(odorLightBiasedSFC(:)))];
barwitherr(meanSEMsfc(2,:),meanSEMsfc(1,:))

[~,pSFC] = ttest(odorSFCCorrected(:), odorLightSFCCorrected(:));

box off;
ylabel('Mean SFC');
xticklabels({'Shuffled', 'Odor', 'Shuffled' '+ChR2'});xtickangle(45)
title(['All data: corrected: P = ' num2str(pSFC)])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
%% SIGNIFICANT RESPONSE + THRESHOLD CROSSING ANALYSIS
%(Only neurons that responded to odor and that after shuft predictor
%correction had a synchronization value > 'minPLvalue'.

%------- PLV and SFC for significant odor responses ----
%--- Neurons that responded significantly to odor (P < 0.05, paired t-test)
significtOdorResponseLocs = find(odorResponse < 0.05);
%-PLV
odorResponsePLVodor = odorPLVCorrected(significtOdorResponseLocs,:);% only P < 0.05
odorResponsePLVodorLight = odorLightPLVCorrected(significtOdorResponseLocs,:);% only P < 0.05

%-SFC
odorResponseSFCodor = odorSFCCorrected(significtOdorResponseLocs,:);
odorResponseSFCodorLight = odorLightSFCCorrected(significtOdorResponseLocs,:);

% Assess the rate of synchronized channels in each neuron
% distributionOfSynchLFPChannels(odorResponsePLVodor, odorResponsePLVodorLight, 0.02)

%Analyze the rate of synchronized MTC neuron-LFP pairs as a function of
%the threshold
% rateOfMTCsynchronyChR2_Function_Of_Threshold(odorResponsePLVodor, odorResponsePLVodorLight, odorResponseSFCodor, odorResponseSFCodorLight)


% reshape to a vector
odorResponsePLVodor = odorResponsePLVodor(:);
odorResponsePLVodorLight = odorResponsePLVodorLight(:);

odorResponseSFCodor = odorResponseSFCodor(:);
odorResponseSFCodorLight = odorResponseSFCodorLight(:);




%-SFC curves
odorResponseSFCcurveOdor = SFCcurveOdor(significtOdorResponseLocs);
odorResponseSFCcurveOdor = [odorResponseSFCcurveOdor{:}]';

odorResponseSFCcurveOdorLight = SFCcurveOdorLight(significtOdorResponseLocs);
odorResponseSFCcurveOdorLight = [odorResponseSFCcurveOdorLight{:}]';

% Shuffled
odorResponseSFCcurveOdorSHUFFLED = SFCcurveOdorSHUFFLED(significtOdorResponseLocs);
odorResponseSFCcurveOdorSHUFFLED = [odorResponseSFCcurveOdorSHUFFLED{:}]';

odorResponseSFCcurveOdorLightSHUFFLED = SFCcurveOdorLightSHUFFLED(significtOdorResponseLocs);
odorResponseSFCcurveOdorLightSHUFFLED = [odorResponseSFCcurveOdorLightSHUFFLED{:}]';


%%%%------ threshold for removeing non-locked neurons follwoing shift
%%%%predictor correction (removed ~zeros)
minPLvalue = 0.02;
%%%%%%%%%%%%----------- PLV----------------------
% values to include
testingValuesOdor = find(odorResponsePLVodor >= minPLvalue);
testingValuesOdorLight = find(odorResponsePLVodorLight >= minPLvalue);

% rate of values of the total.
%Rate out of responding to odor
remainingNeuronsOdor = length(testingValuesOdor)./length(odorResponsePLVodor);
remainingNeuronsOdorLight = length(testingValuesOdorLight)./length(odorResponsePLVodorLight);
% rate out of total # neurons
remainingNeuronsOdorTotal = length(testingValuesOdor)./length(odorPLVCorrected(:));
remainingNeuronsOdorLightTotal = length(testingValuesOdorLight)./length(odorLightPLVCorrected(:));
%-Binomial dist. std
% remainingNeuronsOdorSTD = sqrt((remainingNeuronsOdor) * (1-remainingNeuronsOdor) * length(odorResponsePLVodor));
% remainingNeuronsOdorLightSTD = sqrt((remainingNeuronsOdorLight) * (1-remainingNeuronsOdorLight) * length(odorResponsePLVodorLight));

%%%%%%%%%%%%----------- SFC ----------------------
% values to include
testingSFCValuesOdor = find(odorResponseSFCodor >= minPLvalue);
testingSFCValuesOdorLight = find(odorResponseSFCodorLight >= minPLvalue);

% rate of values of the total.
%Rate out of responding to odor
remainingNeuronsSFCOdor = length(testingSFCValuesOdor)./length(odorResponseSFCodor);
remainingNeuronsSFCOdorLight = length(testingSFCValuesOdorLight)./length(odorResponseSFCodorLight);
% rate out of total # neurons
remainingNeuronsSFCOdorTotal = length(testingSFCValuesOdor)./length(odorSFCCorrected(:));
remainingNeuronsSFCOdorLightTotal = length(testingSFCValuesOdorLight)./length(odorLightSFCCorrected(:));
%-Binomial dist. std
% remainingNeuronsSFCOdorSTD = sqrt((remainingNeuronsSFCOdor) * (1-remainingNeuronsSFCOdor) * length(odorResponseSFCodor));
% remainingNeuronsOdorSFCLightSTD = sqrt((remainingNeuronsSFCOdorLight) * (1-remainingNeuronsSFCOdorLight) * length(odorResponseSFCodorLight));

subplot(3,4,7)
% barwitherr([remainingNeuronsOdorSTD remainingNeuronsOdorLightSTD remainingNeuronsSFCOdorSTD remainingNeuronsOdorSFCLightSTD] ... 
%            ,[remainingNeuronsOdor remainingNeuronsOdorLight remainingNeuronsSFCOdor remainingNeuronsSFCOdorLight].*100)
bar([remainingNeuronsOdor remainingNeuronsOdorLight remainingNeuronsSFCOdor remainingNeuronsSFCOdorLight ; ...
     remainingNeuronsOdorTotal remainingNeuronsOdorLightTotal remainingNeuronsSFCOdorTotal remainingNeuronsSFCOdorLightTotal]'.*100)
box off;
ylabel('% synchronized  neuron-LFP pairs');
xticklabels({'Odor VS', '+ChR2 VS', 'Odor SFC', '+ChR2 SFC'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('% synch. neurons')
legend('Out of responding cells', 'Out of all cells')


%% Percent synchronization change in CORRECTED DATA
% Remove neuron-LFP pairs with ~ synch value
%-PLV
odorResponsePLVodor = odorResponsePLVodor(testingValuesOdor);
odorResponsePLVodorLight = odorResponsePLVodorLight(testingValuesOdorLight);

%-SFC
odorResponseSFCodor = odorResponseSFCodor(testingSFCValuesOdor);
odorResponseSFCodorLight = odorResponseSFCodorLight(testingSFCValuesOdorLight);

%-SFC curves
% odorResponseSFCcurveOdor = odorResponseSFCcurveOdor(testingSFCValuesOdor, :);
% odorResponseSFCcurveOdorLight = odorResponseSFCcurveOdorLight(testingSFCValuesOdorLight, :);


% Quantify the precetn change
%Rate
percentChangeSYNCHsfc = (mean(odorResponseSFCodorLight)-mean(odorResponseSFCodor))./mean(odorResponseSFCodor);
percentChangeSYNCHplv = (mean(odorResponsePLVodorLight)-mean(odorResponsePLVodor))./mean(odorResponsePLVodor);
%binomial dist. std
percentChhangeSYNCHsfcSTD = sqrt(percentChangeSYNCHsfc*(1-percentChangeSYNCHsfc)/length(odorResponseSFCodor));
percentChangeSYNCHplvSTD = sqrt(percentChangeSYNCHplv*(1-percentChangeSYNCHplv)/length(odorResponsePLVodor));

subplot(3,4,8)
barwitherr([percentChhangeSYNCHsfcSTD percentChangeSYNCHplvSTD].*100, [percentChangeSYNCHsfc percentChangeSYNCHplv].*100)
box off;
ylabel('% synchronized change +ChR2');
xticklabels({'SFC', 'VS'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('% synch. change')
%% Real value of synch change. CORRECTED
subplot(3,4,9)
%-PLV
[~,pPLVfrCorrected] = ttest2(odorResponsePLVodor, odorResponsePLVodorLight);
%-SFC
[~,pSFCfrCorrected] = ttest2(odorResponseSFCodor, odorResponseSFCodorLight);

barwitherr([std(odorResponseSFCodor)./sqrt(length(odorResponseSFCodor)), std(odorResponseSFCodorLight)./sqrt(length(odorResponseSFCodorLight)), ... %SE
            std(odorResponsePLVodor)./sqrt(length(odorResponsePLVodor)), std(odorResponsePLVodorLight)./sqrt(length(odorResponsePLVodorLight))], ...
            [mean(odorResponseSFCodor) mean(odorResponseSFCodorLight) mean(odorResponsePLVodor) mean(odorResponsePLVodorLight)]);
           
box off;
ylabel({'Synchronization value'; '(Biased - shuffled)'});
xticklabels({'SFC odor','SFC +ChR2', 'VS odor', 'VS +ChR2'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title(['Synch. values, SFC: ' num2str(pSFCfrCorrected) ', VS: ' num2str(pPLVfrCorrected)])


%% BIASED VALUES
%------- PLV and SFC for significant odor responses (NOT CORRECTED) ----
%--- Neurons that responded significantly to odor (P < 0.05, paired t-test)
%-PLV
odorResponsePLVodorBIASED = odorBiasedPLV(significtOdorResponseLocs,:);odorResponsePLVodorBIASED = odorResponsePLVodorBIASED(:);% only P < 0.05
odorResponsePLVodorLightBIASED = odorLightBiasedPLV(significtOdorResponseLocs,:);odorResponsePLVodorLightBIASED = odorResponsePLVodorLightBIASED(:);
%-SFC
odorResponseSFCodorBIASED = odorBiasedSFC(significtOdorResponseLocs,:);odorResponseSFCodorBIASED = odorResponseSFCodorBIASED(:);
odorResponseSFCodorLightBIASED = odorLightBiasedSFC(significtOdorResponseLocs,:);odorResponseSFCodorLightBIASED = odorResponseSFCodorLightBIASED(:);


% ------Remove neuron-LFP pairs with ~ synch value ------
%-PLV
odorResponsePLVodorBIASED = odorResponsePLVodorBIASED(testingValuesOdor);
odorResponsePLVodorLightBIASED = odorResponsePLVodorLightBIASED(testingValuesOdorLight);

%-SFC
odorResponseSFCodorBIASED = odorResponseSFCodorBIASED(testingSFCValuesOdor);
odorResponseSFCodorLightBIASED = odorResponseSFCodorLightBIASED(testingSFCValuesOdorLight);


%%%%% STATISTICS OVER THE CORRECTED DATA %%%%
% %-PLV
% [~,pPLVfrCorrected] = ttest2(odorResponsePLVodor, odorResponsePLVodorLight);
% %-SFC
% [~,pSFCfrCorrected] = ttest2(odorResponseSFCodor, odorResponseSFCodorLight);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,10)
%-Plot the biased SFC mean+-SEM
barwitherr([std(odorResponseSFCodorBIASED)./sqrt(length(odorResponseSFCodorBIASED)), std(odorResponseSFCodorLightBIASED)./sqrt(length(odorResponseSFCodorLightBIASED))], ...
            [mean(odorResponseSFCodorBIASED) mean(odorResponseSFCodorLightBIASED)]);
           
box off;
ylabel('Spike-field coherence');
xticklabels({'odor','+ChR2'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title(['Spike-field coherence (biased)'])

subplot(3,4,11)
%-Plot the biased PLV mean+-SEM
barwitherr([std(odorResponsePLVodorBIASED)./sqrt(length(odorResponsePLVodorBIASED)), std(odorResponsePLVodorLightBIASED)./sqrt(length(odorResponsePLVodorLightBIASED))], ...
            [mean(odorResponsePLVodorBIASED) mean(odorResponsePLVodorLightBIASED)]);
           
box off;
ylabel('Vector strength');
xticklabels({'odor','+ChR2'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title(['VS (biased)'])
%% ------------------SYNCH vs FIRING------------------------
%%%%%% QUANTIFING RATIOS %%%%%%
% VS ratio corrected
vsRatioCorrected = (mean(odorLightPLVCorrected')./mean(odorPLVCorrected'));vsRatioCorrected = vsRatioCorrected(:);
% VS ratio biased
vsRatioBiased = (mean(odorLightBiasedPLV')./mean(odorBiasedPLV'));vsRatioBiased = vsRatioBiased(:);

% SFC ratio corrected
sfcRatioCorrected = (mean(odorLightSFCCorrected')./mean(odorSFCCorrected'));sfcRatioCorrected = sfcRatioCorrected(:);
% SFC ratio biased
sfcRatioBiased = (mean(odorLightBiasedSFC')./mean(odorBiasedSFC'));sfcRatioBiased = sfcRatioBiased(:);

% change of firing rate
percentChangeFR = (mean(odorLightFR')./mean(odorFR'));percentChangeFR = percentChangeFR(:);
%%
percentChangeBins = [1 0.85 0.7 0.55];

meanSynchLevel.PLVbiased(1,:) =  [100 0];
meanSynchLevel.PLVcorrected(1,:) =  [100 0];
meanSynchLevel.SFCbiased(1,:) =  [100 0];
meanSynchLevel.SFCcorrected(1,:) =  [100 0];
meanSynchLevel.PLVcorrectedN = {};
meanSynchLevel.SFCcorrectedN = {};
firingWindow = 0.05;% the window size above and below the mean value of reduction

for i = 1:length(percentChangeBins)-1
      tmpPercent = [];
%     %%%% PLV Biased %%%
%     tmpPercent = vsRatioBiased(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%     tmpPercent = tmpPercent(~isinf(tmpPercent));
%     meanSynchLevel.PLVbiased(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%     tmpPercent = [];


    %%%% PLV CORRECTED %%%
    tmpPercent = vsRatioCorrected(find(percentChangeFR < percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%     tmpPercent = vsRatioCorrected(find(percentChangeFR < percentChangeBins(i) & percentChangeFR >= percentChangeBins(i+1)));
    tmpPercent = tmpPercent(~isinf(tmpPercent));
    meanSynchLevel.PLVcorrected(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
    meanSynchLevel.PLVcorrectedN{end+1} = num2str(length(tmpPercent));
    tmpPercent = [];
    
    
%     %%%% SFC Biased %%%
%     tmpPercent = sfcRatioBiased(find(percentChangeFR <= percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%     tmpPercent = tmpPercent(~isinf(tmpPercent));
%     meanSynchLevel.SFCbiased(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
%     tmpPercent = [];


    %%%% SFC CORRECTED %%%
    tmpPercent = [];
    tmpPercent = sfcRatioCorrected(find(percentChangeFR < percentChangeBins(i)+firingWindow & percentChangeFR >= percentChangeBins(i)-firingWindow));
%     tmpPercent = sfcRatioCorrected(find(percentChangeFR < percentChangeBins(i) & percentChangeFR >= percentChangeBins(i+1)));
    tmpPercent = tmpPercent(~isinf(tmpPercent));
    meanSynchLevel.SFCcorrected(end+1,:) = [mean(tmpPercent) std(tmpPercent)./sqrt(length(tmpPercent))].*100;
    meanSynchLevel.SFCcorrectedN{end+1} = num2str(length(tmpPercent));
    tmpPercent = [];

end

% %%%% Biased %%%
% subplot(3,4,11)
% plotBins = [1 percentChangeBins];
% hold all;
% errorbar(plotBins, meanSynchLevel.PLVbiased(:,1), meanSynchLevel.PLVbiased(:,2), 'r-o')
% errorbar(plotBins, meanSynchLevel.SFCbiased(:,1), meanSynchLevel.SFCbiased(:,2), 'b-o')
% plot([min(plotBins) 1], [100 100], 'k--')
% set ( gca, 'xdir', 'reverse' )
% xtickangle(45)
% ylabel('% Synchronization ratio')
% xlabel('% Firing rate change')
% box off;
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% set(gca,'linewidth',0.25)
% title('Biased')
% legend('VS', 'SFC')

%%%% corrected
subplot(3,4,12)
hold all;
% plotBins = [1 percentChangeBins];
plotBins = [percentChangeBins];
errorbar(plotBins, meanSynchLevel.PLVcorrected(:,1), meanSynchLevel.PLVcorrected(:,2), 'b-s', 'MarkerFaceColor', 'b', 'LineWidth', 0.5)
errorbar(plotBins, meanSynchLevel.SFCcorrected(:,1), meanSynchLevel.SFCcorrected(:,2), 'r-s', 'MarkerFaceColor', 'r', 'LineWidth', 0.5)
plot([min(plotBins) 1], [100 100], 'k--')
set ( gca, 'xdir', 'reverse' )
xtickangle(45)
% text(percentChangeBins(2:end), ones(1,length(percentChangeBins(2:end))).*60, meanSynchLevel.PLVcorrectedN)

ylabel('% Synchronization ratio')
xlabel('% Firing rate change')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
title('Corrected')
legend('VS', 'SFC')
hold off

%% Data Points and Histograms of corrected phase locking data
figure;
subplot(2,2,1)
scatter(ones(1, length(odorResponsePLVodor)), odorResponsePLVodor, 'k')
hold on;
scatter(2.*ones(1, length(odorResponsePLVodorLight)), odorResponsePLVodorLight, 'b')
xlim([0 3])
%Plot the distibution of values for the VS analysis
subplot(2,2,2)
histogram(odorResponsePLVodor,0:0.02:0.5, 'Normalization', 'probability', 'EdgeColor', 'k', 'DisplayStyle', 'stairs')
hold on;
histogram(odorResponsePLVodorLight,0:0.02:0.5, 'Normalization', 'probability', 'EdgeColor', 'b', 'DisplayStyle', 'stairs')
set(gca,'view',[90 -90])
ylim([0 0.35])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
box off;
[~,pKScorrectedVS] = kstest2(odorResponsePLVodor, odorResponsePLVodorLight);
title(['PLV, P = ' num2str(pKScorrectedVS)])
legend(['Odor, N = ' num2str(length(odorResponsePLVodor))], ['+ChR2, N = ' num2str(length(odorResponsePLVodor))])
set(gca, 'YColor', 'k')
set(gca, 'XColor', 'k')

subplot(2,2,2)
histogram(odorResponseSFCodor,0:0.04:0.4, 'Normalization', 'probability', 'EdgeColor', 'k', 'DisplayStyle', 'stairs')
hold on;
histogram(odorResponseSFCodorLight,0:0.04:0.4, 'Normalization', 'probability', 'EdgeColor', 'b', 'DisplayStyle', 'stairs')
ylim([0 0.6])
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
box off;
[~,pKScorrectedSFC] = kstest2(odorResponseSFCodor, odorResponseSFCodorLight);
title(['SFC, P = ' num2str(pKScorrectedSFC)])
legend(['Odor, N = ' num2str(length(odorResponseSFCodor))], ['+ChR2, N = ' num2str(length(odorResponseSFCodorLight))])
set(gca, 'YColor', 'k')
set(gca, 'XColor', 'k')


figure;
subplot(1,2,1)
ecdf(odorResponsePLVodor)
hold on;
ecdf(odorResponsePLVodorLight)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
box off;
set(gca, 'YColor', 'k')
set(gca, 'XColor', 'k')
ylim([0 1.02])
legend('Odor', '+ChR2')
title('VS')

subplot(1,2,2)
ecdf(odorResponseSFCodor)
hold on;
ecdf(odorResponseSFCodorLight)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
box off;
set(gca, 'YColor', 'k')
set(gca, 'XColor', 'k')
ylim([0 1.02])
legend('Odor', '+ChR2')
title('SFC')
%% ---------------- additional plots -------
%% ---- Phase histogram for significant odor responses ----
odorResponsePhaseodor = odorMeanPhase(find(odorResponse < 0.05),:);odorResponsePhaseodor = odorResponsePhaseodor(:);
odorResponsePhaseodorLight = odorLightMeanPhase(find(odorResponse < 0.05),:);odorResponsePhaseodorLight = odorResponsePhaseodorLight(:);

odorResponsePhaseodor = odorResponsePhaseodor(testingValuesOdor);
odorResponsePhaseodorLight = odorResponsePhaseodorLight(testingValuesOdorLight);

% phaseCycleBins = 0:0.1667*pi:2*pi;
phaseCycleBins = 0:0.2*pi:2*pi;%-pi:0.2*pi:pi;
figure;
subplot(1,2,1)
h = histogram(wrapTo2Pi(odorResponsePhaseodor),phaseCycleBins, 'Normalization', 'probability','FaceColor','k','FaceAlpha',.3);
h.DisplayStyle = 'stairs';h.LineWidth = 1.5;h.EdgeColor = 'k';
hold all;
h = histogram(wrapTo2Pi(odorResponsePhaseodorLight),phaseCycleBins, 'Normalization', 'probability','FaceColor','b','FaceAlpha',.3);
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
allPhaseDistOdor = allPhaseDistOdor(significtOdorResponseLocs,:);allPhaseDistOdor = allPhaseDistOdor(:);
allPhaseDistOdorLight = allPhaseDistOdorLight(significtOdorResponseLocs,:);allPhaseDistOdorLight = allPhaseDistOdorLight(:);

allPhaseDistOdor = allPhaseDistOdor(testingValuesOdor);
allPhaseDistOdorLight = allPhaseDistOdorLight(testingValuesOdorLight);
spikePhaseBins = deg2rad(15:30:360);%-165:30:180;
subplot(1,2,2)
plot(spikePhaseBins , circshift(mean(cell2mat(allPhaseDistOdor)),length(spikePhaseBins)/2),'-ks', 'MarkerFaceColor',[.5 .5 .5])
hold all;
plot(spikePhaseBins, circshift(mean(cell2mat(allPhaseDistOdorLight)),length(spikePhaseBins)/2),'-bs', 'MarkerFaceColor','b')
title('Spike-dist phase-locked')
% plot the mean phase
spikeDistMeanPhaseOdor = circ_mean(spikePhaseBins', circshift(mean(cell2mat(allPhaseDistOdor)),length(spikePhaseBins)/2)');
spikeDistMeanPhaseOdorLight = circ_mean(spikePhaseBins', circshift(mean(cell2mat(allPhaseDistOdorLight)),length(spikePhaseBins)/2)');
plot([spikeDistMeanPhaseOdor spikeDistMeanPhaseOdor], [0.087 0.09], 'k')
plot([spikeDistMeanPhaseOdorLight spikeDistMeanPhaseOdorLight], [0.087 0.09], 'b')


xlabel('Mean phase preference (radians)');
ylabel('Rate')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
xlim([0 2*pi])
hold off;

suptitle('Spike-phase analyses')
%% PSD analysis
%Each experiment has its own power spectrum - take only one recording from
%each experiment.
psdExperiments = [1 9 14 18 22 28 31 36 38 54];
uniquePSDstruct = allPLV(psdExperiments);
fPSD = 0:.5:2000;%Hz

spectrumOdor = [];spectrumOdorLight = [];
spectrumOdorRelative = [];spectrumOdorLightRelative = [];
relativeGammaPowerOdor = [];relativeGammaPowerOdorLight = [];
bgSpectrumOdor = [];bgSpectrumOdorLight = [];
for psdInd = 1:size(uniquePSDstruct, 2)
    spectrumOdor(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.FGspectrum.odor);
    spectrumOdorLight(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.FGspectrum.odorLight);
    
    bgSpectrumOdor(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.BGspectrum.odor);
    bgSpectrumOdorLight(end+1, :) = 10*log10(uniquePSDstruct(psdInd).data.PSD.BGspectrum.odorLight);
    
    spectrumOdorRelative(end+1, :) = uniquePSDstruct(psdInd).data.PSD.FGspectrum.odor./uniquePSDstruct(psdInd).data.PSD.BGspectrum.odor;
    spectrumOdorLightRelative(end+1, :) = uniquePSDstruct(psdInd).data.PSD.FGspectrum.odorLight./uniquePSDstruct(psdInd).data.PSD.BGspectrum.odorLight;
    
    relativeGammaPowerOdor(end+1,1) = uniquePSDstruct(psdInd).data.PSD.relativePower.odor;
    relativeGammaPowerOdorLight(end+1,1) = uniquePSDstruct(psdInd).data.PSD.relativePower.odorLight;

end
%------- PLOT -------------
%%%% PLOT EXAMPLES
examplesBarLines = [];
figure;
subplot(1,3,1)
examplesTrials = 7;
plot(fPSD, spectrumOdor(examplesTrials,:), 'k')
hold all;
plot(fPSD, mean([bgSpectrumOdor(examplesTrials,:) ; bgSpectrumOdorLight(examplesTrials,:)]), '--', 'Color', [0.5 .5 .5])
xlim([0 120])
plot(fPSD, spectrumOdorLight(examplesTrials,:), 'b')
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
examplesBarLines = [examplesBarLines examplesTrials];

%examples #2
subplot(1,3,2)
examplesTrials = 6;
plot(fPSD, spectrumOdor(examplesTrials,:), 'k')
hold all;
plot(fPSD, mean([bgSpectrumOdor(examplesTrials,:) ; bgSpectrumOdorLight(examplesTrials,:)]), '--', 'Color', [0.5 .5 .5])
xlim([0 120])
plot(fPSD, spectrumOdorLight(examplesTrials,:), 'b')
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
examplesBarLines = [examplesBarLines examplesTrials];

%examples #2
subplot(1,3,3)
examplesTrials = 10;%Awake
plot(fPSD, spectrumOdor(examplesTrials,:), 'k')
hold all;
plot(fPSD, mean([bgSpectrumOdor(examplesTrials,:) ; bgSpectrumOdorLight(examplesTrials,:)]), '--', 'Color', [0.5 .5 .5])
xlim([0 120])
plot(fPSD, spectrumOdorLight(examplesTrials,:), 'b')
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
examplesBarLines = [examplesBarLines examplesTrials];


% PLOT MEAN
figure;
suptitle('PSD')
% plot the power in dB
subplot(2,3,1)
errorbar_patch(fPSD, mean(spectrumOdor), std(spectrumOdor)./sqrt(size(spectrumOdor,1)), 'k')
hold all;
errorbar_patch(fPSD, mean(spectrumOdorLight), std(spectrumOdorLight)./sqrt(size(spectrumOdorLight,1)), 'b')

plot(fPSD, mean([bgSpectrumOdorLight ; bgSpectrumOdor]), 'r--');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
title('abs power')
xlim([0 100])
set(gcf, 'Renderer', 'painters')
% legend('+ChR2', 'Odor', 'BG')


% plot the power during the trial divided by the bg power
subplot(2,3,2)
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
subplot(2,3,3)
barwitherr(std([relativeGammaPowerOdor relativeGammaPowerOdorLight])./sqrt(length(relativeGammaPowerOdor)), ...
    mean([relativeGammaPowerOdor relativeGammaPowerOdorLight]))
hold all;
plot(1:2, [relativeGammaPowerOdor relativeGammaPowerOdorLight], 'Color', [0.5 .5 .5])
plot(1:2, [relativeGammaPowerOdor(examplesBarLines) relativeGammaPowerOdorLight(examplesBarLines)], 'Color', 'g')

plot(1:2, [1 1], 'k--')
box off;
ylabel('Relative power');
xticklabels({'odor','+ChR2'});xtickangle(45)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
[~,pPSD] = ttest(relativeGammaPowerOdor, relativeGammaPowerOdorLight);
title(['P = ' num2str(pPSD)])

%Plot the SFC curves
subplot(2,3,4)
plot(SFCcurvesFrequencyVector, smooth(mean(odorResponseSFCcurveOdor), 5), 'k', 'LineWidth', 1.5);
hold all;
plot(SFCcurvesFrequencyVector,smooth(mean(odorResponseSFCcurveOdorLight), 5), 'b', 'LineWidth', 1.5);

box off;
ylabel('Spike-field coherence');
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
xlabel('Frequency (Hz)')
hold off;
xlim([10 100])


subplot(2,3,5)
errorbar_patch(SFCcurvesFrequencyVector, smooth(mean(odorResponseSFCcurveOdorLight-odorResponseSFCcurveOdor),5), std(odorResponseSFCcurveOdorLight-odorResponseSFCcurveOdor)./sqrt(size(odorResponseSFCcurveOdorLight-odorResponseSFCcurveOdor,1)))
xlim([20 100])
box off;
ylabel('\Delta Spike-field coherence');
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
xlabel('Frequency (Hz)')
hold on;
plot([20 100], [0 0], 'k--')
hold off;
set(gcf,'Renderer', 'painter');

suptitle('PSD and SFC curves')
    
%%%%% Quantify SFC in other frequencies------
otherFrequencyRange = [10 20];
SFCcurvesFrequencyVector = ceil(SFCcurvesFrequencyVector);
betaRangeIndices = find(SFCcurvesFrequencyVector >= otherFrequencyRange(1) & SFCcurvesFrequencyVector <= otherFrequencyRange(2));
%----ODOR
% SFC odor real
betaValuesOdor = mean(odorResponseSFCcurveOdor(:,betaRangeIndices)');

% SFC odor shuffled
betaValuesOdorSHUFFLED = mean(odorResponseSFCcurveOdorSHUFFLED(:,betaRangeIndices)');

%----ODOR+ChR2
% SFC real
betaValuesOdorLight = mean(odorResponseSFCcurveOdorLight(:,betaRangeIndices)');

% SFC  shuffled
betaValuesOdorLightSHUFFLED = mean(odorResponseSFCcurveOdorLightSHUFFLED(:,betaRangeIndices)');
%---------delta
%----ODOR
deltaBetaOdor = (betaValuesOdor - betaValuesOdorSHUFFLED);
deltaBetaOdor = deltaBetaOdor(find(deltaBetaOdor > minPLvalue));
%-----ODOR+ChR2
deltaBetaOdorLight = (betaValuesOdorLight - betaValuesOdorLightSHUFFLED);
deltaBetaOdorLight = deltaBetaOdorLight(find(deltaBetaOdorLight > minPLvalue));

[~,pBeta] = ttest2(deltaBetaOdorLight, deltaBetaOdor);
display(['N_odor: ' num2str(length(deltaBetaOdor)) '; N_odor+light: ' num2str(length(deltaBetaOdorLight)) '; P = ' num2str(pBeta)])

% Quantify the precent change
%Rate
percentChangeSYNCHbeta = (mean(deltaBetaOdorLight)-mean(deltaBetaOdor))./mean(deltaBetaOdor);
%binomial dist. std
percentChangeSYNCHbetaSTD = sqrt(percentChangeSYNCHbeta*(1-percentChangeSYNCHbeta)/length(deltaBetaOdor));


subplot(2,3,5)
barwitherr(percentChangeSYNCHbetaSTD.*100, percentChangeSYNCHbeta .*100)
ylabel('%');xlabel(num2str(otherFrequencyRange))
    
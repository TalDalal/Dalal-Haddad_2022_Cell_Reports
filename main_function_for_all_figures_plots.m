%% FIGURE 1-MAIN
%Data example
plotWaveletTrialExample

%MTC firing rates
scatterMTCfiringBeforeAfterLight('chr2')

% Gamma synchrony
% compute for a single neuron
shiftPredictorSpikeLFPsynch(matName, 'hilbert')
% run over all neurons
runShiftPredictorSynchAllNeurons(mouseState, 'hilbert', savingPath)
%plot the results
plotNormPLVchr2

%% SUPPLEMENTAL FIGURE 1

%MTC firing rates
scatterMTCfiringBeforeAfterLight('arch')

%Synchrony
shiftPredictorSpikeLFPsynchDAQexps(matName, methodFlag);
runShiftPredictorSynchAllNeuronsDAQexps('ArchT', 'hilbert', []);
% plot phase-locking values
plotNormPLV_DAQexps

%% FIGURE 2-MAIN
%CHR2
%%% PC ChR2 scatter plot, average PSTH, percent change, awake aneth
%%% comparison
scatterAllDataPCchr2(2, 'chr2', 'single');
%%% ChR2 light inetnsity - PC and MTC
chr2LightIntensityAnalysis

%ArchT
%%% PC ArchT scatter plot and percent change
scatterAllDataPCchr2(2, 'arch', 'single');
%%% ChR2 light inetnsity - PC and MTC
ArchTLightIntensityAnalysis

%% SUPPLEMENTAL FIGURE 2
%All plots are included in the main figure functions.

%% FIGURE 3-MAIN
%PCA analysis
pcaSpheresPCneurons

%%% Correlation - including lifetime and population sparseness analysis
corrStimDuration

%%% Decoding - delta
decodingDeltaSuccRate;

%%% CV - PC neurons
cvBothConditionsPercent;

% CV of MTC
cvMTC

%%% Within correlation
compareBetweenWithinCorrelation

%% SUPPLEMENTAL FIGURE 3
% chnage in # of effective odors
computeTuningCurvePCchr2('chr2')
computeTuningCurvePCchr2('arch')


%% FIGURE 4-MAIN

%%%MTC - firing rate
scatterMTCfiringBeforeAfterLight('tbet-halo')

%Synchrony
shiftPredictorSpikeLFPsynchDAQexps(matName, methodFlag);
runShiftPredictorSynchAllNeuronsDAQexps('Tbet_ArchT', 'hilbert', []);
% All phase-locking values
plotNormPLV_DAQexps

%%% PC - scatter firing rate
scatterAllDataPCchr2(2,'Tbet-NpHR', 'single');

%%% Tbet-NpHR light intensity
TbetNpHRLightIntensityAnalysis


%% SUPPLEMENTAL FIGURE 4
%%% PC - scatter firing rate
scatterAllDataPCchr2(2,'JGC', 'single');
%JGC light intensity analysis
JGchr2LightIntensityAnalysis

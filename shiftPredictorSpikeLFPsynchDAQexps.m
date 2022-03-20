function [outputPLV] = shiftPredictorSpikeLFPsynchDAQexps(matName, methodFlag)
%% load
load(matName)
computePSD = 1;
%% Params
params.bgTimeBounds = [-2 0];
params.fgTimeBounds = [0 2];
params.decimateFs = expLFP_SPIKES(1).decimateFs;
params.lfpRange = [40 70];
params.stimType = expLFP_SPIKES(1).stimuliTyps;%0=odor, 1 = odor+light, 2 = light.
params.numSpikesThreshold = 5;
params.t = expLFP_SPIKES(1).t;% time
params.winSize = 0.05;%50ms
params.phaseBins = deg2rad(-165:30:180);
params.SFCcurveFlag = 0;


%% %%%% FILTERS %%%%
lpFiltLFP = designfilt('lowpassiir','FilterOrder',10, ...
    'PassbandFrequency',params.lfpRange(2),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

hpFiltLFP = designfilt('highpassiir','FilterOrder',8, ...
    'PassbandFrequency',params.lfpRange(1),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

odorTrials = find(params.stimType==0);
odorLightTrials = find(params.stimType==1);

%%

for trialsLoc = 1:length(odorTrials)
    %-------------- spikes per condition ---------
    odorSpikes{trialsLoc} = expLFP_SPIKES(odorTrials(trialsLoc)).spikes;
    odorLightSpikes{trialsLoc} = expLFP_SPIKES(odorLightTrials(trialsLoc)).spikes;
    
    %%%-------------- divide LFP by conditions ------------------
    odorLFP(trialsLoc,:) = expLFP_SPIKES(odorTrials(trialsLoc)).LFP';
    odorlightLFP(trialsLoc,:) = expLFP_SPIKES(odorLightTrials(trialsLoc)).LFP';
end

%%%%% FLIP THE LFP BACK (IT WAS FLIPED IN THE MATLAB ANALYSIS)
odorLFP = odorLFP.*-1;
odorlightLFP = odorlightLFP.*-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init values
firingRateOdor = nan(1,length(odorSpikes));
firingRateOdorLight = firingRateOdor;
bgFiringOdor = firingRateOdor;
bgFiringOdorLight = firingRateOdor;

%%
%%%-------------- filter the LFPs------------------------
%low pass filter
odorLightLFPbpf = filtfilt(lpFiltLFP, odorlightLFP');
odorLFPbpf = filtfilt(lpFiltLFP, odorLFP');

% high pass filter
odorLightLFPbpf = filtfilt(hpFiltLFP, odorLightLFPbpf);
odorLFPbpf = filtfilt(hpFiltLFP, odorLFPbpf);


%%%------------- extract the trial duration signals------------------
FGodorLightLFPbpf = odorLightLFPbpf(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);
FGodorLFPbpf = odorLFPbpf(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);

%%%------------- Shifted signals------------------
FGodorLightLFPbpfSHIFTED = circshift(FGodorLightLFPbpf,1,2);% shift of 1 trial
FGodorLFPbpfSHIFTED = circshift(FGodorLFPbpf, 1, 2);% shift of 1 trial

%%%% init
allShiftedOdorSpikePhases = [];
allRealOdorSpikePhases = [];
allShiftedOdorLightSpikePhases = [];
allRealOdorLightSpikePhases = [];

realOdorSFC = nan(1,size(FGodorLFPbpf,2));
realOdorLightSFC = realOdorSFC;
shiftedOdorSFC = realOdorSFC;
shiftedOdorLightSFC = realOdorSFC;

%------------------------ PSD ------------------------------------
if computePSD == 1
    %- trial duration UN-FILTERED
    unFiteredLFP.FGunFilteredLFPodor = odorLFP(:,find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)));
    unFiteredLFP.FGunFilteredLFPodorLight = odorlightLFP(:,find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)));
    
    %- baseline UN-FILTERED
    unFiteredLFP.BGunFilteredLFPodor = odorLFP(:,find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)));
    unFiteredLFP.BGunFilteredLFPodorLight = odorlightLFP(:,find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)));
    
    psdStructure = multiTaperPSD(unFiteredLFP, params); 
end



for j = 1:size(FGodorLFPbpf,2)%%%%% run over trials %%%%

    %% curr trial data
    
    spikeLFPdata.spikes.odor = odorSpikes{j};
    spikeLFPdata.spikes.odorLight = odorLightSpikes{j};
    
    % baseline firing rate (Hz)
    bgFiringOdor(j) = length(spikeLFPdata.spikes.odor(find(spikeLFPdata.spikes.odor >= params.bgTimeBounds(1) ...
        & spikeLFPdata.spikes.odor < params.bgTimeBounds(2))))./diff(params.fgTimeBounds);
    
    bgFiringOdorLight(j) = length(spikeLFPdata.spikes.odorLight(find(spikeLFPdata.spikes.odorLight >= params.bgTimeBounds(1) ...
        & spikeLFPdata.spikes.odorLight < params.bgTimeBounds(2))))./diff(params.fgTimeBounds);
    
    %%%------------- extract the trial duration of spikes------------------
    spikeLFPdata.spikes.odor = spikeLFPdata.spikes.odor(find(spikeLFPdata.spikes.odor >= params.fgTimeBounds(1) & spikeLFPdata.spikes.odor < params.fgTimeBounds(2)));
    spikeLFPdata.spikes.odorLight = spikeLFPdata.spikes.odorLight(find(spikeLFPdata.spikes.odorLight >= params.fgTimeBounds(1) & spikeLFPdata.spikes.odorLight < params.fgTimeBounds(2)));
    
    spikeLFPdata.realLFP.odor = FGodorLFPbpf(:,j);
    spikeLFPdata.realLFP.odorLight = FGodorLightLFPbpf(:,j);
    
    spikeLFPdata.shiftedLFP.odor = FGodorLFPbpfSHIFTED(:,j);
    spikeLFPdata.shiftedLFP.odorLight = FGodorLightLFPbpfSHIFTED(:,j);
    
    % FG firing rate (Hz)
    firingRateOdor(j) = length(spikeLFPdata.spikes.odor)./diff(params.fgTimeBounds);
    firingRateOdorLight(j) = length(spikeLFPdata.spikes.odorLight)./diff(params.fgTimeBounds);

    %% compute the VS (and shift predictor)
    
    if strcmpi(methodFlag, 'peakDetection')
        [spikePhase] = computeVSpeakDetection(spikeLFPdata, params);
        
    elseif strcmpi(methodFlag, 'hilbert')
        [spikePhase] = computeVShilbertTransform(spikeLFPdata, params);
    end
    
    %%%------------ assign the output phases----------------------
    allRealOdorSpikePhases = [allRealOdorSpikePhases spikePhase.realOdor];
    
    allShiftedOdorSpikePhases = [allShiftedOdorSpikePhases spikePhase.shiftedOdor];
    
    allRealOdorLightSpikePhases = [allRealOdorLightSpikePhases spikePhase.realOdorLight];
    
    allShiftedOdorLightSpikePhases = [allShiftedOdorLightSpikePhases spikePhase.shiftedOdorLight];
    
    %% Compute the SFC (and Shift predictor).
    
    [SFCdata] = computeSpikeFieldCoherenceShiftPredictor(spikeLFPdata, params);
    
    realOdorSFC(j) = SFCdata.RealSFC.odor;%biased
    realOdorLightSFC(j) = SFCdata.RealSFC.odorLight;
    
    shiftedOdorSFC(j) = SFCdata.ShuffledSFC.odor;%Biased
    shiftedOdorLightSFC(j) = SFCdata.ShuffledSFC.odorLight;
    
end
%% compute the VS
    %---
    %VS - Odor
    realOdorPLV = circ_r(allRealOdorSpikePhases');
    shiftedOdorPLV = circ_r(allShiftedOdorSpikePhases');
    %Rayleigh test
    pOdor = circ_rtest(allRealOdorSpikePhases');
    
    % Mean phase
    meanPhaseOdor = circ_mean(allRealOdorSpikePhases');
    % spike phase dist.
    tmpPhaseHist = hist(allRealOdorSpikePhases, params.phaseBins);
    phaseDistOdor = tmpPhaseHist./length(allRealOdorSpikePhases);
    %mean phase shuffle
    meanPhaseOdorSHUFF = circ_mean(allShiftedOdorSpikePhases');
    
    %VS - Odor+Light
    realOdorLightPLV = circ_r(allRealOdorLightSpikePhases');
    shiftedOdorLightPLV = circ_r(allShiftedOdorLightSpikePhases');
    %Rayleigh test
    pOdorLight = circ_rtest(allRealOdorLightSpikePhases');
    
    % Mean phase
    meanPhaseOdorLight = circ_mean(allRealOdorLightSpikePhases');
     % spike phase dist.
    tmpPhaseHist = hist(allRealOdorLightSpikePhases, params.phaseBins);
    phaseDistOdorLight = tmpPhaseHist./length(allRealOdorLightSpikePhases);
    %mean phase shuffle
    meanPhaseOdorLightSHUFF = circ_mean(allShiftedOdorLightSpikePhases');
    %---
    
%% normalize the PLV
%%%%% Subtract the shift predictor from the biased values %%%
normOdorPLV = (realOdorPLV-shiftedOdorPLV);
normOdorLightPLV = (realOdorLightPLV-shiftedOdorLightPLV);

%%%%% truncate the negative values to zero %%%%
normOdorPLV(find(normOdorPLV < 0)) = 0;
normOdorLightPLV(find(normOdorLightPLV < 0)) = 0;

%% normalized the SFC
%%%%% Subtract the shift predictor from the biased values %%%
normOdorSFC = nanmean(realOdorSFC-shiftedOdorSFC);
normOdorLightSFC = nanmean(realOdorLightSFC-shiftedOdorLightSFC);

%%%%% truncate the negative values to zero %%%%
normOdorSFC(find(normOdorSFC < 0)) = 0;
normOdorLightSFC(find(normOdorLightSFC < 0)) = 0;

allRealOdorSFC = nanmean(realOdorSFC);
allRealOdorLightSFC = nanmean(realOdorLightSFC);

realOdorSFC = [];realOdorLightSFC = [];
shiftedOdorSFC = [];shiftedOdorLightSFC = [];
%% Test if there was a significant odor response
[~,pFROdor] = ttest(firingRateOdor, bgFiringOdor);
[~,pFROdorLight] = ttest(firingRateOdorLight, bgFiringOdorLight);

evokedFiringOdor = mean(firingRateOdor-bgFiringOdor);
evokedFiringOdorLight = mean(firingRateOdorLight-bgFiringOdorLight);

%% %% assign to outputStruct
% PLV
outputPLV.realPLV.odor = realOdorPLV;
outputPLV.realPLV.odorLight = realOdorLightPLV;

outputPLV.normPLV.odor = normOdorPLV;
outputPLV.normPLV.odorLight = normOdorLightPLV;

outputPLV.pRayleighTest.odor = pOdor;
outputPLV.pRayleighTest.odorLight = pOdorLight;

% Firing 
outputPLV.FR.odor = firingRateOdor;
outputPLV.FR.odorLight = firingRateOdorLight;

outputPLV.baselineFR.odor = bgFiringOdor;
outputPLV.baselineFR.odorLight = bgFiringOdorLight;

outputPLV.odorResponsePvalue.odor = pFROdor;
outputPLV.odorResponsePvalue.odorLight = pFROdorLight;

outputPLV.evokedResponse.odor = evokedFiringOdor;
outputPLV.evokedResponse.odorLight = evokedFiringOdorLight;

%Phases
outputPLV.meanPhase.odor = meanPhaseOdor;
outputPLV.meanPhase.odorLight = meanPhaseOdorLight;

outputPLV.meanPhaseShuffled.odor = meanPhaseOdorSHUFF;
outputPLV.meanPhaseShuffled.odorLight = meanPhaseOdorLightSHUFF;

outputPLV.phaseDist.odor = phaseDistOdor;
outputPLV.phaseDist.odorLight = phaseDistOdorLight;

% SFC values
outputPLV.realSFC.odor = allRealOdorSFC;
outputPLV.realSFC.odorLight = allRealOdorLightSFC;

outputPLV.normSFC.odor = normOdorSFC;
outputPLV.normSFC.odorLight = normOdorLightSFC;

% PSD
if computePSD == 1
    outputPLV.PSD = psdStructure;
end



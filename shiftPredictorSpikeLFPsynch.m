function [outputPLV] = shiftPredictorSpikeLFPsynch(matName, methodFlag)
%% load
load(matName)
computePSD = 0;
%% Params
params.bgTimeBounds = [-2 0];
params.fgTimeBounds = [0 2];%Sec
params.decimateFs = expLFP_SPIKES(1).decimateFs;
params.lfpRange = [40 70];%Hz
params.SFCcurveRange = [10 120];%Hz
params.stimType = expLFP_SPIKES(1).stimuliTyps;%0=odor, 1 = odor+light, 2 = light.
params.numChannels = size(expLFP_SPIKES,2);% number of recording channels.
params.numSpikesThreshold = 5;
params.t = expLFP_SPIKES(1).t;% time
params.winSize = 0.05;%50ms window of SFC.
params.phaseBins = deg2rad(-165:30:180);
params.SFCcurveFlag = 1;


%% %%%% FILTERS %%%%
lpFiltLFP = designfilt('lowpassiir','FilterOrder',10, ...
    'PassbandFrequency',params.lfpRange(2),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

hpFiltLFP = designfilt('highpassiir','FilterOrder',8, ...
    'PassbandFrequency',params.lfpRange(1),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

%Filters for SFC curve
lpFiltSFCcurve = designfilt('lowpassiir','FilterOrder',10, ...
    'PassbandFrequency',params.SFCcurveRange(2),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

hpFiltSFCcurve = designfilt('highpassiir','FilterOrder',8, ...
    'PassbandFrequency',params.SFCcurveRange(1),'PassbandRipple',0.2, ...
    'SampleRate',params.decimateFs);

odorTrials = find(params.stimType==0);
odorLightTrials = find(params.stimType==1);


%%
%-------------- spikes per condition ---------
odorLightSpikes = expLFP_SPIKES(1).spikes(odorLightTrials);
odorSpikes = expLFP_SPIKES(1).spikes(odorTrials);

%% % init values
firingRateOdor = nan(1,length(odorSpikes));
firingRateOdorLight = firingRateOdor;

realOdorPLV = firingRateOdor;
realOdorLightPLV = firingRateOdor;

shiftedOdorPLV = firingRateOdor;
shiftedOdorLightPLV = firingRateOdor;

pOdor = firingRateOdor;
pOdorLight = firingRateOdor;

meanPhaseOdor = nan(1,params.numChannels);
meanPhaseOdorLight = meanPhaseOdor;
meanPhaseOdorSHUFF = meanPhaseOdor;
meanPhaseOdorLightSHUFF = meanPhaseOdorLight;

phaseDistOdor = cell(1,params.numChannels);
phaseDistOdorLight = cell(1,params.numChannels);

allNormOdorLightSFC = nan(1,params.numChannels);
allNormOdorSFC = allNormOdorLightSFC;

allNegativeNormOdorSFC = allNormOdorLightSFC;
allNegativeNormOdorLightSFC = allNormOdorLightSFC;

allRealOdorLightSFC = nan(1,params.numChannels);
allRealOdorSFC = allRealOdorLightSFC;

% Init SFC curves
allRealOdorSFCcurve = nan(params.numChannels, params.decimateFs*(params.winSize)+1);
allRealOdorLightSFCcurve = allRealOdorSFCcurve;

% Shuffled values
allShiftedOdorSFCcurve = allRealOdorSFCcurve;
allShiftedOdorLightSFCcurve = allRealOdorSFCcurve;

allSpikePhasesOdor = cell(1,params.numChannels);
allSpikePhasesOdorLight = cell(1,params.numChannels);
%%
for i = 1:params.numChannels %%%% run over all channels
    
    %%%-------------- divide LFP by conditions ------------------
    %Mutliply in -1 to flip the LFP signal again.
    odorlightLFP = expLFP_SPIKES(i).LFP(odorLightTrials,:).*-1;
    odorLFP = expLFP_SPIKES(i).LFP(odorTrials,:).*-1;
    
    
    %--------------- PSD ------------------------------------
    if computePSD == 1
        %- trial duration UN-FILTERED
        unFiteredLFP.FGunFilteredLFPodor = odorLFP(:,find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)));
        unFiteredLFP.FGunFilteredLFPodorLight = odorlightLFP(:,find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)));
        
        %- baseline UN-FILTERED
        unFiteredLFP.BGunFilteredLFPodor = odorLFP(:,find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)));
        unFiteredLFP.BGunFilteredLFPodorLight = odorlightLFP(:,find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)));
        
        [psdStructure] = multiTaperPSD(unFiteredLFP, params);
        
        allChannelsPSD(i).PSD = psdStructure;
    end
    
        
    
    %%%-------------- Phase-locking - filter the LFPs------------------------
    %low pass filter
    odorLightLFPbpf = filtfilt(lpFiltLFP, odorlightLFP');
    odorLFPbpf = filtfilt(lpFiltLFP, odorLFP');
    
    % high pass filter
    odorLightLFPbpf = filtfilt(hpFiltLFP, odorLightLFPbpf);
    odorLFPbpf = filtfilt(hpFiltLFP, odorLFPbpf);
    
    %%%-----FILTER for SFC curve (20-100Hz)---------
    %low pass filter
    odorLightLFPSFCcurve = filtfilt(lpFiltSFCcurve, odorlightLFP');
    odorLFPSFCcurve = filtfilt(lpFiltSFCcurve, odorLFP');
    
    % high pass filter
    odorLightLFPSFCcurve = filtfilt(hpFiltSFCcurve, odorLightLFPSFCcurve);
    odorLFPSFCcurve = filtfilt(hpFiltSFCcurve, odorLFPSFCcurve);
    
    
    %%%------------- extract the trial duration signals------------------
    FGodorLightLFPbpf = odorLightLFPbpf(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);
    FGodorLFPbpf = odorLFPbpf(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);
    
    % ----------for SFC curve (20-100Hz)--------------
    FGodorLightLFPSFCcurve = odorLightLFPSFCcurve(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);
    FGodorLFPSFCcurve = odorLFPSFCcurve(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)),:);
    
    %     %%%------------- extract the baseline duration signals------------------
%     BGodorLightLFPbpf = odorLightLFPbpf(find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)),:);
%     BGodorLFPbpf = odorLFPbpf(find(params.t >= params.bgTimeBounds(1) & params.t < params.bgTimeBounds(2)),:);
    
%     %%%----- Z-score the fg signal %%%%
%     zScoreFGodor = (FGodorLFPbpf-mean(BGodorLFPbpf))./std(BGodorLFPbpf);
%     zScoreFGodorLight = (FGodorLightLFPbpf-mean(BGodorLightLFPbpf))./std(BGodorLightLFPbpf);

    
    %%%------------- Shifted signals (SHUFFLED)------------------
    FGodorLightLFPbpfSHIFTED = circshift(FGodorLightLFPbpf,1,2);% shift of 1 trial
    FGodorLFPbpfSHIFTED = circshift(FGodorLFPbpf, 1, 2);% shift of 1 trial
    
    % ----------for SFC curve (20-100Hz)--------------
    FGodorLightLFPbpfSHIFTEDSFCcurve = circshift(FGodorLightLFPSFCcurve,1,2);% shift of 1 trial
    FGodorLFPbpfSHIFTEDSFCcurve = circshift(FGodorLFPSFCcurve, 1, 2);% shift of 1 trial
    
    %%%% init
    allShiftedOdorSpikePhases = [];
    allRealOdorSpikePhases = [];
    allShiftedOdorLightSpikePhases = [];
    allRealOdorLightSpikePhases = [];
    
    %Values
    realOdorSFC = nan(1,size(FGodorLFPbpf,2));
    realOdorLightSFC = realOdorSFC;
    shiftedOdorSFC = realOdorSFC;
    shiftedOdorLightSFC = realOdorSFC;
    
    %Curves
    realOdorSFCcurve = nan(size(FGodorLFPbpf,2), params.decimateFs*(params.winSize)+1);
    realOdorLightSFCcurve = realOdorSFCcurve;
    shiftedOdorSFCcurve = realOdorSFCcurve;
    shiftedOdorLightSFCcurve = realOdorSFCcurve;
    
    % STA
    odorSTA = nan(size(FGodorLFPbpf,2), params.decimateFs*(params.winSize*2)+1);
    odorLightSTA = nan(size(FGodorLFPbpf,2), params.decimateFs*(params.winSize*2)+1);
   
    
    for j = 1:size(FGodorLFPbpf,2)%%%%% run over trials %%%%
        
        %% curr trial data
        
        spikeLFPdata.spikes.odor = odorSpikes{j};
        spikeLFPdata.spikes.odorLight = odorLightSpikes{j};
        
        %------------- compute baseline firing rate ---------
        if i==1
            % baseline firing rate
            bgFiringOdor(j) = length(spikeLFPdata.spikes.odor(find(spikeLFPdata.spikes.odor >= params.bgTimeBounds(1) ...
                & spikeLFPdata.spikes.odor < params.bgTimeBounds(2))))./diff(params.fgTimeBounds);
            
            bgFiringOdorLight(j) = length(spikeLFPdata.spikes.odorLight(find(spikeLFPdata.spikes.odorLight >= params.bgTimeBounds(1) ...
                & spikeLFPdata.spikes.odorLight < params.bgTimeBounds(2))))./diff(params.fgTimeBounds);
        end
        
        %%%------------- extract the trial duration of spikes------------------
        spikeLFPdata.spikes.odor = spikeLFPdata.spikes.odor(find(spikeLFPdata.spikes.odor >= params.fgTimeBounds(1) & spikeLFPdata.spikes.odor < params.fgTimeBounds(2)));
        spikeLFPdata.spikes.odorLight = spikeLFPdata.spikes.odorLight(find(spikeLFPdata.spikes.odorLight >= params.fgTimeBounds(1) & spikeLFPdata.spikes.odorLight < params.fgTimeBounds(2)));
        
        spikeLFPdata.realLFP.odor = FGodorLFPbpf(:,j);
        spikeLFPdata.realLFP.odorLight = FGodorLightLFPbpf(:,j);
        
        spikeLFPdata.shiftedLFP.odor = FGodorLFPbpfSHIFTED(:,j);
        spikeLFPdata.shiftedLFP.odorLight = FGodorLightLFPbpfSHIFTED(:,j);
        
        % ----------for SFC curve (20-100Hz)--------------
        spikeLFPdata.LFPforSFCcurve.odor = FGodorLFPSFCcurve(:,j);
        spikeLFPdata.LFPforSFCcurve.odorLight = FGodorLightLFPSFCcurve(:,j);
        
        spikeLFPdata.shiftedLFPforSFCcurve.odor = FGodorLFPbpfSHIFTEDSFCcurve(:,j);
        spikeLFPdata.shiftedLFPforSFCcurve.odorLight = FGodorLightLFPbpfSHIFTEDSFCcurve(:,j);
        

        if i==1  
            % firing during stimuli presentation
            firingRateOdor(j) = length(spikeLFPdata.spikes.odor)./diff(params.fgTimeBounds);
            firingRateOdorLight(j) = length(spikeLFPdata.spikes.odorLight)./diff(params.fgTimeBounds);
        end
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
        
        % Assign SFC values
        realOdorSFC(j) = SFCdata.RealSFC.odor;%biased
        realOdorLightSFC(j) = SFCdata.RealSFC.odorLight;
        
        shiftedOdorSFC(j) = SFCdata.ShuffledSFC.odor;%Biased
        shiftedOdorLightSFC(j) = SFCdata.ShuffledSFC.odorLight;
        
        % Save the STA per trial
        odorSTA(j,:) = SFCdata.RealSTA.odor;
        odorLightSTA(j,:) = SFCdata.RealSTA.odorLight;
        
        % Assign SFC curves
        % Real values
        realOdorSFCcurve(j,:) = SFCdata.RealSFCcurve.odor;
        realOdorLightSFCcurve(j,:) = SFCdata.RealSFCcurve.odorLight;
        
        % Shuffled values
        shiftedOdorSFCcurve(j,:) = SFCdata.ShuffledSFCcurve.odor;
        shiftedOdorLightSFCcurve(j,:) = SFCdata.ShuffledSFCcurve.odorLight;
        
        %Frequency vector
        SFCcurvesFrequencyVector = SFCdata.SFCcurveFrequency;
    end
    %% compute the VS
    if isempty(allRealOdorSpikePhases) | isempty(allRealOdorLightSpikePhases)
        break;
    else
        %VS
        realOdorPLV(i) = circ_r(allRealOdorSpikePhases');
        shiftedOdorPLV(i) = circ_r(allShiftedOdorSpikePhases');
        %Rayleigh test
        pOdor(i) = circ_rtest(allRealOdorSpikePhases');
        
        % Mean phase
        meanPhaseOdor(i) = circ_mean(allRealOdorSpikePhases');
        tmpPhaseHist = hist(allRealOdorSpikePhases, params.phaseBins);
        phaseDistOdor{i} = tmpPhaseHist./length(allRealOdorSpikePhases);
        
        allSpikePhasesOdor{i} = allRealOdorSpikePhases;
        
        meanPhaseOdorSHUFF(i) = circ_mean(allShiftedOdorSpikePhases');
        
        %VS
        realOdorLightPLV(i) = circ_r(allRealOdorLightSpikePhases');
        shiftedOdorLightPLV(i) = circ_r(allShiftedOdorLightSpikePhases');
        %Rayleigh test
        pOdorLight(i) = circ_rtest(allRealOdorLightSpikePhases');
        % Mean phase
        meanPhaseOdorLight(i) = circ_mean(allRealOdorLightSpikePhases');
        tmpPhaseHist = hist(allRealOdorLightSpikePhases, params.phaseBins);
        phaseDistOdorLight{i} = tmpPhaseHist./length(allRealOdorLightSpikePhases);
        
        allSpikePhasesOdorLight{i} = allRealOdorLightSpikePhases;
        
        meanPhaseOdorLightSHUFF(i) = circ_mean(allShiftedOdorLightSpikePhases');
    end
    
    %% normalized the SFC
    %%%%% Subtract the shift predictor from the biased values %%%
    normOdorSFC = nanmean(realOdorSFC-shiftedOdorSFC);
    normOdorLightSFC = nanmean(realOdorLightSFC-shiftedOdorLightSFC);
    
    negativeNormOdorSFC = normOdorSFC;
    negativeNormOdorLightSFC = normOdorLightSFC;
    
    %%%%% truncate the negative values to zero %%%%
    normOdorSFC(find(normOdorSFC < 0)) = 0;
    normOdorLightSFC(find(normOdorLightSFC < 0)) = 0;
    
    
    allNormOdorSFC(i) = normOdorSFC;
    allNormOdorLightSFC(i) = normOdorLightSFC;
    
    allNegativeNormOdorSFC(i) = negativeNormOdorSFC;
    allNegativeNormOdorLightSFC(i) = negativeNormOdorLightSFC;
    
    allRealOdorSFC(i) = nanmean(realOdorSFC);
    allRealOdorLightSFC(i) = nanmean(realOdorLightSFC);
    
    %----Average all SFC curves across trials -----
    % Real values
    allRealOdorSFCcurve(i,:) = nanmean(realOdorSFCcurve);
    allRealOdorLightSFCcurve(i,:) = nanmean(realOdorLightSFCcurve);
    
    % Shuffled values
    allShiftedOdorSFCcurve(i,:) = nanmean(shiftedOdorSFCcurve);
    allShiftedOdorLightSFCcurve(i,:) = nanmean(shiftedOdorLightSFCcurve);
    
    % clear variables
    realOdorSFC = [];realOdorLightSFC = [];
    shiftedOdorSFC = [];shiftedOdorLightSFC = [];
    realOdorSFCcurve = [];realOdorLightSFCcurve = [];
    shiftedOdorSFCcurve = [];shiftedOdorLightSFCcurve = [];
    


end
%% normalize the PLV
%%%%% Subtract the shift predictor from the biased values %%%
normOdorPLV = (realOdorPLV-shiftedOdorPLV);
normOdorLightPLV = (realOdorLightPLV-shiftedOdorLightPLV);

negativeNormOdorPLV = normOdorPLV;
negativeNormOdorLightPLV = normOdorLightPLV;

%%%%% truncate the negative values to zero %%%%
normOdorPLV(find(normOdorPLV < 0)) = 0;
normOdorLightPLV(find(normOdorLightPLV < 0)) = 0;

%% test the response to odor
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

outputPLV.negativeNormPLV.odor = negativeNormOdorPLV;
outputPLV.negativeNormPLV.odorLight = negativeNormOdorLightPLV;

outputPLV.pRayleighTest.odor = pOdor;
outputPLV.pRayleighTest.odorLight = pOdorLight;

outputPLV.meanPhase.odor = meanPhaseOdor;
outputPLV.meanPhase.odorLight = meanPhaseOdorLight;

outputPLV.meanPhaseShuffled.odor = meanPhaseOdorSHUFF;
outputPLV.meanPhaseShuffled.odorLight = meanPhaseOdorLightSHUFF;

outputPLV.phaseDist.odor = phaseDistOdor;
outputPLV.phaseDist.odorLight = phaseDistOdorLight;

outputPLV.allSpikePhases.odor = allSpikePhasesOdor;
outputPLV.allSpikePhases.odorLight = allSpikePhasesOdorLight;


% SFC values
outputPLV.realSFC.odor = allRealOdorSFC;
outputPLV.realSFC.odorLight = allRealOdorLightSFC;

outputPLV.normSFC.odor = allNormOdorSFC;
outputPLV.normSFC.odorLight = allNormOdorLightSFC;

outputPLV.negativeNormSFC.odor = allNegativeNormOdorSFC;
outputPLV.negativeNormSFC.odorLight = allNegativeNormOdorLightSFC;


% SFC curves
% Real values
outputPLV.realSFCcurve.Odor = allRealOdorSFCcurve;
outputPLV.realSFCcurve.OdorLight = allRealOdorLightSFCcurve;

% Shuffled values
outputPLV.shuffledSFCcurve.Odor = allShiftedOdorSFCcurve;
outputPLV.shuffledSFCcurve.OdorLight = allShiftedOdorLightSFCcurve;
outputPLV.SFCcurvesFrequencyVector = SFCcurvesFrequencyVector;


% Firing
outputPLV.FR.odor = firingRateOdor;
outputPLV.FR.odorLight = firingRateOdorLight;

outputPLV.baselineFR.odor = bgFiringOdor;
outputPLV.baselineFR.odorLight = bgFiringOdorLight;

outputPLV.odorResponsePvalue.odor = pFROdor;
outputPLV.odorResponsePvalue.odorLight = pFROdorLight;

outputPLV.evokedResponse.odor = evokedFiringOdor;
outputPLV.evokedResponse.odorLight = evokedFiringOdorLight;

% PSD
if computePSD == 1
    outputPLV.PSD = allChannelsPSD;
end


% trialNum = 1;
% %%%%------Theta range filters------------
% ThetaRange = [2 10];
% bpFiltTheta = designfilt('bandpassiir','FilterOrder',20, ...
%          'HalfPowerFrequency1',ThetaRange(1),'HalfPowerFrequency2',ThetaRange(2), ...
%          'SampleRate',params.decimateFs);
% odorLightTheta = filtfilt(bpFiltTheta, odorlightLFP');
% odorTheta = filtfilt(bpFiltTheta, odorLFP');
% 
% %%%%------Gamma range filters------------
% odorLightLFPbpf = filtfilt(lpFiltLFP, odorlightLFP');
% odorLFPbpf = filtfilt(lpFiltLFP, odorLFP');
% 
% % high pass filter
% odorLightLFPbpf = filtfilt(hpFiltLFP, odorLightLFPbpf);
% odorLFPbpf = filtfilt(hpFiltLFP, odorLFPbpf);
% 
% figure;
% %-----ODOR plots----
% subplot(3,2,1)
% plot(params.t, odorLFP(trialNum, :), 'k')
% title('Odor')
% xlim([-2 2])
% subplot(3,2,3)
% plot(params.t, odorTheta(:,trialNum), 'k')
% xlim([-2 2])
% subplot(3,2,5)
% plot(params.t, odorLFPbpf(:,trialNum), 'k')
% xlim([-2 2])
% %-----ODOR+Light plots----
% subplot(3,2,2)
% plot(params.t, odorlightLFP(trialNum, :), 'b')
% title('Odor')
% xlim([-2 2])
% subplot(3,2,4)
% plot(params.t, odorLightTheta(:,trialNum), 'b')
% xlim([-2 2])
% subplot(3,2,6)
% plot(params.t, odorLightLFPbpf(:,trialNum), 'b')
% xlim([-2 2])


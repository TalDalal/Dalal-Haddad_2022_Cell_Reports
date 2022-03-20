function [SFCdata] = computeSpikeFieldCoherenceShiftPredictor(spikeLFPdata, params)

%---------------- ODOR CONDITION -------------------
%---------------- Real LFP ODOR ---------------
[STAodor, cohCurveOdorReal,~,~,~,~,f, ~] = ...
    getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.realLFP.odor, params.decimateFs, params.winSize,0 , 0);

f = ceil(f);% To have integer frequencies
% average SFC in the range of interest.
SFCdata.RealSFC.odor =  mean(cohCurveOdorReal(find(f >= params.lfpRange(1) & f <=  params.lfpRange(2))));
%Save the STA
SFCdata.RealSTA.odor = STAodor;

%----------------- SFC curve - wide frequency range ----------------------
if params.SFCcurveFlag > 0
    [~, cohCurveWideFreqRangeOdor,~,~,~,~,fSFCcurve, ~] = ...
        getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.LFPforSFCcurve.odor, params.decimateFs, params.winSize,0 , 0);% double the window size for better F resolution
    
    SFCdata.RealSFCcurve.odor = cohCurveWideFreqRangeOdor;
end


%---------------- Shuffled LFP ODOR ---------------
[~, cohCurveOdorShuffled,~,~,~,~,f, ~] = ...
    getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.shiftedLFP.odor, params.decimateFs, params.winSize,0 , 0);
f = ceil(f);% To have integer frequencies
% average SFC in the range of interest.
SFCdata.ShuffledSFC.odor =  mean(cohCurveOdorShuffled(find(f >= params.lfpRange(1) & f <=  params.lfpRange(2))));


%--For shuffled SFC curve
if params.SFCcurveFlag > 0
    [~, cohCurveWideFreqRangeOdorShuffled,~,~,~,~,fSFCcurve, ~] = ...
        getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.shiftedLFPforSFCcurve.odor, params.decimateFs, params.winSize,0 , 0);% double the window size for better F resolution
    
    SFCdata.ShuffledSFCcurve.odor = cohCurveWideFreqRangeOdorShuffled;
end
%% ---------------- ODOR + Light CONDITION -------------------
%---------------- Real LFP ---------------
[STAodorLight, cohCurveOdorLightReal,~,~,~,~,f, ~] = ...
    getCoherence(spikeLFPdata.spikes.odorLight, spikeLFPdata.realLFP.odorLight, params.decimateFs, params.winSize,0 , 0);
f = ceil(f);% To have integer frequencies
% average SFC in the range of interest.
SFCdata.RealSFC.odorLight =  mean(cohCurveOdorLightReal(find(f >= params.lfpRange(1) & f <=  params.lfpRange(2))));
%Save the STA
SFCdata.RealSTA.odorLight = STAodorLight;

%----------------- SFC curve - wide frequency range ----------------------
if params.SFCcurveFlag > 0
    [~, cohCurveWideFreqRangeOdorLight,~,~,~,~,fSFCcurve, ~] = ...
        getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.LFPforSFCcurve.odorLight, params.decimateFs, params.winSize,0 , 0);% double the window size for better F resolution
    
    SFCdata.RealSFCcurve.odorLight = cohCurveWideFreqRangeOdorLight;
end

%---------------- Shuffled LFP ---------------

[~, cohCurveOdorLightShuffled,~,~,~,~,f, ~] = ...
    getCoherence(spikeLFPdata.spikes.odorLight, spikeLFPdata.shiftedLFP.odorLight, params.decimateFs, params.winSize,0 , 0);
f = ceil(f);% To have integer frequencies
% average SFC in the range of interest.
SFCdata.ShuffledSFC.odorLight =  mean(cohCurveOdorLightShuffled(find(f >= params.lfpRange(1) & f <=  params.lfpRange(2))));


%--For shuffled SFC curve
if params.SFCcurveFlag > 0
    [~, cohCurveWideFreqRangeOdorLightShuffled,~,~,~,~,fSFCcurve, ~] = ...
        getCoherence(spikeLFPdata.spikes.odor, spikeLFPdata.shiftedLFPforSFCcurve.odorLight, params.decimateFs, params.winSize,0 , 0);% double the window size for better F resolution
    
    SFCdata.ShuffledSFCcurve.odorLight = cohCurveWideFreqRangeOdorLightShuffled;
    SFCdata.SFCcurveFrequency = fSFCcurve;
end

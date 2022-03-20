
function [psdStructure] = multiTaperPSD(unFiteredLFP, params)

%The output of this function is the average power spectrum across trials
%for both conditions, together with the relative power FG/BG in the
%requested frequency band.


%%%%% multi taper power estimation %%%%
TW = 5;%TW = (N(seconds)*dFreq(spectral resolution))/2
%L = number of tapers, L = (2*TW) -1.
%The output 'Freq' is 1/signalLength in Hz.

%In single elctrode experiments where we suffered line noise
removeLineFrequencies = [40 53];%Hz
%% ---------------------------- ODOR ----------------------------
%Multi Taper method
%-------Baseline ------------
[PxxBG,Freq] = pmtm(unFiteredLFP.BGunFilteredLFPodor',TW,length(unFiteredLFP.BGunFilteredLFPodor),params.decimateFs);
%integrate the bg PSD
if exist('removeLineFrequencies', 'var')
    freqNoLineNoise = Freq(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)));
    %Remove the infected frequencies
    freqNoLineNoise = setdiff(freqNoLineNoise, removeLineFrequencies(1):Freq(2)-Freq(1):removeLineFrequencies(2));
    bgBandVals = trapz(PxxBG(find(Freq>=freqNoLineNoise(1) & Freq < freqNoLineNoise(end)),:));
else
    bgBandVals = trapz(PxxBG(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)),:));
end

%------trial duration-----------
[PxxFG, Freq] = pmtm(unFiteredLFP.FGunFilteredLFPodor',TW,length(unFiteredLFP.FGunFilteredLFPodor), params.decimateFs);
%integrate the fg PSD
if exist('removeLineFrequencies', 'var')
    freqNoLineNoise = Freq(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)));
    %Remove the infected frequencies
    freqNoLineNoise = setdiff(freqNoLineNoise, removeLineFrequencies(1):Freq(2)-Freq(1):removeLineFrequencies(2));
    fgBandVals = trapz(PxxFG(find(Freq>=freqNoLineNoise(1) & Freq < freqNoLineNoise(end)),:));
else
    fgBandVals = trapz(PxxFG(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)),:));
end
psdStructure.Freq = Freq;

% Relative power (divide the fg by the bg).
 psdStructure.relativePower.odor = mean(fgBandVals./bgBandVals);

%%%%%%%%%%% save the power spectrum
% mean PSC curve across trials. FG
psdStructure.FGspectrum.odor = mean(PxxFG');

% mean PSC curve across trials. BG
psdStructure.BGspectrum.odor = mean(PxxBG');

% clear the power spectrum
fgBandVals = [];bgBandVals = [];
PxxFG = [];PxxBG = [];

%% ---------------------------- +ChR2 ----------------------------
%Multi Taper method
%-------Baseline ------------
[PxxBG,Freq] = pmtm(unFiteredLFP.BGunFilteredLFPodorLight',TW,length(unFiteredLFP.BGunFilteredLFPodorLight),params.decimateFs);
%integrate the bg PSD
if exist('removeLineFrequencies', 'var')
    freqNoLineNoise = Freq(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)));
    %Remove the infected frequencies
    freqNoLineNoise = setdiff(freqNoLineNoise, removeLineFrequencies(1):Freq(2)-Freq(1):removeLineFrequencies(2));
    bgBandVals = trapz(PxxBG(find(Freq>=freqNoLineNoise(1) & Freq < freqNoLineNoise(end)),:));
else
    bgBandVals = trapz(PxxBG(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)),:));
end

%------trial duration-----------
[PxxFG, Freq] = pmtm(unFiteredLFP.FGunFilteredLFPodorLight',TW,length(unFiteredLFP.FGunFilteredLFPodorLight), params.decimateFs);
%integrate the fg PSD
if exist('removeLineFrequencies', 'var')
    freqNoLineNoise = Freq(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)));
    %Remove the infected frequencies
    freqNoLineNoise = setdiff(freqNoLineNoise, removeLineFrequencies(1):Freq(2)-Freq(1):removeLineFrequencies(2));
    fgBandVals = trapz(PxxFG(find(Freq>=freqNoLineNoise(1) & Freq < freqNoLineNoise(end)),:));
else
    fgBandVals = trapz(PxxFG(find(Freq>=params.lfpRange(1) & Freq < params.lfpRange(2)),:));
end

% Relative power (divide the fg by the bg).
 psdStructure.relativePower.odorLight = mean(fgBandVals./bgBandVals);

%%%%%%%%%%% save the power spectrum
% mean PSC curve across trials. FG
psdStructure.FGspectrum.odorLight = mean(PxxFG');

% mean PSC curve across trials. BG
psdStructure.BGspectrum.odorLight = mean(PxxBG');

% clear the power spectrum
fgBandVals = [];bgBandVals = [];
PxxFG = [];PxxBG = [];


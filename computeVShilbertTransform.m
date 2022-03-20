function [spikePhaseStruct] = computeVShilbertTransform(spikeLFPdata, params)


%%%%%% odor  spike train %%%%
tFG = params.t(find(params.t >= params.fgTimeBounds(1) & params.t < params.fgTimeBounds(2)));
stOdor = zeros(1,length(tFG));
for i = 1:length(stOdor)-1
    if ~isempty(find(spikeLFPdata.spikes.odor >= tFG(i) & spikeLFPdata.spikes.odor <= tFG(i+1)))
        stOdor(i) = 1;
    end
end

%%----------------- REAL LFP --------------------------------
%%%% Hilbert transform for inst. phase %%
phaseTMP = hilbert(spikeLFPdata.realLFP.odor);
realPhaseOdor = angle(phaseTMP);
%%%% assign the phases
spikePhaseStruct.realOdor = realPhaseOdor(logical(stOdor))';
%%%% assign the power of each phase
spikePhaseStruct.realOdorPhasePower = abs(phaseTMP(logical(stOdor)));
phaseTMP = [];


%%----------------- SHIFTED VS --------------------------------
%%%% Hilbert transform for inst. phase %%
phaseTMP = hilbert(spikeLFPdata.shiftedLFP.odor);
shiftedPhaseOdor = angle(phaseTMP);

%%%% assign the phases
spikePhaseStruct.shiftedOdor = shiftedPhaseOdor(logical(stOdor))';
%%%% assign the power of each phase
spikePhaseStruct.shiftedOdorPhasePower = abs(phaseTMP(logical(stOdor)));
phaseTMP = [];

%% ----------------------- odor + Light condition -------------------------

%%%%%% odor + light  spike train %%%%
stOdorLight = zeros(1,length(tFG));
for i = 1:length(stOdorLight)-1
    if ~isempty(find(spikeLFPdata.spikes.odorLight >= tFG(i) & spikeLFPdata.spikes.odorLight <= tFG(i+1)))
        stOdorLight(i) = 1;
    end
end

%%----------------- REAL LFP --------------------------------
%%%% Hilbert transform for inst. phase %%
phaseTMP = hilbert(spikeLFPdata.realLFP.odorLight);
realPhaseOdorLight = angle(phaseTMP);

%%%% assign the phases
spikePhaseStruct.realOdorLight = realPhaseOdorLight(logical(stOdorLight))';
%%%% assign the power of each phase
spikePhaseStruct.realOdorLightPhasePower = abs(phaseTMP(logical(stOdorLight)));
phaseTMP = [];

%%----------------- SHIFTED VS --------------------------------
%%%% Hilbert transform for inst. phase %%
phaseTMP = hilbert(spikeLFPdata.shiftedLFP.odorLight);
shiftedPhaseOdorLight = angle(phaseTMP);

%%%% assign the phases
spikePhaseStruct.shiftedOdorLight = shiftedPhaseOdorLight(logical(stOdorLight))';
%%%% assign the power of each phase
spikePhaseStruct.shiftedOdorLightPhasePower = abs(phaseTMP(logical(stOdorLight)));
phaseTMP = [];


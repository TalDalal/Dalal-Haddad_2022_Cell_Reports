function [] = runShiftPredictorSynchAllNeurons(mouseState, methodFlag, savingPath)

%**** 'methodFlag': 'peakDetection' or 'hilbert'
%% CD
if strcmpi(mouseState, 'Anesth')
    
    % data is found in '/Volumes/My
    % Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2/Neuronexus_firstExp'.
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2/Neuronexus_firstExp')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2\Neuronexus_firstExp')
    end
    
    %load the PSD structure
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data/Awake_data_updated/vsShiftPredictorAnalysis/ChR2/anesthPSDstruct.mat');
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data\Awake_data_updated\vsShiftPredictorAnalysis\ChR2\anesthPSDstruct.mat');
    end
    
elseif strcmpi(mouseState, 'Awake')
    
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data')
    end
    
    %load the PSD structure
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/Neuronexus_GC_ChR2_Awake_02122020/Neuronexus_GC_ChR2_Awake_02122020_data/Awake_data_updated/vsShiftPredictorAnalysis/ChR2/awakePSDstruct.mat');
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\Neuronexus_GC_ChR2_Awake_02122020\Neuronexus_GC_ChR2_Awake_02122020_data\Awake_data_updated\vsShiftPredictorAnalysis\ChR2\awakePSDstruct.mat');
    end
    
end
    
%% Load  
load('allNeuronsSFCstruct')
%% Run over the data
for i = 1:size(allNeuronsSFCstruct,2)
    
    matName = allNeuronsSFCstruct(i).name;
        
    % analyze the SFC and other parameters
    [outputPLV] = shiftPredictorSpikeLFPsynch(matName, methodFlag);
    
    %save the data in struct
    allPLV(i).data = outputPLV;
    allPLV(i).expName = matName;
    
    
    %% Add the PSD data to the structure
    if strcmpi(mouseState, 'Anesth')
        allPLV(i).data.PSD = anesthPSDstruct(i).data.PSD.PSD;
        allPLV(i).recordingElectrodeIndx = 'nan';
        
    elseif strcmpi(mouseState, 'Awake')
        allPLV(i).data.PSD = awakePSDstruct(i).data.PSD.PSD;
        % the electrode that recorded the current spikes.
        allPLV(i).recordingElectrodeIndx = awakePSDstruct(i).recordingElectrodeIndx;
    end
    
    %%
    clc
    fprintf([' Analyzed ' num2str(i) '/' num2str(size(allNeuronsSFCstruct,2)) 'neurons \n'])

end
%% save
if ~isempty(savingPath)
    fprintf([' Finished the analysis, now saving the data to: \n ' savingPath ])
    %%% save the data structure.
    save(savingPath, 'allPLV');
end


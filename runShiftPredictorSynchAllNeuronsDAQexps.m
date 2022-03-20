function [] = runShiftPredictorSynchAllNeuronsDAQexps(mouseType, methodFlag, savingPath)

%**** 'methodFlag': 'peakDetection' or 'hilbert'
%% CD
% ^^^^^^^^^^^^^^GC_ArchT data is found in ^^^^^^^^^^
% :'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M48GADCre_01092020\LFP_SPIKES_DATA'
%(mouseType=='ArchT')

% ^^^^^^^^^^^^^^Tbet_ArchT data is found in ^^^^^^^^^^
% :'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M4TbetNpHr_ 11082020\LFP_SPIKES_DATA'
%(mouseType=='Tbet_ArchT')

if strcmpi(mouseType,'ArchT')
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/M48GADCre_01092020/LFP_SPIKES_DATA')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M48GADCre_01092020\LFP_SPIKES_DATA')
    end
elseif strcmpi(mouseType, 'Tbet_ArchT')
    if ismac
        cd('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/M4TbetNpHr_ 11082020/LFP_SPIKES_DATA')
    else
        cd('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\M4TbetNpHr_ 11082020\LFP_SPIKES_DATA')
    end
end
    
%% Load  
load('allNeuronsSFCstruct')
%% Run over the data
for i = 1:size(allNeuronsSFCstruct,2)
    matName = allNeuronsSFCstruct(i).name;
        
    % analyze the SFC and other parameters
    [outputPLV] = shiftPredictorSpikeLFPsynchDAQexps(matName, methodFlag);
    
    %save the data in struct
    allPLV(i).data = outputPLV;
    allPLV(i).expName = matName;
    
    fprintf([' Analyzed ' num2str(i) '/' num2str(size(allNeuronsSFCstruct,2)) ' neurons \n'])
    
end
%% save
if ~isempty(savingPath)
    fprintf([' Finished the analysis, now saving the data to: \n ' savingPath ])
    %%% save the data structure.
    save(savingPath, 'allPLV');
end


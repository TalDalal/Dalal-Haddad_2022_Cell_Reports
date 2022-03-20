%%%% Description
%VirusType - ChR2 or ArchT or tbet-halo.
% It plots the firing rate during baseline, 1st sniff and odor duration.

function [] = scatterMTCfiringBeforeAfterLight(virusType)
%% LOAD %%%%%% OB data %%%%%%%%%

%%%%%%%%%%%% 1=add the firing rate from the single-electrode LFP experimets (M48)%%%%%%%%%%%%%%%
combineFRdataLFPexperiments = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(virusType,'chr2')%%%%%%chr2
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOBChR2_12032020.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOBChR2_12032020.mat')
    end
    data = dataStructOBChR2;
    
elseif strcmpi(virusType,'arch')%%%%%%arch

    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOB_09032020.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOB_09032020.mat')
    end
    
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOBInhibited_12032020.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOBInhibited_12032020.mat')
    end
    
    data = [dataStructOB dataStructOBInhibited];% combine both data sets.
    
elseif strcmpi(virusType,'tbet-halo')%%%%%%arch
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOB_MTarchT.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOB_MTarchT.mat')
    end
    data = dataStructOB_MTarchT;
end
%init params



allOdorDurationEVOKED = [];
allBaselineDurationEVOKED = [];
deltaFR = [];
%% run over the db and compute the FR

for i = 1:size(data,2)
    
    %%%%odor onset%%%%

    %EVOKED FIRING
    evOdor = mean(data(i).odorOnset(2).BgFg(:,2)-mean(data(i).odorOnset(2).BgFg(:,1)));
    evOdorLight = mean(data(i).odorOnset(1).BgFg(:,2)-mean(data(i).odorOnset(1).BgFg(:,1)));
    allOdorDurationEVOKED(i,:) = [evOdor evOdorLight];
    
    deltaFR(i) = mean(data(i).odorOnset(1).BgFg(:,2)) - mean(data(i).odorOnset(2).BgFg(:,2));
    
    
    %%%% baseline
    if size(data(i).odorOnset,2)>2
        allBaselineDurationEVOKED(end+1,:) = [mean(data(i).odorOnset(3).BgFg(:,1)) mean(data(i).odorOnset(3).BgFg(:,2))];
    end
    

end

%% pool the firing rate with those from the LFP experiments.
%% Load the struct
if combineFRdataLFPexperiments==1 & strcmpi(virusType,'arch')
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureArchTMTC_for_CV_analysis.mat')
    else
        
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureArchTMTC_for_CV_analysis.mat')
    end

    
elseif combineFRdataLFPexperiments==1 & strcmpi(virusType,'tbet-halo')
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureNpHRMTC_for_CV_analysis.mat')
    else
        
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureNpHRMTC_for_CV_analysis.mat')
    end

        
elseif combineFRdataLFPexperiments==1 & strcmpi(virusType,'chr2')

    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureChR2MTC_for_CV_analysis.mat')
    else
        
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureChR2MTC_for_CV_analysis.mat')
    end
    
    % a vector indicating anesth (0) or awake recording (1)
    allMouseStateIndex = [zeros(1, size(allOdorDurationEVOKED,1)) stateIndex];
end

%% Run over the data to compute the evoked firing
for lfpExpInd = 1:size(FRstructMTC,2)
    allOdorDurationEVOKED(end+1,:) = [mean(FRstructMTC(lfpExpInd).odorOnset(2).BgFg(:,2)-mean(FRstructMTC(lfpExpInd).odorOnset(2).BgFg(:,1))) ... 
        mean(FRstructMTC(lfpExpInd).odorOnset(1).BgFg(:,2)-mean(FRstructMTC(lfpExpInd).odorOnset(1).BgFg(:,1)))];
    
    deltaFR(end+1) = mean(FRstructMTC(lfpExpInd).odorOnset(1).BgFg(:,2)) - mean(FRstructMTC(lfpExpInd).odorOnset(2).BgFg(:,2)); 
end
    
%% plot
limeLine = [-25 50];
figure;
subplot(1,3,2)
 hold all
%For both awake and anesh recordings
if exist('allMouseStateIndex', 'var')
    % scatter anesth recordings
%     scatter(allOdorDurationEVOKED(find(allMouseStateIndex==0),1), allOdorDurationEVOKED(find(allMouseStateIndex==0),2), 12,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
scatter(allOdorDurationEVOKED(find(allMouseStateIndex==0),1), allOdorDurationEVOKED(find(allMouseStateIndex==0),2), 70,[.5 .5 .5]);
    %scatter the awake recordings
   
    scatter(allOdorDurationEVOKED(find(allMouseStateIndex==1),1), allOdorDurationEVOKED(find(allMouseStateIndex==1),2), 70,'MarkerEdgeColor','k');
else
    scatter(allOdorDurationEVOKED(:,1), allOdorDurationEVOKED(:,2), 70,[.5 .5 .5]);
end
set(gca,'fontSize',14)
box off;
xlabel('Evoked odor response (Hz)')
ylabel('Evoked odor+light response (Hz)')
set(gca, 'YColor', 'k');set(gca, 'XColor', 'k');

plot(limeLine, limeLine, 'k--')
title('Odor stimulation')
xlim(limeLine);ylim(limeLine)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'YColor', 'k');set(gca, 'XColor', 'k')
hold off;

subplot(1,3,1)
histogram(diff(allBaselineDurationEVOKED'),10, 'Normalization', 'probability', 'FaceColor', 'm')
set(gca,'fontSize',14)
box off;
ylabel('Rate of MTC')
xlabel('\Delta firing rate (light only - baseline)')
hold on;
title('Baseline')
set(gcf, 'Renderer', 'painters')
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'YColor', 'k');set(gca, 'XColor', 'k')

if combineFRdataLFPexperiments==1 & strcmpi(virusType,'chr2')
    subplot(1,3,3)
    awakeDeltaFR = deltaFR(allMouseStateIndex==1)';
    anesthDeltaFR = deltaFR(allMouseStateIndex==0)';
    
    barwitherr([std(awakeDeltaFR)./sqrt(length(awakeDeltaFR)) std(anesthDeltaFR)./sqrt(length(anesthDeltaFR))], ...
        [mean(awakeDeltaFR) mean(anesthDeltaFR)])
    
    box off;
ylabel('\Delta firing rate (Hz)')
set(gcf, 'Renderer', 'painters')
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'YColor', 'k');set(gca, 'XColor', 'k')
set(gca, 'XTickLabel', {'Awake', 'Anesth.'}); xtickangle(45)

end
[~,pOdorStim] = ttest(deltaFR);
[~,pBaseline] = ttest(allBaselineDurationEVOKED(:,1), allBaselineDurationEVOKED(:,2));
disp(['MTC odor stimulation, P = : ' num2str(pOdorStim) ', Baseline, P = ' num2str(pBaseline)])

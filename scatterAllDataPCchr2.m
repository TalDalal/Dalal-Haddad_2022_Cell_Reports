function [dataInh, dataExc] = scatterAllDataPCchr2(thType, virusType, singleORpanelOdors)
%If thType==0, when odor or odor+light were significant
%If thType==1, only when odor was significant
%If thType==2, all data.
%% load
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% chr2 %%%%%%%%%%%%
if strcmpi(virusType, 'chr2')
    %%%% single odors
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_singleOdor.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_singleOdor.mat')
    end
    % awake-head fixed cell-odor pairs
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/PC_ChR2_AWAKE_FIRING_RATE_20_12_20')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\PC_ChR2_AWAKE_FIRING_RATE_20_12_20')
    end
%     if ismac
%         load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCChR2Awake_UPDATED_06_12_21')
%     else
%         load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCChR2Awake_UPDATED_06_12_21')
%     end
    %%%% 6 odors
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_27Neurons.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_27Neurons.mat')
    end
    %     if ismac
    %         load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_panelOdors.mat')
    %     else
    %         load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_panelOdors.mat')
    %     end
    if strcmpi(singleORpanelOdors, 'single')
        %         data = dataStructPCX_ChR2_singleOdor;
        % pool anesth and awake data.
        data = [dataStructPCChR2Awake dataStructPCX_ChR2_singleOdor];
        
        %%%------ indexing awake and anesth experiments ----%%%
        % 1 = 'awake'; 0 = 'anesth'
        mouseStateIndex = [ones(1,size(dataStructPCChR2Awake,2)) zeros(1,size(dataStructPCX_ChR2_singleOdor,2))];
        for stateInd = 1:size(data,2)
            data(stateInd).mouseStateIndex = mouseStateIndex(stateInd);
        end
        %-------------------------------------------------------
        trialNum = 20;%%% max trial number
        
    elseif strcmpi(singleORpanelOdors, 'panel')
        data = dataStructPCX_ChR2;
        trialNum = 7;
    end
    
    titStr = 'chr2';
    
elseif strcmpi(virusType, 'arch')
    
    %%%%%%%%%%% ArchT %%%%%%%%
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataPCarch_6OdorsAllNeurons.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataPCarch_6OdorsAllNeurons.mat')
    end
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ArchT_singleOdor.mat')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ArchT_singleOdor.mat')
    end
    if strcmpi(singleORpanelOdors, 'single')
        data = dataStructPCX_ArchT_singleOdor;
        trialNum = 20;
    elseif strcmpi(singleORpanelOdors, 'panel')
        data = dataPCarch_6OdorsAllNeurons;
        data(7:18) = [];%no valid respiration
        trialNum = 7;
    end
    titStr = 'ArchT';
    
elseif strcmpi(virusType, 'Tbet-NpHR')
    %%%%%%%%%%%%%%%%%%% Tbet-NpHR %%%%%%%%%%%
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ArchInMT','dataStructPCX_ArchInMT')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ArchInMT','dataStructPCX_ArchInMT')
    end
    data = dataStructPCX_ArchInMT;
    titStr = 'Tbet-NpHR';
    trialNum = 20;
    
elseif strcmpi(virusType, 'JGC')
    %%%%%%%%%%%%%%%%%%% JGC %%%%%%%%%%%
    if ismac
        load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPC_JG_ChR2','dataStructPC_JG_ChR2')
    else
        load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPC_JG_ChR2','dataStructPC_JG_ChR2')
    end
    data = dataStructPC_JG_ChR2;
    titStr = 'JGC';
    trialNum = 20;
end
%% compute
odor = [];light = [];inhLoc = [];
fgLightAll = [];fgOdorAll = [];
finalStateIndex = [];
numOdors = 6;
allExc = [];allInh = [];
c = 0;odorLightDiff = [];
diffSignIndx = [];bgLightAll = [];bgOdorAll = [];
inhibitionIndx = [];excitationIndx = [];
thAlpha = .05;%/(numOdors);% correct for mult. comparisons
for i = 1:size(data,2)
    maxTrialNum = min(min(length(data(i).odorOnset(2).BgFg(:,2)) ...
        ,length(data(i).odorOnset(1).BgFg(:,2))), trialNum);
    
    %odor condition
    fgOdor = data(i).odorOnset(2).BgFg(1:maxTrialNum,2);
    bgOdor = data(i).odorOnset(2).BgFg(1:maxTrialNum,1);
    %odor+light condition
    fgLight = data(i).odorOnset(1).BgFg(1:maxTrialNum,2);
    bgLight = data(i).odorOnset(1).BgFg(1:maxTrialNum,1);
    % subtract baseline
    evOdor = fgOdor-mean(bgOdor);
    evLight = fgLight-mean(bgLight);
    % avoid means with different sign between conditions.
    if (mean(evOdor) > 0 & mean(evLight) < 0) | (mean(evOdor) < 0 & mean(evLight) > 0)
        diffSignIndx(end+1) = i;
                continue
    end
    
    [~,pTH1] = ttest(fgLight, bgLight);% excitatory response
    [~,pTH2] = ttest(fgOdor, bgOdor);% excitatory response
    if thType==0
        if (pTH1 < thAlpha) | ...
                (pTH2 < thAlpha)
            odor(end+1) = mean(fgOdor-mean(bgOdor));
            light(end+1) = mean(fgLight-mean(bgLight));
            bgLightAll(end+1) = mean(bgLight);
            bgOdorAll(end+1) = mean(bgOdor);
            inhLoc(end+1) = i;% find the location.
            c = c+1;
            % separate stat over exc and inh.
            if mean(evOdor) >= 0 & mean(evLight) >= 0
                allExc(end+1,:) = [mean(evOdor) mean(evLight)];
                excitationIndx(end+1) = i;
            elseif mean(evOdor) < 0 & mean(evLight) < 0
                allInh(end+1,:) = [mean(evOdor) mean(evLight)];
                inhibitionIndx(end+1) = i;
            end
        end
    elseif thType==1
        if (pTH2 < thAlpha)
            odor(end+1) = mean(fgOdor-mean(bgOdor));
            light(end+1) = mean(fgLight-mean(bgLight));
            bgLightAll(end+1) = mean(bgLight);
            bgOdorAll(end+1) = mean(bgOdor);
            inhLoc(end+1) = i;% find the location.
            c = c+1;
            if mean(evOdor) >= 0
                allExc(end+1,:) = [mean(evOdor) mean(evLight)];
            elseif mean(evOdor) < 0
                allInh(end+1,:) = [mean(evOdor) mean(evLight)];
            end
        end
    elseif thType==2
        if isfield(data, 'mouseStateIndex')
            finalStateIndex(end+1) = data(i).mouseStateIndex;
        end
        odor(end+1) = mean(fgOdor-mean(bgOdor));
        light(end+1) = mean(fgLight-mean(bgLight));
        bgLightAll(end+1) = mean(bgLight);
        bgOdorAll(end+1) = mean(bgOdor);
        inhLoc(end+1) = i;% find the location.
        c = c+1;
%                 if mean(evOdor) >= 0
%                     allExc(end+1,:) = [mean(evOdor) mean(evLight)];
%                     excitationIndx(end+1) = i;
%                 elseif mean(evOdor) < 0
%                     allInh(end+1,:) = [mean(evOdor) mean(evLight)];
%                     inhibitionIndx(end+1) = i;
%                 end
        if mean(evOdor) >= 0 & mean(evLight) >= 0
            allExc(end+1,:) = [mean(evOdor) mean(evLight)];
            excitationIndx(end+1) = i;
        elseif mean(evOdor) < 0 & mean(evLight) < 0
            allInh(end+1,:) = [mean(evOdor) mean(evLight)];
            inhibitionIndx(end+1) = i;
        end
    end
    
    
    %%%%%%%% check for signif. difference between conditions %%%%%%%%%%
    [~,p] = ttest(evOdor,evLight);
    if p < 0.05
        odorLightDiff(end+1) = c; % cases with significant difference between the condition`s fr.
    end
end
%% plot
%%%%% baseline analysis
% [~,pBG] = ttest(bgLightAll, bgOdorAll);
% disp(['p-value baseline: ' num2str(pBG)])
% meanBar = mean([bgLightAll' bgOdorAll']);
% seBar = std([bgLightAll' bgOdorAll'])./sqrt(length(bgOdorAll));
% figure;
% barwitherr(seBar, meanBar)
% set(gca, 'XTickLabel', {'Odor+Light', 'Odor'})
% ylabel('Evoked firing rate(spikes/sec)')
% set(gca,'fontSize',14)
% box off;
% set(gca,'linewidth',0.25)
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);

figure;
% plot the COP with significant difference between conditions.
% scatter(odor(odorLightDiff),light(odorLightDiff),30,'g','MarkerFaceColor','g','MarkerEdgeColor','none','LineWidth',4);%%cases with significant difference between conditions
% lsline
hold all;
if ~isempty(finalStateIndex)
    %Plot anesth data
    scatter(odor(find(finalStateIndex==0)), light (find(finalStateIndex==0)) , 20,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
    hold on;
    %Plot awake data
    scatter(odor(find(finalStateIndex==1)), light (find(finalStateIndex==1)) , 20,'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',4);
else% in case of only anesth recordings
    scatter(odor, light,30,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
end
lineLim = [round(min(min(odor),min(light))-2) round(max(max(odor),max(light))+2)];
xlim(lineLim);ylim(lineLim)
plot(lineLim,lineLim,'--k')
plot([0 0], lineLim, 'k--')
plot( lineLim,[0 0], 'k--')
xlabel('Evoked odor(Spikes/Sec)','fontSize',14);ylabel('Evoked odor+light(Spikes/Sec)','fontSize',14)
set(gca,'fontSize',14)
box off;
hold off;
set(gca,'linewidth',0.25)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
title(titStr)
hold off

%% stat
%exc
[~,pExc] = ttest(allExc(:,1),allExc(:,2));
disp(['p-value excitatory COPs ' num2str(pExc) ', N = ' num2str(size(allExc,1)) ' out of ' ...
    num2str(length(odor))])
%inh

[~,pInh] = ttest(allInh(:,1),allInh(:,2));
disp(['p-value inhibitory COPs ' num2str(pInh) ', N = ' num2str(size(allInh,1)) ' out of ' ...
    num2str(length(odor))])

%% Plot the awake data bar graph
figure;
deltaAwake = light(find(finalStateIndex==1)) - odor(find(finalStateIndex==1));
barwitherr(std(deltaAwake)./sqrt(length(deltaAwake)), mean(deltaAwake));
hold on;
 a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(deltaAwake)) + a;
scatter(r.*ones(1, length(deltaAwake)), deltaAwake , 40,'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',4);

ylabel('\Delta firing rate (light on - light off)','fontSize',14)
set(gca,'fontSize',14)
box off;
hold off;
set(gca,'linewidth',0.25)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);

[~,deltaAwakePvalue] = ttest(deltaAwake);
title(['P = ' num2str(deltaAwakePvalue)])

%% bar graphs
dataInh = data(inhibitionIndx);
dataExc = data(excitationIndx);

barPCchr2ExcInh(dataInh, dataExc, singleORpanelOdors)
%% Percent change in response
percentChnagePCresponse(dataInh, dataExc)

%%% percent change of PC odor response comparison - Awake vs Aneth %%%
if strcmpi(virusType, 'chr2') & strcmpi(singleORpanelOdors, 'single')
    percentChangeFRpcAwakeAnesth(dataExc)
end

%% Population mean PSTH
if strcmpi(virusType, 'Tbet-NpHR') | strcmpi(virusType, 'JGC')
%     populationPSTH([dataExc dataInh], virusType, 'inhalationOnset', [-.5 2])
      populationPSTH(dataExc, virusType, 'inhalationOnset', [-.5 2])
else
    populationPSTH(dataExc, virusType, 'inhalationOnset', [-.5 2])
end

%% Response latency
PCresponseLatency(dataExc, virusType)


%% params
clear;clc
sCalcTime = [0 2];
odorVec  = [1:4 6 7];
maxTrialNum = 5;
virType = [0, 1]; % ArchT  and chr2.
meanWithin = []; seWithin = [];
allOdorsMeanWithin = [];allOdorsMeanWithinSE = [];allValuesMeanWithin = [];
ChR2ArchTWithinCorr = [];
%%
for jj = 1:length(virType)
    %% load data
    %%%%%%%%%%%%%%%%% PC - GC chr2 %%%%%%%%%%%%%%%%%%%%%%%
    if virType(jj)==1
        
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_27Neurons.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_27Neurons.mat')
        end
        data = dataStructPCX_ChR2;
       
        % if ismac
        %     load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_panelOdors.mat')
        % else
        %     load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_panelOdors.mat')
        % end
        
        
    elseif virType(jj)==0
        %%%%%%%%%%%%%%%%% PC - GC Arch %%%%%%%%%%%%%%%%%%%%%%%
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataPCarch_6OdorsAllNeurons.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataPCarch_6OdorsAllNeurons.mat')
        end
        data = dataPCarch_6OdorsAllNeurons;
       
    end
    
    %% compute correlation
    labels = [];dbOL = [];dbOdor = [];
    
    %% odor+light database
    dbOL = [];
    for i = 1:size(data,2)% run over all cell-odor pairs
        if ~isempty(data(i).InhOnset)
            fg =[];bg = [];
            for j = 1:min(maxTrialNum, ...
                    min(size(data(i).InhOnset(1).rasterPlot,2), size(data(i).InhOnset(2).rasterPlot,2)))
                % run over trials.
                spikes = data(i).InhOnset(1).rasterPlot(j).times;
                
                bg(j) = length(find(spikes >= -sCalcTime(2) ...
                    & spikes < sCalcTime(1)))/sCalcTime(2);% baseline
                
                fg(j) = length(find(spikes >= sCalcTime(1) ...
                    & spikes < sCalcTime(2)))/sCalcTime(2); % forground
                
            end
            labels(end+1) = data(i).odorNum;
            % save the evoked response
            dbOL(end+1,:) = fg-mean(bg);
        end
    end
    %% odor only database
    fg = [];bg = [];dbOdor = [];labels = [];
    for i = 1:size(data,2)% run over all cell-odor pairs
        if ~isempty(data(i).InhOnset)
            fg = [];bg = [];
            for j = 1:min(maxTrialNum, ...
                    min(size(data(i).InhOnset(1).rasterPlot,2), size(data(i).InhOnset(2).rasterPlot,2)))
                spikes = data(i).InhOnset(2).rasterPlot(j).times;% run over trials.
                
                bg(j) = length(find(spikes >= -sCalcTime(2) ...
                    & spikes < sCalcTime(1)))/sCalcTime(2);% baseline
                
                fg(j) = length(find(spikes >= sCalcTime(1) ...
                    & spikes < sCalcTime(2)))/sCalcTime(2); % forground
            end
            labels(end+1) = data(i).odorNum;
            % save the evoked response
            dbOdor(end+1,:) = fg-mean(bg);
        end
    end
    %% build the odor-neuron matrix
    odorLightNeuronMat = [];
    odorNeuronMat = [];
    for i = 1:length(odorVec)
        loc = find(labels==odorVec(i));
        for j = 1:maxTrialNum
            odorLightNeuronMat(end+1, :) = dbOL(loc,j)';
            odorNeuronMat(end+1, :) = dbOdor(loc, j)';
        end
    end
    %% plot
    %%%%%%%%%%%%%% odor+light %%%%%%%%%%%%%%%%
    % figure;
    Rlight = corr(odorLightNeuronMat');
    %     imagesc(Rlight);
    % title('odor+chr2')
    % caxis([.2 1])
    Rbetween = Rlight;
    withCorr = [];meanWithCorr = [];
    seWithCorr = [];
    for i = 1:maxTrialNum:size(Rlight,1)
        res = triu(Rlight(i:i+maxTrialNum-1,i:i+maxTrialNum-1),1);
        Rbetween(i:i+maxTrialNum-1,i:i+maxTrialNum-1) = 0;
        upperTriangleCorr = res(find(res~=0))';
        withCorr = [withCorr upperTriangleCorr];
        meanWithCorr(end+1) = mean(upperTriangleCorr);
        seWithCorr(end+1) = std(upperTriangleCorr)./sqrt(length(upperTriangleCorr));
        upperTriangleCorr = [];
    end
    
    BetweenCorr = triu(Rbetween,1);
    BetweenCorr = BetweenCorr(find(BetweenCorr~=0));
    
    allWithinCorrLight = withCorr;
    allMeanWithinCorrLight = meanWithCorr;
    allSEWithinCorrLight = seWithCorr;
    allBetweenCorrLight = BetweenCorr';
    %
    % %%%%%%%%%%%%%% odor only %%%%%%%%%%%%%%%%
    % figure;
    Rodor = corr(odorNeuronMat');
    %     imagesc(Rodor);
    % caxis([0.1 1])
    % title('odor')
    Rbetween = Rodor;
    withCorr = [];meanWithCorr = [];
    seWithCorr = [];
    for i = 1:maxTrialNum:size(Rodor,1)
        res = triu(Rodor(i:i+maxTrialNum-1,i:i+maxTrialNum-1),1);
        Rbetween(i:i+maxTrialNum-1,i:i+maxTrialNum-1) = 0;
        upperTriangleCorr = res(find(res~=0))';
        withCorr = [withCorr upperTriangleCorr];
        meanWithCorr(end+1) = mean(upperTriangleCorr);
        seWithCorr(end+1) = std(upperTriangleCorr)./sqrt(length(upperTriangleCorr));
        upperTriangleCorr = [];
    end
    
    BetweenCorr = triu(Rbetween,1);
    BetweenCorr = BetweenCorr(find(BetweenCorr~=0));
    
    allWithinCorrOdor = withCorr;
    allMeanWithinCorrOdor = meanWithCorr;
    allSEWithinCorrOdor = seWithCorr;
    allBetweenCorrOdor = BetweenCorr';
    
    ChR2ArchTWithinCorr = [ChR2ArchTWithinCorr (allWithinCorrLight-allWithinCorrOdor)'];
    %% plot
    %%% errorbars
    %------within using all 5 choose 2 values
    meanWithin(1,end+1) = mean(allWithinCorrLight)/mean(allWithinCorrOdor);
    seWithin(1,end+1) = (std(allWithinCorrLight)/std(allWithinCorrOdor))./sqrt(length(allWithinCorrOdor));
%     %stat
%     [~,pWithin] = ttest(allWithinCorrLight, allWithinCorrOdor);
    
    % mean within values using the mean from each odor
    allOdorsMeanWithin(1,end+1) = (mean((allMeanWithinCorrLight)./(allMeanWithinCorrOdor))).*100;
    allOdorsMeanWithinSE(1,end+1) = (std((allMeanWithinCorrLight)./(allMeanWithinCorrOdor))./sqrt(length(allMeanWithinCorrLight))).*100;
    
    % the relative values
    allValuesMeanWithin(end+1,:) = ((allMeanWithinCorrLight)./(allMeanWithinCorrOdor)).*100;
    %stat
    [~,pWithin] = ttest(allMeanWithinCorrLight, allMeanWithinCorrOdor);
    
    if virType(jj)==1
        disp(['N = ' num2str(length(allMeanWithinCorrLight)) ', Within ChR2, p = ' num2str(pWithin)])
    else
        disp(['N = ' num2str(length(allMeanWithinCorrLight)) ', Within ArchT, p = ' num2str(pWithin)])
    end
end
%========================= PLOT ======================================
% Using all within values 6* 5choose2
figure;
barwitherr(seWithin, meanWithin)% plot [left: ArchT, right: ChR2]
hold all;
plot([.5 2.5], [1 1], 'k--')
box off;
ylabel('Correlation (% of odor only)')
set(gca,'fontSize',14)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
ylim([.5 1.5])
hold off;
title('All 6 * 5choose2 values')

figure;
%for jitter
a = 0.9;
b = 1.1;
r = (b-a).*rand(1,length(allMeanWithinCorrOdor)) + a;

barwitherr(allOdorsMeanWithinSE, allOdorsMeanWithin)
hold all;
scatter(r.*ones(1, length(allMeanWithinCorrOdor)), allValuesMeanWithin(1,:));%ArchT
scatter(r.*2.*ones(1, length(allMeanWithinCorrOdor)), allValuesMeanWithin(2,:));%+ChR2
plot([.5 2.5], [100 100], 'k--')
box off;
ylabel('Change in mean correlation (%)')
set(gca,'fontSize',14)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
hold off;
xticklabels({'+ArchT', '+ChR2'});xtickangle(45)
title('Mean values per odor')

%------------------- PLOT THE DISTRIBUTION OF DELTA VALUES ------------
figure;
deltaCorrBins = [-.8:.1:.8];%delta corr
%+ArchT
histogram(ChR2ArchTWithinCorr(:,1), deltaCorrBins, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'r')
hold all;

%+ChR2
histogram(ChR2ArchTWithinCorr(:,2), deltaCorrBins,  'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'b')
legend('+ArchT', '+Chr2')

%Plot the means+SEM
ciWithinCorrArchT = std(ChR2ArchTWithinCorr(:,1))./sqrt(length(ChR2ArchTWithinCorr(:,1)))*tinv(.975, size(ChR2ArchTWithinCorr,1)-1);
ciWithinCorrChR2 = (std(ChR2ArchTWithinCorr(:,2))./sqrt(length(ChR2ArchTWithinCorr(:,2))))*tinv(.975, size(ChR2ArchTWithinCorr,1)-1);

errorbar(mean(ChR2ArchTWithinCorr(:,1)),.25,ciWithinCorrArchT ,'horizontal', 'ro')
errorbar(mean(ChR2ArchTWithinCorr(:,2)),.25, ciWithinCorrChR2 ,'horizontal', 'bo')

[~,pArchWithinAllValues] = ttest(ChR2ArchTWithinCorr(:,1));
[~,pChR2WithinAllValues] = ttest(ChR2ArchTWithinCorr(:,2));
xlabel('\Delta correlation within trials (light on - light off)');ylabel('Proportion')
hold off;
box off;
title(['Mean + 95% CI: ArchT: P=' num2str(pArchWithinAllValues) ', ChR2: P=' num2str(pChR2WithinAllValues)])

%------------------- PLOT THE DISTRIBUTION OF DELTA VALUES ------------
meanDeltaCorrPerOdor = [mean(reshape(ChR2ArchTWithinCorr(:,1)',length(ChR2ArchTWithinCorr(:,1))./length(odorVec), length(odorVec))); ... %ArchT
                       mean(reshape(ChR2ArchTWithinCorr(:,2)',length(ChR2ArchTWithinCorr(:,2))./length(odorVec), length(odorVec)))];%ChR2
                   
seDeltaCorrPerOdor = [std(reshape(ChR2ArchTWithinCorr(:,1)',length(ChR2ArchTWithinCorr(:,1))./length(odorVec), length(odorVec)))./sqrt(length(odorVec)); ... 
                       std(reshape(ChR2ArchTWithinCorr(:,2)', length(ChR2ArchTWithinCorr(:,1))./length(odorVec), length(odorVec)))./sqrt(length(odorVec))];
     
[meanDeltaCorrPerOdorSORTED,meanDeltaCorrPerOdorINDEX] = sort(meanDeltaCorrPerOdor', 'descend');
plotBins = .75:.1:1.25;% six bins for 6 odors.
                   
figure;
% errorbar(1+plotBins, meanDeltaCorrPerOdor(1,:), seDeltaCorrPerOdor(1,:), 'or')%ArchT
% hold on;
% errorbar(plotBins, meanDeltaCorrPerOdor(2,:), seDeltaCorrPerOdor(2,:), 'ob')%ChR2
errorbar(1+plotBins, meanDeltaCorrPerOdorSORTED(:,1), seDeltaCorrPerOdor(1, meanDeltaCorrPerOdorINDEX(:,1)), 'or', 'LineWidth', .25 , 'MarkerFaceColor', 'r')%ArchT
hold on;
errorbar(plotBins, meanDeltaCorrPerOdorSORTED(:, 2), seDeltaCorrPerOdor(2,meanDeltaCorrPerOdorINDEX(:,2)), 'ob', 'LineWidth', .25, 'MarkerFaceColor', 'b')%ChR2

box off;
ylabel({'\Delta correlation within trials' ; '(light on - light off)'})
set(gca,'fontSize',14)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'linewidth',0.25)
xlim([.5 2.5])
hold off;

           


% %between
% meanWithin = mean([allBetweenCorrOdor./allBetweenCorrLight]);
% seWithin = std([allBetweenCorrOdor./allBetweenCorrLight])./sqrt(length(allBetweenCorrOdor));
% hold on;
% errorbar(1:2, meanWithin, seWithin, 'g')
% xlim([0 3])
% box off;
% set(gca,'fontSize',14)
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% set(gca,'linewidth',0.25)
% %stat
% [~,pBetween] = ttest(allBetweenCorrLight, allBetweenCorrOdor);
% title(['Between, p = ' num2str(pBetween)])
% ylim([.2 .8])

% plot within\between corr distribution
% within
% k = -1:.1:1;
% figure;
% [a,b] = hist(allWithinCorrOdor,k);
% plot(b,a./sum(a),'k')
% hold on;
% [a,b] = hist(allWithinCorrLight,k);
% plot(b,a./sum(a),'b')
% %stat
% [~,pWithin] = kstest2(allWithinCorrLight, allWithinCorrOdor);
% title(['N = ' num2str(length(allWithinCorrLight)) ' pairs, Within correlation, KS: ' num2str(pWithin)])
% % legend('Odor', 'Odor+ChR2')
% xlabel('Correlation');ylabel('% of pairs')
% box off;
% set(gca,'fontSize',14)
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% set(gca,'linewidth',0.25)
% %between
% figure;
% [a,b] = hist(allBetweenCorrOdor,k);
% plot(b,a./sum(a),'k')
% hold on;
% [a,b] = hist(allBetweenCorrLight,k);
% plot(b,a./sum(a),'b')
% %stat
% [~,pBetween] = kstest2(allBetweenCorrOdor, allBetweenCorrLight);
% title(['N = ' num2str(length(allBetweenCorrOdor)) ' pairs, Between correlation, KS: ' num2str(pBetween)])
% % legend('Odor', 'Odor+ChR2')
% xlabel('Correlation');ylabel('% of pairs')
% box off;
% set(gca,'fontSize',14)
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% set(gca,'linewidth',0.25)


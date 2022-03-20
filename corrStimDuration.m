%%%  Description
% This function copute the correlation between PC odor responses under
% expression of ChR2 or ArchT in the GCL, for odor/odor+light conditions.
%The computation of spikes is during the odor presentation duration.
%% similarity analysis
%% params
clear;clc
sCalcTime = [0 2];
odorVec  = [1:4 6 7];
inh=1;virType=[1 0];%1=chr2, 0 = arch.
allMean = [];allSEs = [];
allVals = [];
maxTrialNum = 5;
c = 0;

%%
for jj = 1:length(virType)% run over chr2 and arch data
    %% load data
    %%%%%%%%%%%%%%%%% PC - GC chr2 %%%%%%%%%%%%%%%%%%%%%%%
    if virType(jj)==1

        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_27Neurons.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_27Neurons.mat')
        end
%                 if ismac
%                     load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_panelOdors.mat')
%                 else
%                     load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_panelOdors.mat')
%                 end
        
        data = dataStructPCX_ChR2;
    elseif virType(jj)==0
        %%%%%%%%%%%%%%%%% PC - GC Arch %%%%%%%%%%%%%%%%%%%%%%%
        %         if ismac
        %             load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataPCXArchdecoder_210420.mat')
        %         else
        %             load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataPCXArchdecoder_210420.mat')
        %         end
        %         data = dataPCXArchdecoder_210420;
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
            dbOL(end+1) = mean(fg-mean(bg));
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
            dbOdor(end+1) = mean(fg-mean(bg));
        end
    end
    
    %% build the odor-neuron matrix
    odorLightNeuronMat = [];
    odorNeuronMat = [];
    for i = 1:length(odorVec)
        loc = find(labels==odorVec(i));
        odorLightNeuronMat(:,i) = [dbOL(loc)];
        odorNeuronMat(:,i) = [dbOdor(loc)];
    end
    
    %% Compute population sparseness
    sparseOdor = [];sparseOdorLight = [];
    for sparseInd = 1:size(odorNeuronMat,1)
        sparseOdor(sparseInd) = (1-(sum(odorNeuronMat(sparseInd,:)./size(odorNeuronMat,2))).^2./sum(odorNeuronMat(sparseInd,:).^2./size(odorNeuronMat,2)))./(1-1/size(odorNeuronMat,2));
        sparseOdorLight(sparseInd) = (1-(sum(odorLightNeuronMat(sparseInd,:)./size(odorLightNeuronMat,2))).^2./sum(odorLightNeuronMat(sparseInd,:).^2./size(odorLightNeuronMat,2)))./(1-1/size(odorLightNeuronMat,2));
    end
    [~,pLifeTimeSparse] = ttest(sparseOdor, sparseOdorLight);
    deltaLTS = sparseOdorLight-sparseOdor;
    allDeltaLTSValues{jj} = deltaLTS;
    allLifeTimeSparse(jj, :) = [std(deltaLTS)./sqrt(length(deltaLTS)), mean(deltaLTS) pLifeTimeSparse];
    
    %% Compute lifetime sparseness
    spOdor = [];spOdorLight = [];
    spOdorMat = odorNeuronMat';
    spOdorLightMat = odorLightNeuronMat';
    
    for sparseInd = 1:size(spOdorMat,1)
        spOdor(sparseInd) = (1-(sum(spOdorMat(sparseInd,:)./size(spOdorMat,2))).^2./sum(spOdorMat(sparseInd,:).^2./size(spOdorMat,2)))./(1-1/size(spOdorMat,2));
        spOdorLight(sparseInd) = (1-(sum(spOdorLightMat(sparseInd,:)./size(spOdorLightMat,2))).^2./sum(spOdorLightMat(sparseInd,:).^2./size(spOdorLightMat,2)))./(1-1/size(spOdorLightMat,2));
    end
    [~,pSP] = ttest(spOdor, spOdorLight);
    deltaSP = spOdorLight-spOdor;
    allDeltaSPValues{jj} = deltaSP;
    allPopulationSparseness(jj, :) = [std(deltaSP)./sqrt(length(deltaSP)), mean(deltaSP) pSP];
    %% compute the distance
    rLight = pdist(odorLightNeuronMat', 'correlation');
    rLight = 1-rLight;
    
    rOdor = pdist(odorNeuronMat', 'correlation');
    rOdor = 1-rOdor;
    c = c+1;
    %     figure;
    subplot(1,4,c)
    R = corr(odorNeuronMat);
    imagesc(R);colormap jet;title('Odor')
    caxis([0 1])
    c = c+1;
    subplot(1,4,c)
    R = corr(odorLightNeuronMat);
    imagesc(R);colormap jet;title('Odor+Light')
    caxis([0.4 1])
    %% stat
    [~,p] = ttest(rOdor,rLight);
    disp(['p = ' num2str(p)])
    
    %% mean+-SEM
    meanRodor = mean(rOdor);
    seRodor = std(rOdor)./sqrt(length(rOdor));
    
    meanRlight = mean(rLight);
    seRLight = std(rLight)./sqrt(length(rLight));
    
    allMean = [allMean meanRodor meanRlight];
    allSEs = [allSEs seRodor seRLight];
    allVals = [allVals rOdor' rLight'];
end
%%
%COMPUTE RATIOS
chr2CorrValues = (allVals(:,2)./allVals(:,1)).*100;
archCorrValues =  (allVals(:,4)./allVals(:,3)).*100;
%COMPUTE DELTAS
chr2CorrDelta = (allVals(:,2)-allVals(:,1));
archCorrDelta =  (allVals(:,4)-allVals(:,3));

%--- MEAN+-SEM
%Ratios
chr2Corr = mean(chr2CorrValues);
archCorr =  mean(archCorrValues);
seCorr = [std(chr2CorrValues)./sqrt(size(allVals,1)) ...
    std(archCorrValues)./sqrt(size(allVals,1))];

%Deltas
chr2CorrDeltaMean = mean(chr2CorrDelta);
archCorrDeltaMean =  mean(archCorrDelta);
seCorrDelta = [std(chr2CorrDelta)./sqrt(size(allVals,1)) ...
    std(archCorrDelta)./sqrt(size(allVals,1))];

%% plot both conditions


% Plot actual values
figure;
barwitherr(allSEs, allMean)
ylabel('Correlation coefficient')
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xtickangle(45)
set(gca, 'XTickLabel', {'Odor', 'Odor+ChR2', 'Odor', 'Odor+Arch'})

% plot the bar graph of correlATION RATIOS
figure;
barwitherr(seCorr, [chr2Corr archCorr])
hold all;
plot([0 3], [100 100], 'k--');
% scatter the data points
scatter(ones(1,length(chr2CorrValues)), chr2CorrValues,20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
scatter(2.*ones(1,length(chr2CorrValues)), archCorrValues,20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
ylabel('Correlation (% from odor only)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xtickangle(45)
set(gca, 'XTickLabel', {'Odor+ChR2' 'Odor+Arch'})
ylim([minmax([chr2CorrValues archCorrValues])])
hold off;

%for jitter
a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(chr2CorrDelta)) + a;
% plot the bar graph of correlATION Deltas
figure;
barwitherr(seCorrDelta, [chr2CorrDeltaMean archCorrDeltaMean])
hold all;
% scatter the data points
scatter(r.*ones(1,length(chr2CorrDelta)), chr2CorrDelta,20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
scatter(1+r.*ones(1,length(archCorrDelta)), archCorrDelta,20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
ylabel({'\Delta correlation coefficient (light on - light off)'})
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xtickangle(45)
set(gca, 'XTickLabel', {'+ChR2' '+Arch'})
hold off;


%%
figure;
%+Chr2
[sortROdorChR2, locCHR2] = sort(allVals(:,1));
sortROdorLightChR2 = allVals(locCHR2, 2);
plot(sortROdorLightChR2./sortROdorChR2, 'bs-')
hold all
%+ArchT
[sortROdorArchT, locArchT] = sort(allVals(:,3));
sortROdorLightArchT = allVals(locArchT, 4);
plot(sortROdorLightArchT./sortROdorArchT, 'rs-')
plot([1 size(allVals, 1)], [1 1], 'k--')
hold off;
box off;
ylabel('Sorted correlation ratio (%)','fontSize',14)
xlabel('Correlation pairs','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);

%% Lifetime sparsness
figure;
barwitherr([allLifeTimeSparse(1,1) allLifeTimeSparse(2,1)], [allLifeTimeSparse(1,2) allLifeTimeSparse(2,2)])
%Plot the data points in a random scatter
a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(allDeltaLTSValues{1})) + a;
hold all;
%ChR2
scatter(r.*ones(1,length(allDeltaLTSValues{1})), allDeltaLTSValues{1},20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
%ArchT
r = (b-a).*rand(1,length(allDeltaLTSValues{2})) + a;
scatter(1+r.*ones(1,length(allDeltaLTSValues{2})), allDeltaLTSValues{2},20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
%Stat
title(['ChR2: P = ' num2str(allLifeTimeSparse(1,3)) ', ArchT: P = ' num2str(allLifeTimeSparse(2,3))])
ylabel('\Delta Lifetime sparseness','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xtickangle(45)
set(gca, 'XTickLabel', {'+ChR2' '+Arch'})

%% populatione sparsness
figure;
barwitherr([allPopulationSparseness(1,1) allPopulationSparseness(2,1)], [allPopulationSparseness(1,2) allPopulationSparseness(2,2)])
%Plot the data points in a random scatter
a = 0.85;
b = 1.15;
r = (b-a).*rand(1,length(allDeltaSPValues{1})) + a;
hold all;
%ChR2
scatter(r.*ones(1,length(allDeltaSPValues{1})), allDeltaSPValues{1},20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
%ArchT
scatter(1+r.*ones(1,length(allDeltaSPValues{2})), allDeltaSPValues{2},20,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none','LineWidth',4);
%Stat
title(['ChR2: P = ' num2str(allPopulationSparseness(1,3)) ', ArchT: P = ' num2str(allPopulationSparseness(2,3))])
ylabel('\Delta population sparseness','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xtickangle(45)
set(gca, 'XTickLabel', {'+ChR2' '+Arch'})
hold off;
%% example of ArchT correlation between two odors - how the correlation increases.
% Figure
% ind = [3,4];%[1 3];
% figure;
% subplot(1,2,1)
% scatter(odorLightNeuronMat(:,ind(1)), odorLightNeuronMat(:,ind(2)), 'r')
% hold on;
% scatter(mean(odorLightNeuronMat(:,ind(1))), mean(odorLightNeuronMat(:,ind(2))), 'r*')
% [r,~] = corr(odorLightNeuronMat(:,ind(1)), odorLightNeuronMat(:,ind(2)));
% xlabel(['Odor #' num2str(ind(1))]);ylabel(['Odor #' num2str(ind(2))]);
% title(['Odor+Light, R = ' num2str(r)])
% 
% 
% subplot(1,2,2)
% scatter(odorNeuronMat(:,ind(1)), odorNeuronMat(:,ind(2)), 'k')
% hold on;
% scatter(mean(odorNeuronMat(:,ind(1))), mean(odorNeuronMat(:,ind(2))), 'k*')
% [r,~] = corr(odorNeuronMat(:,ind(1)), odorNeuronMat(:,ind(2)));
% xlabel(['Odor #' num2str(ind(1))]);ylabel(['Odor #' num2str(ind(2))]);
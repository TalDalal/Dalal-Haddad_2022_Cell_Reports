virusType = {'chr2', 'arch', 'Tbet-ArchT'};allVals = [];thAlpha = 0.05;
for jj = 1:length(virusType)
    %%%% description:
    % This function computes the Cv for each cell odor pair across trials.
    %%%%%%%%%%%%%%%%% ChR2
    if strcmpi(virusType{jj}, 'chr2')
        % Single odors
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_singleOdor.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_singleOdor.mat')
        end
        %%%% 6 odors
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ChR2_27Neurons.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ChR2_27Neurons.mat')
        end
        
        % awake-head fixed cell-odor pairs
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/PC_ChR2_AWAKE_FIRING_RATE_20_12_20')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\PC_ChR2_AWAKE_FIRING_RATE_20_12_20')
        end
        % combine datasets
        data = [dataStructPCX_ChR2_singleOdor dataStructPCX_ChR2 dataStructPCChR2Awake];
        xlabelTitle = {'Odor', '+ChR2'};
        
    elseif strcmpi(virusType{jj}, 'arch')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% ArchT %%%%%%%%%%%%%
        
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
        data = [dataPCarch_6OdorsAllNeurons dataStructPCX_ArchT_singleOdor];
        
        xlabelTitle = {'Odor', '+ArchT'};
    elseif strcmpi(virusType{jj}, 'Tbet-ArchT')
        %%%%%%%%%%%%%%%%%%% Tbet-NpHR %%%%%%%%%%%
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructPCX_ArchInMT','dataStructPCX_ArchInMT')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructPCX_ArchInMT','dataStructPCX_ArchInMT')
        end
        data = dataStructPCX_ArchInMT;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    stat = 'mean';cvOdor = [];cvLight = [];
    for i = 1:size(data,2)
        % odor +light
        fgLight = data(i).odorOnset(1).BgFg(:,2);
        bgLight = data(i).odorOnset(1).BgFg(:,1);
        %odor
        fgOdor = data(i).odorOnset(2).BgFg(:,2);
        bgOdor = data(i).odorOnset(2).BgFg(:,1);
        
        % compute evoked response.
        evLight = fgLight-mean(bgLight);
        evOdor = fgOdor-mean(bgOdor);
        % avoid means with different sign between conditions.
        if (mean(evOdor) > 0 & mean(evLight) < 0) | (mean(evOdor) < 0 & mean(evLight) > 0)
            continue
        end
        
        [~,pTH1] = ttest(fgLight, bgLight);% excitatory response
        [~,pTH2] = ttest(fgOdor, bgOdor);% excitatory response
        if (pTH1 < thAlpha) | (pTH2 < thAlpha)
            
            % compute evoked response.
            light = abs(fgLight-mean(bgLight));
            Odor = abs(fgOdor-mean(bgOdor));
            
            % cv
            cvLight(end+1) = std(light)./mean(light);
            cvOdor(end+1) = std(Odor)./mean(Odor);
        end
        
    end
    allVals{jj} = [cvOdor' cvLight'];
    [~,pTtest] = ttest(cvLight, cvOdor);
    disp([virusType{jj} ': #cell-odor pairs:  ' num2str(size(cvLight,2)) ', p = ' num2str(pTtest)]);
end


chr2CV = allVals{1};
archCV = allVals{2};
nphrCV = allVals{3};

%%%% plot both ChR2 and ArchT in actual values %%%
figure;
barwitherr([std(chr2CV)./sqrt(size(chr2CV,1)) std(archCV)./sqrt(size(archCV,1))], [mean(chr2CV) mean(archCV)])
hold all;
plot(chr2CV', 'Color', [.5 .5 .5], 'LineWidth', .1)
plot([3:4], archCV', 'Color', [.5 .5 .5], 'LineWidth', .1)
hold off;
ylabel('Odor response variability')
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'XTickLabel', {'Odor', '+ChR2', 'Odor', '+Arch'});xtickangle(45);
set(gca, 'XColor', 'k', 'YColor', 'k')

% MEAN+-SEM
% chr2CorrMean = mean(chr2CV);chr2CorrMean = chr2CorrMean(2)/chr2CorrMean(1);
% archCorrMean = mean(archCV);archCorrMean = archCorrMean(2)/archCorrMean(1);
% seCorr = [std(chr2CV(:,2)./chr2CV(:,1))./sqrt(size(chr2CV,1)) ...
%     std(archCV(:,2)./archCV(:,1))./sqrt(size(archCV,1))];
% %%%% plot both ChR2 and ArchT in relative plot %%%
% figure;
% barwitherr(seCorr*100, [chr2CorrMean archCorrMean].*100)
% ylabel('Coefficient of variation (% of odor only)','fontSize',14)
% set(gca,'fontSize',14)
% box off;
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
% set(gca, 'XTickLabel', {'Odor+ChR2' 'Odor+Arch'});xtickangle(45);
% ylim([80 120])
% hold on;
% plot([0 3], [100 100], 'k--')

%%%% plot Tbet-NpHR in another figure %%%
figure;
barwitherr(std(nphrCV)./sqrt(size(nphrCV,1)), mean(nphrCV))
ylabel('Coefficient of variation','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'XTickLabel', {'Odor' 'Tbet-NpHR'});xtickangle(45);
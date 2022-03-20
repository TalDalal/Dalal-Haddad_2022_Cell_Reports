%% Compute the odor response variability of MTC responses (Cv)
virusType = {'chr2', 'arch', 'Tbet-ArchT'};allVals = [];thAlpha = 0.05;
for jj = 1:length(virusType)
    %%%% description:
    % This function computes the Cv for each cell odor pair across trials.
    %%%%%%%%%%%%%%%%% ChR2
    if strcmpi(virusType{jj}, 'chr2')
        % For single electrode recordings
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOBChR2_12032020.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOBChR2_12032020.mat')
        end
        % For neuronexus recordings
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureChR2MTC_for_CV_analysis.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureChR2MTC_for_CV_analysis.mat')
        end
        dataStructOBChR2 = rmfield(dataStructOBChR2, {'BgFg', 'pValue', 'InhOnset', 'odorNum', 'stimDur', 'ISI', 'neuronNum'});
        %Combine datasets.
        data = [dataStructOBChR2 FRstructMTC];
        FRstructMTC = [];
        
        xlabelTitle = {'Odor', '+ChR2'};
        
    elseif strcmpi(virusType{jj}, 'arch')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% ArchT %%%%%%%%%%%%%
        
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
        
        % For LFP recordings
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureArchTMTC_for_CV_analysis.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureArchTMTC_for_CV_analysis.mat')
        end
        dataStructOB = rmfield(dataStructOB, {'BgFg', 'pValue', 'InhOnset', 'odorNum', 'stimDur', 'ISI', 'neuronNum'});
        dataStructOBInhibited = rmfield(dataStructOBInhibited, {'BgFg', 'pValue', 'InhOnset', 'odorNum', 'stimDur', 'ISI', 'neuronNum'});
        
        %Combine datasets.
        data = [dataStructOB dataStructOBInhibited FRstructMTC];% combine both data sets.
        FRstructMTC = [];
        
        xlabelTitle = {'Odor', '+ArchT'};
    elseif strcmpi(virusType{jj}, 'Tbet-ArchT')
        %%%%%%%%%%%%%%%%%%% Tbet-NpHR %%%%%%%%%%%
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/dataStructOB_MTarchT.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\dataStructOB_MTarchT.mat')
        end
        
        % For LFP recordings
        if ismac
            load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/FRstructureNpHRMTC_for_CV_analysis.mat')
        else
            load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\FRstructureNpHRMTC_for_CV_analysis.mat')
        end
        dataStructOB_MTarchT = rmfield(dataStructOB_MTarchT, {'expID', 'InhOnset', 'odorNum', 'stimDur', 'ISI', 'neuronNum'});
        %Combine datasets.
        data = [dataStructOB_MTarchT FRstructMTC];
        FRstructMTC = [];

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
    disp([virusType{jj} ': N = ' num2str(size(cvLight,2)) ' #cell-odor pairs , p = ' num2str(pTtest)]);
end


chr2CV = allVals{1};
archCV = allVals{2};
nphrCV = allVals{3};

%%%% plot both ChR2 and ArchT in actual values %%%
figure;
barwitherr([std(chr2CV)./sqrt(size(chr2CV,1)) std(archCV)./sqrt(size(archCV,1)) std(nphrCV)./sqrt(size(nphrCV,1))], ...
    [mean(chr2CV) mean(archCV) mean(nphrCV)])
hold all;
plot(chr2CV', 'Color', [.5 .5 .5], 'LineWidth', .1)
plot([3:4], archCV', 'Color', [.5 .5 .5], 'LineWidth', .1)
plot([5:6], nphrCV', 'Color', [.5 .5 .5], 'LineWidth', .1)
hold off;
ylabel('Odor response variability')
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca, 'XTickLabel', {'Odor', '+ChR2', 'Odor', '+Arch', 'Odor', 'Tbet-NpHR'});xtickangle(45);
set(gca, 'XColor', 'k', 'YColor', 'k')


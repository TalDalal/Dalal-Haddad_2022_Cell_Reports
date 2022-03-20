%% Description:
% This function plots the decoding success rate as a function of the #
% neurons, for both arch and chr2 experiemnts. we subtract the odor only
% success rate from the odor+light success rate and plot the delta SR.
%% run the decoder for both virus types
[matChR2_Light,matChR2_Odor, numOfNeurons_ChR2] = decodingAnalysisPCchr2('chr2');
[matArchT_Light,matArchT_Odor, numOfNeurons_Arch] = decodingAnalysisPCchr2('arch');

%% compute mean+-SEM
chr2Dec = matChR2_Light-matChR2_Odor;
mean_ChR2 = mean(chr2Dec);ste_ChR2 = std(chr2Dec);

archDec = matArchT_Light-matArchT_Odor;
mean_Arch = mean(archDec);ste_Arch = std(archDec);
%% plot
figure;
errorbar_patch(numOfNeurons_ChR2,mean_ChR2,ste_ChR2,'b')%% odor+light
hold all;
errorbar_patch(numOfNeurons_Arch,mean_Arch,ste_Arch,'y')%% odor+light
% zero marks the odor sucees rate.
chance = 0;
plot([0 numOfNeurons_Arch(end)],[chance chance],'k--')
box off;
xlabel('# Neurons', 'FontSize',14)
ylabel('Classification success rate (On - off)','FontSize',14);
box off;
set(gca,'linewidth',.25)
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
set(gca,'fontSize',14)
set(gcf,'Renderer', 'painter');


%========intensity analysis PC-OB chr2 ==================
[allIntPCchr2] = PCchr2LightIntensityAnalysis(1);
[allIntOBchr2] = OBchr2LightIntensityAnalysis(1);
% load the VS and SFC mean intensity values
load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\SFCsynchLevel');
load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\VSsynchLevel');
VSsynchLevel(1,:) = [1 0];
SFCsynchLevel(1,:) = [1 0];

lightInt = [0 10 25 30];%./30;%mW/mm^2
figure;
errorbar(lightInt,mean(allIntPCchr2),std(allIntPCchr2)./sqrt(size(allIntPCchr2,1)),'-om', 'linewidth',0.5, 'markerfacecolor', 'm')%PC
hold all;
errorbar(lightInt,mean(allIntOBchr2),std(allIntOBchr2)./sqrt(size(allIntOBchr2,1)),'-og', 'linewidth',0.5, 'markerfacecolor', 'g')%OB
errorbar(lightInt, SFCsynchLevel(:,1), SFCsynchLevel(:,2), '-or', 'linewidth',0.5, 'markerfacecolor', 'r')%
errorbar(lightInt, VSsynchLevel(:,1), VSsynchLevel(:,2), '-or', 'linewidth',0.5, 'markerfacecolor', 'r')%
plot([0 30], [1 1], 'k--')
ylabel('Normalized firing rate','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xlim([0 30])
hold off;
legend(['PC (N=' num2str(size(allIntPCchr2,1)) ')'], ['OB (N=' num2str(size(allIntOBchr2,1)) ')'])
ylim([0.4 2.2])

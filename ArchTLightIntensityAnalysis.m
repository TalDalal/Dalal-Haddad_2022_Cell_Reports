[allIntPCarch] = PCarchLightIntensityAnalysis;
[allIntOBarch] = OBarchLightIntensityAnalysis;

lightInt = [0 10 25 30];%mW/mm^2
figure;
errorbar(lightInt,mean(allIntPCarch),std(allIntPCarch)./sqrt(size(allIntPCarch,1)),'-og', 'linewidth',2)
hold on;
errorbar(lightInt,mean(allIntOBarch),std(allIntOBarch)./sqrt(size(allIntOBarch,1)),'-ok', 'linewidth',2)

plot([0 30], [1 1], 'k--')
ylabel('% firing rate change','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xlim([0 30])
hold off;
legend(['PC (N=' num2str(size(allIntPCarch,1)) ')'], ['OB (N=' num2str(size(allIntOBarch,1)) ')'])



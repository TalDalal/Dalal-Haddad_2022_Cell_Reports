%%%%% M/T suppressing intensity figure%%%%%%%
[allIntPCchr2] = MTarchPCLightIntensityAnalysis(1);
[allIntOBchr2] = MTarchOBLightIntensityAnalysis(1);


lightInt = [0 0.3 0.7 1].*18;% normalized intensity (1==~20mW/mm2)
figure;
errorbar(lightInt,mean(allIntPCchr2),std(allIntPCchr2)./sqrt(size(allIntPCchr2,1)),'-og', 'linewidth',2)
hold all;
errorbar(lightInt,mean(allIntOBchr2),std(allIntOBchr2)./sqrt(size(allIntOBchr2,1)),'-ok', 'linewidth',2)
plot([0 max(lightInt)+1], [1 1], 'k--')

ylabel('% firing rate change','fontSize',14)
xlabel('Light intensity (mW/mm^2)','fontSize',14)
set(gca,'fontSize',14)
box off;
set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
xlim(minmax(lightInt))
hold off;
legend(['PC (N=' num2str(size(allIntPCchr2,1)) ')'], ['OB (N=' num2str(size(allIntOBchr2,1)) ')'])

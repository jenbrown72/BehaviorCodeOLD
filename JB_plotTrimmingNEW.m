percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
positionGraph3 = [   680   110   756   868];
dataToAnalyse = {'performance';'dprime'};
sessionTypes = {'full pad';'None'};
thresholdLine = [0.5 0.5;1 1]
f=figure(22);clf
set(f,'Position',positionGraph3);

mrk = {'o','s','x'}';
for hh = 1:length(dataToAnalyse)
    
    dataAN = eval(dataToAnalyse{hh})';
    
    subplot(1,2,hh)
    h = plot(dataAN,'k');
    hold on
    set(h,{'marker'},mrk)
    ylabel(dataToAnalyse{hh});
    set(gca,'XTick',[1:length(sessionTypes)]);
    set(gca,'XTickLabel',sessionTypes)
    plot([0.5 2.5],thresholdLine(hh,:),'k--','LineWidth',2);
    
    if strcmp(dataToAnalyse{hh},'performance')
        ylimit = [0 1];
        set(gca,'ylim',ylimit)
    end
end

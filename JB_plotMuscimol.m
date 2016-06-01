percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
positionGraph3 = [   680   110   756   868];
dataToAnalyse = {'performance';'dprime'};
sessionTypes = {'control';'aCSF S1';'muscimol S1';'muscimol V1'};
thresholdLine = [0.5 0.5;1 1]
f=figure(22);clf
set(f,'Position',positionGraph3);

mrk = {'o','s','x'}';
for hh = 1:length(dataToAnalyse)
    
    subplot(1,3,hh)
    h = plot(eval(dataToAnalyse{hh}),'k');
    hold on
    set(h,{'marker'},mrk)
    ylabel(dataToAnalyse{hh});
    set(gca,'XTick',[1:length(sessions)]);
    set(gca,'XTickLabel',sessions)
    plot([min(xlim) max(xlim)],thresholdLine(hh,:),'k--','LineWidth',2);
    
    if strcmp(dataToAnalyse{hh},'performance')
        ylimit = [0 1];
        set(gca,'ylim',ylimit)
    end
end

markerLineWidth = 2;

for hh = 1:length(sessions)
    subplot(1,3,3)
    if strcmp(sessions{hh},'control')
        for kk = 1:size(FA,2)
       plot(FA(hh,kk),HIT(hh,kk),'k','marker',mrk{kk},'LineWidth',markerLineWidth)
        hold on
        end
    elseif strcmp(sessions{hh},'aCSF S1')
              for kk = 1:size(FA,2)
        plot(FA(hh,kk),HIT(hh,kk),'b','marker',mrk{kk},'LineWidth',markerLineWidth)
        hold on
              end
    elseif strcmp(sessions{hh},'muscimol S1')
              for kk = 1:size(FA,2)
        plot(FA(hh,kk),HIT(hh,kk),'r','marker',mrk{kk},'LineWidth',markerLineWidth)
        hold on
              end
    elseif strcmp(sessions{hh},'muscimol V1')
              for kk = 1:size(FA,2)
        plot(FA(hh,kk),HIT(hh,kk),'g','marker',mrk{kk},'LineWidth',markerLineWidth)
        hold on
              end
    end
    hline = refline(1,0);
    hline.Color = 'k';
    xlabel('false alarm rate')
    ylabel('hit rate')
end

figure(2)

for hh = 1:length(dataToAnalyse)
    
    data = eval(dataToAnalyse{hh});
    allDATA = nan(10,length(sessionTypes))
    
    for jj = 1:length(sessionTypes)
        idx = find(strcmp(sessions,sessionTypes{jj}))
        allDATA(1:(length(idx)*size(data,2)),jj) = reshape(data(idx,:),[(length(idx)*size(data,2)),1]);
        
    end
    
    average = nanmean(allDATA);
    std = nanstd(allDATA);
    sem = std./sqrt(sum(~isnan(allDATA),1));
    
    subplot(2,1,hh)
    bar(average,0.5);
    hold on
    errorbar(average,sem);
    set(gca,'XTickLabel',sessionTypes);
    ylabel(dataToAnalyse{hh})
end









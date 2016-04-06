function [analysedDATA] = JB_plotGroupStats(analysedDATA,plotON)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

positionGraph2 = [1321 548 587 448];
plotRows = 1;
plotCols = length(analysedDATA{1}.parameter) ;
currPlot = 1;
Labels = {'Full', 'row/column','double','single','None'};

if (plotON==1)
    ffff=figure;clf
    set(ffff,'Position',positionGraph2);
else
    figure('Visible','off');clf;
end

for i = 1:length(analysedDATA{1}.parameter)
    k=1;
    clear p h temp x sigPairsPval sigPairs;
    for kk = 1:size((analysedDATA{1}.parameter{i}.tempDATA),1)
    
    [p(kk),h(kk)] = ranksum(analysedDATA{1}.parameter{i}.tempDATA(kk,:),analysedDATA{2}.parameter{i}.tempDATA(kk,:))
    
       if h(kk)==1
                sigPairs{k} = [kk,kk];
                sigPairsPval(k) = [p(kk)];
                k = k+1;
       end

    end

    subplot(plotRows,plotCols,currPlot);

    for j = 1:length(analysedDATA)

        toPlot(:,j) = (mean(analysedDATA{j}.parameter{i}.tempDATA,2));
         toPlotsem(:,j) = (std(analysedDATA{j}.parameter{i}.tempDATA,0,2));
    end
    
    h = bar(toPlot);
    colormap(gray)
    set(h, 'BarWidth',1);
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'XTickLabel',Labels);
    set(get(gca,'YLabel'),'String',(analysedDATA{1}.parameter{i}.name))
    numbers = size(toPlot,2);
    numgroups = size(toPlot,1);
    groupWidth = min(0.8,numbers/(numbers+1.5));
    hold on
    for i=1:numbers
        x(i,:) = (1:numgroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numbers);
        errorbar(x(i,:),toPlot(:,i),toPlotsem(:,i),'k','Linestyle','none');
    end
    
    temp = sigPairs;
    
    for gg = 1:length(sigPairs)
        
      temp{gg} = [x(1,sigPairs{gg}(1)) x(2,sigPairs{gg}(1))]
    end
    
      sigstar(temp,sigPairsPval,0,0.3);
    
    currPlot = currPlot+1;
    
end




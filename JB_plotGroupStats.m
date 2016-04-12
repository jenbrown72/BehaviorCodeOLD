function [] = JB_plotGroupStats(analysedDATA,normdata,diffdata,plotON)
%UNTITLED6 Summary of this function goes here
%   analysedDATA matrix is generated from : [analysedDATA] = JB_plotGroupAverages(AllDATA,plotON)
%   norm = 1; normalise data norm = 0; raw data
%   stats on the difference
%   plotON =1; plot, plotON = 0, no plot

positionGraph2 = [1095         548         813         448];
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

for i = 1:length(analysedDATA{1}.parameter);

    tempTitle = ['AverageData', ' ', 'Norm', num2str(normdata),' ', 'Diff', num2str(diffdata)];
    set(gcf,'name',tempTitle,'numbertitle','off')
    
    for ff=1:length(analysedDATA)
        data{ff} = analysedDATA{ff}.parameter{i}.tempDATA;
    end
    
    %normalised data
    
    if normdata==1;
        for ff = 1:size((data),2)
            data{ff} = bsxfun(@rdivide,data{ff}, max(data{ff}));
        end
    end

    if (diffdata==1)
          for ff = 1:size((data),2)
                for ffg = 1:size((data{ff}),2)
              for fg = 1:size((data{ff}),1)-1
                
                      
            tempdata{ff}(fg,ffg)           = diff([data{ff}(1,ffg) data{ff}(fg+1,ffg)]);
          %  data{ff}(fg,ffg) = diff([data{ff}(1,ffg) data{ff}(fg+1,ffg)]);
                  end
              end
          end
          
        Labels = {'Full:R/C', 'Full:double','Full:single','Full:None'};
    end
    data = tempdata;
    
    k=1;
    clear p h temp x sigPairsPval sigPairs;
    p = nan(1,size(data{1},1));
    h = nan(1,size(data{1},1));
    sigPairsPval = nan(1,size(data{1},1));
    sigPairs = cell(1,size(data{1},1));
    %for each experimental condition (e.g. Full, 4 whiskers etc, find the
    %ranksum to see if the performance/d prime is different between the two groups
    
    for kk = 1:size(data{1},1)
        [p(kk),h(kk)] = ranksum(data{1}(kk,:),data{2}(kk,:));
        if h(kk)==1;
            sigPairs{k} = [kk,kk];
            sigPairsPval(k) = p(kk);
            k = k+1;
        end
    end
    
    sigPairs = sigPairs(~cellfun('isempty',sigPairs));
    sigPairsPval(isnan(sigPairsPval(1,:))) =[];
    
    subplot(plotRows,plotCols,currPlot);
    averageData = nan(size(data{1},1),length(analysedDATA));
    semData = nan(size(data{1},1),length(analysedDATA));
    for j = 1:length(analysedDATA);
        averageData(:,j) = (mean(data{j},2));
        semData(:,j) = (std(data{j},0,2));
    end
    
    h = bar(averageData);
    colormap(gray)
    set(h, 'BarWidth',1);
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'XTickLabel',Labels);
    set(get(gca,'YLabel'),'String',(analysedDATA{1}.parameter{i}.name))
    numbers = size(averageData,2);
    numgroups = size(averageData,1);
    groupWidth = min(0.8,numbers/(numbers+1.5));
    hold on
    x = nan(numbers,numgroups);
    for ii=1:numbers
        x(ii,:) = (1:numgroups) - groupWidth/2 + (2*ii-1) * groupWidth / (2*numbers);
        errorbar(x(ii,:),averageData(:,ii),semData(:,ii),'k','Linestyle','none');
    end
    
    temp = sigPairs;
    for gg = 1:length(sigPairs)
        temp{gg} = [x(1,sigPairs{gg}(1)) x(2,sigPairs{gg}(1))];
    end
    sigstar(temp,sigPairsPval,0,0.2);
    currPlot = currPlot+1;
    
end




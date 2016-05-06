function [] = JB_plotGroupLearningT(AllDATA,plotON,neg)
%UNTITLED5 Summary of this function goes here
%   AllDATA matrix is generated from : [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)
%   norm = 1; normalise data norm = 0; raw data
%   exclude = 1 - exclude mice depending on threshold set below
% trimSession=1 (first timetrimmed), =2 (secondtimeTrimmed)
%   plotON =1; plot, plotON = 0, no plot

%load('AllDATA.mat');

positionGraph1 = [14    42   353   954];
positionGraph2 = [380    42   353   954];
dprimeThreshold = [1 1];
percentCorrectChance = [0.5 0.5];
markerSizePlots=6;

if (plotON==1)
    f1 = figure(1);clf
    set(f1,'Position',positionGraph1);
    set(f1,'name','learningPerformance','numbertitle','off');
    
    f2 = figure(2);clf
    set(f2,'Position',positionGraph2);
    set(f2,'name','learningDprime','numbertitle','off');
else
    figure('Visible','off');clf;
end

thresholds = [0.5 0.5;1 1];
dataToAnalyse = {'sessionperformance';'dprime'};

noSubPlots = length(AllDATA);

for gg = 1:length(dataToAnalyse)
    for kk = 1:length(AllDATA)
        basicPropertiesToPlot = AllDATA{1,kk};
        numPoints = 1:1:length(basicPropertiesToPlot);
        
        for j = 1:length(basicPropertiesToPlot);
            if strcmp('S1auto', basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 1;
            elseif strcmp('S1',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 2;
            elseif strcmp('S2',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 3;
            elseif strcmp('S6',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 4;
            elseif strcmp('S8',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 5;
            elseif strcmp('S10',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 6;
            elseif strcmp('S12',basicPropertiesToPlot{j,1}.sessionType)
                colorCode(j,1) = 7;
            end
        end
        
        for j = 1:length(numPoints);
            negReinforcer(j,1) = basicPropertiesToPlot{j,1}.negReinforcer;
            waterSchedule(j,1) = basicPropertiesToPlot{j,1}.waterSchedule;
            optogenetics(j,1) = basicPropertiesToPlot{j,1}.optogenetics;
        end
        
        ColorCodeColors = [0.9;0.8;0.6;0.4;0.2;0];
        
        for jj = 1:length(basicPropertiesToPlot);
            dataTemp(jj,:)  = [basicPropertiesToPlot{jj,1}.dprime colorCode(jj,1) negReinforcer(jj,1)];
        end
        
        sessionNoStop = (find((dataTemp(:,2)==5),1))+1;
        
        %plot performance
        diffData = diff(dataTemp(:,2));
        figure(gg);
        subplot(noSubPlots,1,kk)
        set(gca,'xtick',[])
        for jj = 1:sessionNoStop-1;
            if diffData(jj)==1;
                if ((negReinforcer(jj)==1) && (neg==1));
                    plot(numPoints(jj),basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}),'o','MarkerEdgeColor',[negReinforcer(jj),0,0],'MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
                    hold on
                else
                    plot(numPoints(jj),basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}),'o','MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
                    hold on
                end
            else
                if ((negReinforcer(jj)==1) && (neg==1));
                    plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}) basicPropertiesToPlot{jj+1,1}.(dataToAnalyse{gg})],'o-','MarkerEdgeColor',[negReinforcer(jj),0,0],'MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
                    hold on
                else
                    plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}) basicPropertiesToPlot{jj+1,1}.(dataToAnalyse{gg})],'o-','MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
                    hold on
                end
            end
        end
        plot([0 24],thresholds(gg,:),'k--')
        
        if strcmp((dataToAnalyse{gg}),'sessionperformance')
            ylim([0 1])
        elseif strcmp((dataToAnalyse{gg}),'dprime')
            ylim([-1 5])
        end
        xlim([0 24])
        xlimit = xlim;
        text(xlimit(1),max(ylim)+0.1,basicPropertiesToPlot{jj,1}.mouseID)
        if (kk==noSubPlots)
            xlabel('Session Number');
            ylabel((dataToAnalyse{gg}));
        else
            set(gca,'XTickLabel',[])
        end
    end
end
end


function [sessionProperties,sessionsToCriterion] = JB_plotGroupLearning(AllDATA,plotON,neg)
%UNTITLED5 Summary of this function goes here
%   AllDATA matrix is generated from : [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)
%   neg = 1; include marker for when negative reinforment started
%   plotON =1; plot, plotON = 0, no plot

positionGraph1 = [14    42   353   954];
positionGraph2 = [380    42   353   954];
positionGraph3 = [735   741   879   255];
positionGraph4 = [1349          44         557         952];
markerSizePlots=6;
percentCorrectChance = [0.5 0.5];
tally2=1;


if (plotON==1)
    f1 = figure(1);clf
    set(f1,'Position',positionGraph1);
    set(f1,'name','learningPerformance','numbertitle','off');
    
    f2 = figure(2);clf
    set(f2,'Position',positionGraph2);
    set(f2,'name','learningDprime','numbertitle','off');
    
    f3 = figure(3);clf
    set(f3,'Position',positionGraph3);
    set(f3,'name','sessionProperties','numbertitle','off');
    
    f4 = figure(4);clf
    set(f4,'Position',positionGraph4);
    set(f4,'name','sessionProperties','numbertitle','off');
else
    figure('Visible','off');clf;
end

thresholds = [0.5 0.5;1 1];
dataToAnalyse = {'sessionperformance';'dprime'};

noSubPlots = length(AllDATA);

%%Plot Learning

for kk = 1:length(AllDATA)
    clear sessionCode;
    basicPropertiesToPlot = AllDATA{1,kk};
    numPoints = 1:1:length(basicPropertiesToPlot);
    
    for j = 1:length(basicPropertiesToPlot);
        if strcmp('S1auto', basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 1;
        elseif strcmp('S1',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 2;
        elseif strcmp('S2',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 3;
        elseif strcmp('S6',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 4;
        elseif strcmp('S8',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 5;
        elseif strcmp('S10',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 6;
        elseif strcmp('S12',basicPropertiesToPlot{j,1}.sessionType)
            sessionCode(j,1) = 7;
        end
        
        blockPerformance = basicPropertiesToPlot{j,1}.temptrialTypes;
        
        blockPerformance(blockPerformance==4)=1; %CR
        blockPerformance(blockPerformance==1)=1; %Hit
        blockPerformance(blockPerformance==2)=0; %FA
        blockPerformance(blockPerformance==3)=0; %Miss
        
        tally = 0;
        basicPropertiesToPlot{j,1}.rolloingperformance = [];
        for kj = 1:length(blockPerformance)
            basicPropertiesToPlot{j,1}.rolloingperformance(kj,1) = (blockPerformance(kj)+tally)/kj;
            tally = tally+blockPerformance(kj);
        end
        
    end
    
    for j = 1:length(numPoints);
        negReinforcer(j,1) = basicPropertiesToPlot{j,1}.negReinforcer;
        waterSchedule(j,1) = basicPropertiesToPlot{j,1}.waterSchedule;
        optogenetics(j,1) = basicPropertiesToPlot{j,1}.optogenetics;
    end
    
    ColorCodeColors = [1;0.9;0.8;0.6;0.4;0.2;0];
    sessionNoStop = (find((sessionCode(:,1)==5),1))+1;
    
    %plot each parameter
    for gg = 1:length(dataToAnalyse)
        diffData = diff(sessionCode(:,1));
        figure(gg);
        subplot(noSubPlots,1,kk)
        set(gca,'xtick',[])
        for jj = 1:sessionNoStop-1;
            if diffData(jj)==1;
                if ((negReinforcer(jj)==1) && (neg==1));
                    plot(numPoints(jj),basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}),'o','MarkerEdgeColor',[negReinforcer(jj),0,0],'MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'Color', [0.4,ColorCodeColors(sessionCode(jj)),0.8]);
                    hold on
                else
                    plot(numPoints(jj),basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}),'o','MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'Color', [0.4,ColorCodeColors(sessionCode(jj)),0.8]);
                    hold on
                end
            else
                if ((negReinforcer(jj)==1) && (neg==1));
                    plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}) basicPropertiesToPlot{jj+1,1}.(dataToAnalyse{gg})],'o-','MarkerEdgeColor',[negReinforcer(jj),0,0],'MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'Color', [0.4,ColorCodeColors(sessionCode(jj)),0.8]);
                    hold on
                else
                    plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.(dataToAnalyse{gg}) basicPropertiesToPlot{jj+1,1}.(dataToAnalyse{gg})],'o-','MarkerSize', markerSizePlots, 'linewidth',1, 'MarkerFaceColor',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'Color', [0.4,ColorCodeColors(sessionCode(jj)),0.8]);
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
    
    spaceBetweenSessions = 40;
    figure(f4);
    subplot(noSubPlots,1,kk)
    
    cummX = 1;
    for jj = 1:sessionNoStop;
        %            for jj = 1:length(basicPropertiesToPlot)
        if ~isempty(basicPropertiesToPlot{jj,1}.rolloingperformance) && (length(basicPropertiesToPlot{jj,1}.rolloingperformance)>100);
            movingAvg=movingmean( basicPropertiesToPlot{jj,1}.rolloingperformance,61);
            plot(cummX:cummX+length(movingAvg)-1,movingAvg,'-', 'Color',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'LineWidth',4)
            hold on
            xlimit = xlim;
            cummX = cummX + length(movingAvg) + spaceBetweenSessions;
            text(xlimit(1),max(ylim)+0.1,basicPropertiesToPlot{jj,1}.mouseID)
        end
    end
    ylim([0 1])
    xlim([0 5000])
    plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',1);
    xlimit = xlim;
    if (kk==noSubPlots)
        ylabel('Moving fraction Correct');
        xlabel('Trial No.');
    else
        set(gca,'XTickLabel',[])
    end
    %     ylabel('Moving fraction Correct')
    %     xlabel('Trial No.')
    
    sessionToinclude = (find(sessionCode(:,1)==5));
    sessionsToCriterion(kk,1) = sessionToinclude(1);
    
    for j = 1:length(sessionToinclude);
        sessionProperties(tally2,:) = [basicPropertiesToPlot{sessionToinclude(j),1}.noTrials basicPropertiesToPlot{sessionToinclude(j),1}.sessionDuration];
        tally2 = tally2+1;
    end
    
end

figure(3)
subplot(1,3,1)
[N,X] = hist(sessionProperties(:,1),15);
bar(X,N,1,'FaceColor','w','EdgeColor','k')
xlabel('number of Trials/session');
ylabel('session count');

subplot(1,3,2)
[N,X] = hist(sessionProperties(:,2),15);
bar(X,N,1,'FaceColor','w','EdgeColor','k')
xlabel('duration of session (mins)');
ylabel('session count');
xlimit = xlim;
strText = ['No. session',' ' , num2str(tally2)];
text(((xlimit(2)/5)*2.5),((max(ylim)/10)*9),strText)
strText = ['No. mice',' ' , num2str(kk)];
text(((xlimit(2)/5)*2.5),((max(ylim)/10)*8),strText)

subplot(1,3,3)
sortedDATA = sort(sessionsToCriterion);
fractionMice = cumsum(sort(sessionsToCriterion))/max(cumsum(sort(sessionsToCriterion)));
[~,~,H] = unique(sortedDATA);
idx = find(diff(H));
idx(end+1) = length(H);
plot(sortedDATA(idx),fractionMice(idx),'o-k','MarkerFaceColor','w');
ylim([0 1])
xlabel('sessions to criterion');
ylabel('cum. frac. of mice');

end


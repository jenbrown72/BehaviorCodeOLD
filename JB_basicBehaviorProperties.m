function [] = JB_basicBehaviorProperties(plotfig)

dprimeChance = [1 1];
dprimeThreshold = [1 1];
percentCorrectChance = [0.5 0.5];
percentCorrectThreshold = [0.7 0.7];
positionGraph1 = [14   83   503   869];
positionGraph2 = [602   83   644   869];
positionGraph3 = [1268   83   644   869];
positionGraph4 = [32   515   354   438];
positionGraph5 = [1214  230  580 736];
positionGraph6 = [602    83   644   869];
positionGraph7 = [754    78   327   869];
plotON=1;

if nargin==1;
    plotON = plotfig;
end
%%
load('DATA.mat')
sessionType = 1;
basisPropertiesID{sessionType,1} = 'sessionType';
date = 2;
basisPropertiesID{date,1} = 'date';
datenum = 3;
basisPropertiesID{datenum,1} = 'datenum';
sessionDuration = 4;
basisPropertiesID{sessionDuration,1} = 'sessionDuration';
totalStepsPerMin = 5;
basisPropertiesID{totalStepsPerMin,1} = 'totalStepsPerMin';
performance = 6;
basisPropertiesID{performance,1} = 'performance';
HIT = 7;
basisPropertiesID{HIT,1} = 'performanceTypeHIT';
FA = 8;
basisPropertiesID{FA,1} = 'performanceTypeFA';
MISS = 9;
basisPropertiesID{MISS,1} = 'performanceTypeMISS';
CR = 10;
basisPropertiesID{CR,1} = 'performanceTypeCR';
possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];

startAngleID = 10;
basisPropertiesID{startAngleID,1} = 'anglePerformanceHIT';
% FAangle = 12;
% basisPropertiesID{FAangle,1} = 'anglePerformanceFA';
% MISSangle = 13;
% basisPropertiesID{MISSangle,1} = 'anglePerformanceMISS';
% CRangle = 14;
% basisPropertiesID{CRangle,1} = 'anglePerformanceCR';

basicProperties = {};
encoder0Pos = 1;
LEDstate = 5;
count = 6;
rawSessionTime = 7;
trialType = 13;
angle = 11;
rewardDuration = 9;

for i=1:length(DATA.allFiles);
    
    metaData = DATA.allFiles{i}.metaData;
    
    %find out if this session included a negative stim
    
    puffIdx = find(strcmp('Negative Reinforcer = ', metaData));
    puffUsed = metaData{puffIdx,2};
    if ischar(puffUsed)==1
        basicProperties{i,1}.negReinforcer = str2num(puffUsed);
    else
        basicProperties{i,1}.negReinforcer = puffUsed;
    end
    
    
    %find out if this session included a water schedule
    waterSchIdx = find(strcmp('Water Schedule = ', metaData));
    waterSchUsed = metaData{waterSchIdx,2};
    if ischar(waterSchUsed)==1
        basicProperties{i,1}.waterSchedule = str2num(waterSchUsed);
        
    else
        
        basicProperties{i,1}.waterSchedule = waterSchUsed;
        
    end
    
    %find out if this session used optogenetics
    optoIdx = find(strcmp('Optogenetics = ', metaData));
    tempOpto = metaData{optoIdx,2};
    if ischar(tempOpto)==1
        basicProperties{i,1}.optogenetics = str2num(tempOpto);
    else
        basicProperties{i,1}.optogenetics = tempOpto;
    end
    
    
    %find out if this session was a trimming session
    
    trimIdx = find(strcmp('Trimming = ', metaData));
    if ~isempty(trimIdx)
        
        tempTrim = metaData{trimIdx,2};
        basicProperties{i,1}.trimming = tempTrim;
        
    else
        basicProperties{i,1}.trimming = [];
    end
    
    %find out if this session had segments
    segIdx = find(strcmp('Segments = ', metaData));
    if ~isempty(segIdx)
        
        tempSeg = metaData{segIdx,2};
        tempSegNum = str2num(tempSeg);
        
        if (tempSegNum(1)==1) && (tempSegNum(2)==1)
            
            tempData = DATA.allFiles{i}.rawData;
            index = find(all(tempData==1,2));
            index = [index; length(tempData)]; %add a one to find the start of the first segment
            
            segments(1,1) = 1; %add first trial
            for kk=1:length(index)-1;
                
                countCummulative = tempData(index(kk)-1,6);
                tempData(index(kk)+1:index(kk+1),6) = tempData(index(kk)+1:index(kk+1),6)+countCummulative;
                basicProperties{i,1}.segments(kk,2) = countCummulative;
                basicProperties{i,1}.segments(kk+1,1) = countCummulative+1;
                
            end
            basicProperties{i,1}.segments(length(basicProperties{i,1}.segments),2) = tempData(length(tempData),6); %add last trial
            
        else
            basicProperties{i,1}.segments = str2num(tempSeg);
        end
    else
        basicProperties{i,1}.segments = [];
    end
    
    basicProperties{i,1}.date = DATA.allFiles{i}.date;
    basicProperties{i,1}.datenum = DATA.allFiles{i}.dateFromFile;
    
    tempName = DATA.allFiles{i}.name;
    startIdx = findstr('_201',tempName);
    endIdx = findstr('_Box',tempName);
    tempName = tempName(startIdx:endIdx);
    tempLocation = findstr(tempName,'_');
    for hh = 1:length(tempLocation)-1;
        if(tempLocation(hh+1)-tempLocation(hh))==2;
            startName = tempName(1:tempLocation(hh));
            endName = tempName(tempLocation(hh)+1:end);
            tempName = strcat(startName,'0',endName);
            tempLocation = tempLocation+1;
        end
    end
    
    toDelete = strfind(tempName,'_');
    tempName(toDelete)='.';
    
    basicProperties{i,1}.namedata = tempName;
    
    tempDATA = DATA.allFiles{i}.rawData;
        basicProperties{i,1}.sessionDuration = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
   basicProperties{i,1}.sessionDuration = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
   basicProperties{i,1}.fileLength = length(tempDATA);
   threshold = 50000;
    aboveThreshold = find(tempDATA(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        tempDATA(aboveThreshold(j),:)=nan;
    end
    
    totalSteps = nansum(diff(tempDATA(:,1))>0);
    basicProperties{i,1}.totalStepsPerMin = totalSteps/basicProperties{i,1}.sessionDuration; % Calculate running velocity
    
    
    %Look at raw data to get stats on performance
    findtrialType = diff(DATA.allFiles{i}.rawData(:,trialType));
    
    [idx, ~] = find(findtrialType>0);
    temptrialTypes = DATA.allFiles{i}.rawData(idx+2,trialType);
    tempCount = DATA.allFiles{i}.rawData(idx,count);
    tempAngleTypes = DATA.allFiles{i}.rawData(idx+2,angle);
    tempRewardDuration = DATA.allFiles{i}.rawData(idx+2,rewardDuration);
    tempOptogenetics = DATA.allFiles{i}.rawData(idx-10,LEDstate);
    % angleUsed = unique(tempAngleTypes);
    % [idx] = find(angleUsed>0); %no 0
    
    %if this was a segmented session, calculate performance in each session
    
    if ~isempty( basicProperties{i,1}.segments)
        
        sessionDivisions = basicProperties{i,1}.segments;
        segmentPerformance = nan(length(sessionDivisions),2);
        
        for j = 1:length(sessionDivisions)
            
            if sessionDivisions(j,1)==0;
                idxStart = 1;
            else
                idxStart = sessionDivisions(j,1);
            end
            idxEnd = sessionDivisions(j,2);
            
            %     idxEnd = find(tempCount-sessionDivisions(j,2)< 0, 1,'last');
            
            tempSegmentData = temptrialTypes(idxStart:idxEnd,:);
            HITangle = sum(tempSegmentData(:,1)==1)+1;
            FAangle = sum(tempSegmentData(:,1)==2)+1;
            MISSangle = sum(tempSegmentData(:,1)==3)+1;
            CRangle = sum(tempSegmentData(:,1)==4)+1;
            basicProperties{i,1}.segmentPerformance(j,1) =  (HITangle+CRangle)/(length(tempSegmentData)+4);
            basicProperties{i,1}.segmentdPrime(j,1)=norminv(( HITangle/( HITangle+MISSangle)),0,1)-norminv((FAangle/(FAangle+CRangle)),0,1);
            
        end
        
    else
        basicProperties{i,1}.segmentPerformance = [];
        basicProperties{i,1}.segmentdPrime = [];
    end
    
    for kk=1:length(possibleAngles)
        %   clear performanceAngleTempSTIM performanceAngleTempNoSTIM;
        
        stimCount = 1;
        nostimCount = 1;
        
        idxAngles = find(tempAngleTypes==possibleAngles(kk));
        performanceAngleTemp = temptrialTypes(idxAngles);
        
        if ~isempty(idxAngles)
            if basicProperties{i,1}.optogenetics==1
                for ff = 1:length(idxAngles)
                    if (tempOptogenetics(idxAngles(ff))==1)
                        performanceAngleTempSTIM(stimCount,1) = temptrialTypes(idxAngles(ff));
                        stimCount = stimCount+1;
                    else
                        performanceAngleTempNoSTIM(nostimCount,1) = temptrialTypes(idxAngles(ff));
                        nostimCount = nostimCount+1;
                    end
                end
                
            else
                performanceAngleTempSTIM = nan;
                performanceAngleTempNoSTIM = nan;
            end
        end
        
        HITangle = sum(performanceAngleTemp==1);
        FAangle = sum(performanceAngleTemp==2);
        MISSangle = sum(performanceAngleTemp==3);
        CRangle = sum(performanceAngleTemp==4);
        %Performance at each angle
        basicProperties{i,1}.performance{kk} =  (HITangle+CRangle)/(length(performanceAngleTemp));
        %Probability of licking
        basicProperties{i,1}.probLicking{kk} =  (HITangle+FAangle)/(MISSangle+CRangle+HITangle+FAangle);
        
        if ~isempty(idxAngles)
            if (basicProperties{i,1}.optogenetics==1)
                HITangleSTIM = sum(performanceAngleTempSTIM==1);
                FAangleSTIM = sum(performanceAngleTempSTIM==2);
                MISSangleSTIM = sum(performanceAngleTempSTIM==3);
                CRangleSTIM = sum(performanceAngleTempSTIM==4);
                basicProperties{i,1}.performanceSTIM{kk} =  (HITangleSTIM+CRangleSTIM)/(length(performanceAngleTempSTIM));
                basicProperties{i,1}.probLickingSTIM{kk} =  (HITangleSTIM+FAangleSTIM)/(MISSangleSTIM+CRangleSTIM+HITangleSTIM+FAangleSTIM);
                
                HITanglenoSTIM = sum(performanceAngleTempNoSTIM==1);
                FAanglenoSTIM = sum(performanceAngleTempNoSTIM==2);
                MISSanglenoSTIM = sum(performanceAngleTempNoSTIM==3);
                CRanglenoSTIM = sum(performanceAngleTempNoSTIM==4);
                basicProperties{i,1}.performanceNoSTIM{kk} =  (HITanglenoSTIM+CRanglenoSTIM)/(length(performanceAngleTempNoSTIM));
                basicProperties{i,1}.probLickingNoSTIM{kk} =  (HITanglenoSTIM+FAanglenoSTIM)/(MISSanglenoSTIM+CRanglenoSTIM+HITanglenoSTIM+FAanglenoSTIM);
                
                
            else
                basicProperties{i,1}.performanceSTIM{kk} = nan;
                basicProperties{i,1}.probLickingSTIM{kk} = nan;
                basicProperties{i,1}.performanceNoSTIM{kk} = nan;
                basicProperties{i,1}.probLickingNoSTIM{kk} = nan;
            end
            
            
        else
            
            basicProperties{i,1}.performanceSTIM{kk} = nan;
            basicProperties{i,1}.probLickingSTIM{kk} = nan;
            basicProperties{i,1}.performanceNoSTIM{kk} = nan;
            basicProperties{i,1}.probLickingNoSTIM{kk} = nan;
        end
        
    end
    
    basisPropertiesID{startAngleID,1};
    
    %Calculate how many of each trial type there were
    basicProperties{i,1}.HIT = sum(temptrialTypes==1)+1; %HitTrialCount
    basicProperties{i,1}.FA = sum(temptrialTypes==2)+1; % FATrialCount
    basicProperties{i,1}.MISS = sum(temptrialTypes==3)+1; %MissTrialCount
    basicProperties{i,1}.CR = sum(temptrialTypes==4)+1; %CRTrialCount
    
    basicProperties{i,1}.sessionperformance = (basicProperties{i,1}.HIT+basicProperties{i,1}.CR)/(length(temptrialTypes)+4);
    
    %        %D' - discrimitability index
    %when H=F, then d' =0
    %Highest Possible d' (greatest sensitivity) is 6.93, the effective limit
    %(using .99 and .01) 4.65, typical values are up to 2 and 69% correct for
    %both different and same trials corresponds to a d' of 1.0.
    
    
    basicProperties{i,1}.dprime=norminv(( basicProperties{i,1}.HIT/( basicProperties{i,1}.HIT+basicProperties{i,1}.MISS)),0,1)-norminv((basicProperties{i,1}.FA/(basicProperties{i,1}.FA+basicProperties{i,1}.CR)),0,1);
    
    %now look at optogenetic sessions
    
    if (basicProperties{i,1}.optogenetics==1)
        
        basicProperties{i,1}.HITStim = sum((temptrialTypes==1)&(tempOptogenetics==1))+1; %HitTrialCount
        basicProperties{i,1}.FAStim = sum((temptrialTypes==2)&(tempOptogenetics==1))+1; % FATrialCount
        basicProperties{i,1}.MISSStim = sum((temptrialTypes==3)&(tempOptogenetics==1))+1; %MissTrialCount
        basicProperties{i,1}.CRStim = sum((temptrialTypes==4)&(tempOptogenetics==1))+1; %CRTrialCount
        
        basicProperties{i,1}.HITNoStim = sum((temptrialTypes==1)&(tempOptogenetics==0))+1; %HitTrialCount
        basicProperties{i,1}.FANoStim = sum((temptrialTypes==2)&(tempOptogenetics==0))+1; % FATrialCount
        basicProperties{i,1}.MISSNoStim = sum((temptrialTypes==3)&(tempOptogenetics==0))+1; %MissTrialCount
        basicProperties{i,1}.CRNoStim = sum((temptrialTypes==4)&(tempOptogenetics==0))+1; %CRTrialCount
        
        basicProperties{i,1}.dprimeSTIM=norminv(( basicProperties{i,1}.HITStim/( basicProperties{i,1}.HITStim+basicProperties{i,1}.MISSStim)),0,1)-norminv((basicProperties{i,1}.FAStim/(basicProperties{i,1}.FAStim+basicProperties{i,1}.CRStim)),0,1);
        basicProperties{i,1}.dprimeNoSTIM=norminv(( basicProperties{i,1}.HITNoStim/( basicProperties{i,1}.HITNoStim+basicProperties{i,1}.MISSNoStim)),0,1)-norminv((basicProperties{i,1}.FANoStim/(basicProperties{i,1}.FANoStim+basicProperties{i,1}.CRNoStim)),0,1);
        basicProperties{i,1}.sessionperformanceSTIM = (basicProperties{i,1}.HITStim+basicProperties{i,1}.CRStim)/(sum(tempOptogenetics==1)+4);
        basicProperties{i,1}.sessionperformanceNoSTIM = (basicProperties{i,1}.HITNoStim+basicProperties{i,1}.CRNoStim)/(sum(tempOptogenetics==0)+4);
    end
    
    %import as metaDATA cell array
    n=1;
    
    
    
    for j=1:length(metaData)
        
        %         match(j,1) = strcmp('Orientation Selected = ', metaData{j,1});
        if(strcmp('Auto reward = ', metaData{j,1}));
            autoSession = metaData{j,2};
        end
        
    end
    
    match = strcmp('Orientation Selected = ', metaData(:,1));
    
    %     angleUsed = unique(DATA.allFiles{1,i}.rawData(:,11));
    %       tempAngles = angleUsed(angleUsed>200);
    clear ind
    %  if iscell(metaData(match,2))==1
    
    for kk=1:length(match)
        if match(kk,1)==1;
            
            if isstr(metaData{kk,2})
                ind(kk,1) = str2num(metaData{kk,2});
                
            else
                
                ind(kk,1) =(metaData{kk,2});
            end
            
            
            
        end
    end
    %  else
    
    %[ind] = str2num(char(metaData(match,2)));
    
    %   end
    tempAngles = ind(ind>0);
    %         if (match(j,1)==1) && (metaData{j,2}>0)
    %             tempAngles(n,1) = metaData{j,2};
    %             n=n+1;
    %         end
    %
    %         %AutoReward
    %         if(strcmp('Auto reward = ', metaData{j,1}));
    %             autoSession = metaData{j,2};
    %         end
    %
    %     end
    
    %Define session type
    tempsessionType = strcat('S',num2str(length(tempAngles)));
    basicProperties{i,1}.stimNumber = length(tempAngles);
    if autoSession==1
        tempsessionType = strcat(tempsessionType,'auto');
        basicProperties{i,1}.stimNumber = 0;
    end
    basicProperties{i,1}.sessionType = tempsessionType;
    
end

for k = 1:length(basicProperties),
    
    days{k,1} = basicProperties{k,1}.datenum(1:8);
    
end

plotNumber = 1;

[un idx_last idx] = unique(days(:,1));
idx = sort(idx);
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});

for k = 1:length(unique_idx)
    tempRunVal=[];
    if length(unique_idx{k})==1
        basicPropertiesToPlot{plotNumber,1} = basicProperties{unique_idx{k}};
        plotNumber = plotNumber+1;
    elseif length(unique_idx{k})>1
        t = unique_idx{k};
        for kk = 1:length(t)
            tempRunVal(kk,1) = basicProperties{t(kk)}.fileLength;
        end
        
        [id di] = max(tempRunVal); %find the data set where the mouse ran the most
        basicPropertiesToPlot{plotNumber,1} = basicProperties{t(di)};
        plotNumber = plotNumber+1;
    end
end
%%

if plotON==1
    
    if (plotON==1)
        f = figure;clf
        set(f,'Position',positionGraph1);
    else
        figure('Visible','off');clf;
    end
    
    tempTitle = DATA.mouseID;
    tempTitle(findstr(tempTitle,'_'))=[];
    set(gcf,'name',tempTitle,'numbertitle','off');
    
    noSubPlots = 4;
    subplot(noSubPlots,1,1)
    numPoints = 1:1:length(basicPropertiesToPlot);
    for j = 1:length(numPoints);
        subplot(noSubPlots,1,1)
        plot(numPoints(j),basicPropertiesToPlot{j}.totalStepsPerMin,'or','MarkerSize', 10,'MarkerFaceColor','r')
        hold on
    end
    
    ylabel('totalStepsPerMin');
    xlabel('Session Number');
    %
    % subplot(noSubPlots,1,2)
    % for j = 1:length(numPoints);
    %     plot(numPoints(j),basicPropertiesToPlot{j}.sessionDuration,'or','MarkerSize', 10,'MarkerFaceColor','r')
    %     hold on
    % end
    %
    % ylabel('sessionDuration');
    % xlabel('Session Number');
    
    %plot primary session type per day
    subplot(noSubPlots,1,2)
    SessionTypes = {'S1auto' ; 'S1'; 'S2'; 'S6'; 'S8'; 'S10'; 'S12'};
    
    for j = 1:length(numPoints);
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
    
    ColorCodeColors = [1.0;0.9;0.8;0.6;0.4;0.2;0];
    
    for j = 1:length(numPoints);
        if waterSchedule(j)==1;
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(j)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        elseif negReinforcer(j)==1;
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(j),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        else
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(j)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        end
        hold on
        if optogenetics(j)==1;
            plot(numPoints(j),colorCode(j,1),'*y','MarkerSize', 5)
        end
    end
    set(gca,'YTickLabel',SessionTypes)
    ylabel('Session Type');
    xlabel('Session Number');
    
    %plot performance
    subplot(noSubPlots,1,3)
    addtoLegend = 1;
    addedWater = 0;
    addedNeg = 0;
    addedOpto = 0;
    
    legendAdd = [];
    legendTab = [];
    for jj = 1:length(numPoints);
        
        if waterSchedule(jj)==1;
            circle1 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(jj)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedWater==0)
                legendAdd(addtoLegend,:) = circle1;
                legendTab{addtoLegend} = 'water schedule';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        elseif negReinforcer(jj)==1;
            circle2 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(jj),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = circle2;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        else
            plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        end
        hold on
        
        if optogenetics(jj)==1;
            circle3 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'*y','MarkerSize', 5);
            hold on
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = circle3;
                legendTab{addtoLegend} = 'optogenetics';
                addedOpto = 1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        plot([0 max(numPoints)],percentCorrectChance,'k--')
       % plot([0 max(numPoints)],percentCorrectThreshold,'k--')
        ylim([0 1])
        xlim([0 length(numPoints)])
        ylabel('Performance');
        xlabel('Session Number');
    end
    
    
    %plot Dprime
    subplot(noSubPlots,1,4)
    addtoLegend = 1;
    addedWater = 0;
    addedNeg = 0;
    addedOpto = 0;
    
    legendAdd = [];
    legendTab = [];
    for jj = 1:length(numPoints);
        
        if waterSchedule(jj)==1;
            circle1 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(jj)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedWater==0)
                legendAdd(addtoLegend,:) = circle1;
                legendTab{addtoLegend} = 'water schedule';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        elseif negReinforcer(jj)==1;
            circle2 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(jj),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = circle2;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        else
            plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        end
        hold on
        
        if optogenetics(jj)==1;
            circle3 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'*y','MarkerSize', 5);
            hold on
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = circle3;
                legendTab{addtoLegend} = 'optogenetics';
                addedOpto = 1;
                addtoLegend=addtoLegend+1;
            end
        end

    %    plot([0 max(numPoints)],dprimeChance,'k--')
         plot([0 max(numPoints)],dprimeThreshold,'k--')
        xlim([0 length(numPoints)])
        ylabel('dprime');
        xlabel('Session Number');
    end
    
    
    if(length(legendAdd))>0
        hLL = legend([legendAdd(1:length(legendAdd))], legendTab{:});
        newPosition = [0 0 0.2 0.1];
        set(hLL, 'Position', newPosition, 'Box', 'off')
    end
    
    baseFileName = strcat(tempTitle); %save old figure
    saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\basicSessionProperties',baseFileName),'jpeg');
    
    %plot sessionType performance to get an idea of how inclinded to lick the
    %mouse is
    
    plotRows = 4;
    plotCols = 3;
    plotTally = (plotRows*plotCols);
    numFigs = 1;
    
    if (plotON==1)
        ff=figure;clf
        set(ff,'Position',positionGraph2);
    else
        figure('Visible','off');clf;
    end
    
    set(gcf,'name',tempTitle,'numbertitle','off')
    currPlot = 1;
    
    addtoLegend = 1;
    addedWater = 0;
    addedNeg = 0;
    addedOpto = 0;
    legendAdd = [];
    
    for h = 1:length(numPoints)
        saved = 0;
        activeAngles = cell2mat(basicPropertiesToPlot{h,1}.performance);
        probLick = cell2mat(basicPropertiesToPlot{h,1}.probLicking);
        plotAngles = possibleAngles-270;
        [~,c] = find(isnan(activeAngles));
        activeAngles(c) = [];
        plotAngles(c) = [];
        probLick(c) = [];
        trialTypecombo = [basicPropertiesToPlot{h,1}.HIT basicPropertiesToPlot{h,1}.MISS;basicPropertiesToPlot{h,1}.FA basicPropertiesToPlot{h,1}.CR];
        
        if(length(activeAngles)>2) %if more than 2 angles were presented
            subplot(plotRows,plotCols,currPlot);
            plot(plotAngles,activeAngles,'o-');
            hold on
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Performance');
            ylim([0 1])
            xlimit = xlim;
            text(xlimit(1),1.05,basicPropertiesToPlot{h,1}.namedata)
            str = strcat('d'' ', num2str(basicPropertiesToPlot{h,1}.dprime));
            text(xlimit(1),0.1,str)
            
            if(negReinforcer(h,1)==1)
                squareLeg1 =  plot(xlimit(2),0.1,'rs','MarkerFaceColor','r', 'MarkerSize',8);
                if (addedNeg==0)
                    legendAdd(addtoLegend,:) = squareLeg1;
                    legendTab{addtoLegend} = 'neg reinforcement';
                    addedNeg=1;
                    addtoLegend=addtoLegend+1;
                end
            end
            
            if(waterSchedule(h,1)==1)
                squareLeg2 = plot(xlimit(2),0.2,'gs','MarkerFaceColor','b', 'MarkerSize',8);
                
                if (addedOpto==0)
                    legendAdd(addtoLegend,:) = squareLeg2;
                    legendTab{addtoLegend} = 'water schedule';
                    addedOpto=1;
                    addtoLegend=addtoLegend+1;
                end
            end
            
            if(optogenetics(h,1)==1)
                squareLeg3 = plot(xlimit(2),0.3,'bs','MarkerFaceColor','g', 'MarkerSize',8);
                if (addedWater==0)
                    legendAdd(addtoLegend,:) = squareLeg3;
                    legendTab{addtoLegend} = 'optogenetics';
                    addedWater=1;
                    addtoLegend=addtoLegend+1;
                end
            end
            
            subplot(plotRows,plotCols,currPlot);
            plot(plotAngles,probLick,'o-');
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Licking Probability');
            ylim([0 1])
            
            subplot(plotRows,plotCols,currPlot)
            bar(trialTypecombo, 'stacked');
            hold on
            XlableAxis = {'HIT/MISS'; 'FA/CR'};
            set(gca,'XTickLabel',XlableAxis);
            currPlot=currPlot+1;
            
            if rem(currPlot,plotTally)==1 %if the value is divisable by 4 - open a new plot
                baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
                if(length(legendAdd))>0
                    hLL = legend([legendAdd(1:length(legendAdd))], legendTab{:});
                    newPosition = [0 0 0.2 0.1];
                    set(hLL, 'Position', newPosition, 'Box', 'off')
                    addtoLegend = 1;
                    addedWater = 0;
                    addedNeg = 0;
                    addedOpto = 0;
                    clear legendTab legendAdd;
                    
                end
                saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
                saved = 1;
                if (plotON==1)
                    ff=figure;clf
                    set(ff,'Position',positionGraph2);
                else
                    figure('Visible','off');clf;
                end
                numFigs = numFigs+1;
                currPlot = 1;
                legendAdd = [];
            end
        end
        
        if saved==0
            if currPlot>1
                baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
                saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
            end
        end
        
    end
    if  (plotON==1)
        if (currPlot==1)
            
            close(ff)
        end
    end
    
    plotRows = 4;
    plotCols = 2;
    plotTally = (plotRows*plotCols);
    numFigs = 1;
    
    if (plotON==1)
        fff=figure;clf
        set(fff,'Position',positionGraph3);
    else
        figure('Visible','off');clf;
    end
    
    set(gcf,'name',tempTitle,'numbertitle','off')
    currPlot = 1;
    tn=1;
    tempdPrime = nan(length(numPoints),2);
    tempdPerformance = nan(length(numPoints),2);
    
    for h = 1:length(numPoints)

        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            activeAnglesSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceSTIM);
            activeAnglesnoSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceNoSTIM);
            probLickSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingSTIM);
            probLickNoSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingNoSTIM);
            plotAngles = possibleAngles;
            [~,c] = find(isnan(activeAnglesSTIM));
            activeAnglesSTIM(c) = [];
            activeAnglesnoSTIM(c) = [];
            plotAngles(c) = [];
            probLickSTIM(c) = [];
            probLickNoSTIM(c) = [];
            
            subplot(plotRows,plotCols,currPlot);
            line1 = plot(plotAngles,activeAnglesSTIM,'o-r');
            hold on;
            line2 = plot(plotAngles,activeAnglesnoSTIM,'o-k');
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Performance');
            xlimit = xlim;
            ylim([0 1])
            text(xlimit(1),1.05,basicPropertiesToPlot{h,1}.namedata)
            str = ['d'' C/S',' ', num2str(basicPropertiesToPlot{h,1}.dprimeNoSTIM),' / ', num2str(basicPropertiesToPlot{h,1}.dprimeSTIM)];
            text(xlimit(1),0.1,str)
            
            tempdPrime(tn,:) = [basicPropertiesToPlot{h,1}.dprimeNoSTIM, basicPropertiesToPlot{h,1}.dprimeSTIM];
            tempdPerformance(tn,:) = [basicPropertiesToPlot{h,1}.sessionperformanceNoSTIM, basicPropertiesToPlot{h,1}.sessionperformanceSTIM];
            tn=tn+1;
            
            subplot(plotRows,plotCols,currPlot);
            plot(plotAngles,probLickSTIM,'o-r');
            hold on;
            plot(plotAngles,probLickNoSTIM,'o-k');
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Licking Probability');
            ylim([0 1])
        end
        
        if (currPlot>1)
            if rem(currPlot,plotTally)==1 %if the value is divisable by 4 - open a new plot
                baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
                hL = legend([line1, line2],{'Stimulated', 'Control'});
                newPosition = [0 0 0.2 0.1];
                set(hL, 'Position', newPosition, 'Box', 'off')
                
                saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
                if (plotON==1)
                    fff=figure;clf
                    set(fff,'Position',positionGraph3);
                else
                    figure('Visible','off');clf;
                end
                numFigs = numFigs+1;
                currPlot = 1;
            end
        end
    end
    
    %PLOT STIM V NONSTIM
    
    if (plotON==1)
        if currPlot==1
            close(fff)
        end
    end
    currPlot=1;
    
   tempdPrime =  tempdPrime(~any(isnan(tempdPrime),2),:);
   tempdPerformance =  tempdPerformance(~any(isnan(tempdPerformance),2),:);
    
    if any(tempdPrime)
        
        plotRows = 2;
        plotCols = 1;
        numFigs = 1;
        
        if (plotON==1)
            ffff=figure;clf
            set(ffff,'Position',positionGraph4);
        else
            figure('Visible','off');clf;
        end
        
        meandPrime = mean(tempdPrime);
        SEMdPrime = std(tempdPrime)/sqrt(length(tempdPrime));
        [~,p] = ttest(tempdPrime(:,1),tempdPrime(:,2));
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempdPrime','k-')
        hold on
        errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
        set(gca,'XTick',[1:2]);
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        set(gca, 'Xlim',[xlimit(1) xlimit(2)])
        ylimit = ylim;
        ylim([0 ylimit(2)]);
        set(gca, 'Ylim',[0 ylimit(2)])
        str = ['P = ', num2str(p)];
        ylabel('d prime');
        text(xlimit(1),0.2,str)
        currPlot=currPlot+1;
        
        
        meandPrime = nanmean(tempdPerformance);
        SEMdPrime = nanstd(tempdPerformance)/sqrt(length(tempdPerformance));
        [~,p] = ttest(tempdPerformance(:,1),tempdPerformance(:,2));
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempdPerformance','k-')
        hold on
        errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
        set(gca,'XTick',[1:2]);
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        set(gca, 'Xlim',[xlimit(1) xlimit(2)])
        ylimit = ylim;
        ylim([0 ylimit(2)]);
        set(gca, 'Ylim',[0 ylimit(2)])
        str = ['P = ', num2str(p)];
        ylabel('Performance');
        text(xlimit(1),0.2,str)
        
    end
    
    trimtally=1;
    indTrim = nan(length(numPoints),1);
    for h = 1:length(numPoints)
        if ~isempty(basicPropertiesToPlot{h,1}.trimming)
            indTrim(h,1)=1;
        else
            indTrim(h,1)=0;
        end
    end
    
    %take the two last sessions before trimming for baseline
    firstIndTrim = find(indTrim,1);
    
    tempTrimPerformance = nan(length(numPoints),1);
    tempTrimdPrime = nan(length(numPoints),1);
    tempTrimWhiskers = cell(length(numPoints),1);
    
    
    if ~isempty(firstIndTrim)
        for h=1:length(numPoints)
            
            if h==firstIndTrim-1 || h==firstIndTrim-2
                tempTrimPerformance(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
                tempTrimdPrime(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
                tempTrimWhiskers{trimtally,1} = 'Full';
                trimtally = trimtally+1;
                
            elseif indTrim(h,1)==1
                
                tempTrimPerformance(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
                tempTrimdPrime(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
                tempTrimWhiskers{trimtally,1} = basicPropertiesToPlot{h,1}.trimming;
                trimtally = trimtally+1;
                
            end
        end
        
        
        
        plotRows = 2;
        plotCols = 1;
        plotTally = (plotRows*plotCols);
        numFigs = 1;
        
        if (plotON==1)
            ffff=figure;clf
            set(ffff,'Position',positionGraph5);
        else
            figure('Visible','off');clf;
        end
        
        set(gcf,'name',tempTitle,'numbertitle','off')
        currPlot = 1;
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempTrimPerformance,'o-k','LineWidth',5);
        hold on
        set(gca,'XTick',[1:length(tempTrimPerformance)]);
        set(gca,'XTickLabel',tempTrimWhiskers);
        ylim([0 1]);
        set(gca, 'Ylim',[0 1])
        ylabel('Performance');
        xlabel('whisker trimming');
        currPlot=currPlot+1;
        

        plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)
       %         plot([1 max(xlim)],percentCorrectThreshold,'k--','LineWidth',2)


        
        subplot(plotRows,plotCols,currPlot);
        plot(tempTrimdPrime,'o-k','LineWidth',5);
        hold on
        set(gca,'XTick',[1:length(tempTrimdPrime)]);
        set(gca,'XTickLabel',tempTrimWhiskers);
            %  plot([1 max(xlim)],dprimeChance,'k--','LineWidth',2)
                plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)

        ylimit = ylim;
        
        if ylimit(1)<0;
            
            yLimitSet = ylimit(1);
            
        else
            
            yLimitSet = 0;
        end
        ylim([yLimitSet ylimit(2)]);
        set(gca, 'Ylim',[yLimitSet max(tempTrimdPrime)])
        ylabel('d prime');
        xlabel('whisker trimming');
        currPlot=currPlot+1;
        % baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
        
        saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
        
    end
    
    
    plotRows = 2;
    plotCols = 1;
    numFigs = 1;
    currPlot=1;
    Segtally=1;

    if (plotON==1)
        fV = figure;clf
        set(fV,'Position',positionGraph7);
    else
        figure('Visible','off');clf;
    end
    
    set(gcf,'name',tempTitle,'numbertitle','off');
    
    indSeg = nan(length(numPoints),1);
    for h = 1:length(numPoints)
        if ~isempty(basicPropertiesToPlot{h,1}.segmentPerformance)
            indSeg(h,1)=1;
            Segtally = Segtally+1;
        else
            indSeg(h,1)=0;
        end
    end
    
    

    for h = 1:length(numPoints)
        if (indSeg(h,1)==1)
            subplot(plotRows,plotCols,currPlot);
            plot(basicPropertiesToPlot{h,1}.segmentPerformance,'ok-','LineWidth',5)
            hold on
            
        plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)
         %       line([1 max(xlim)],percentCorrectThreshold,'k--','LineWidth',2)
            % SessionTypes = {'Cntr' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %SessionTypes = {'Block 1' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %NumTicks = length(SessionTypes);
            NumTicks = length(basicPropertiesToPlot{h,1}.segmentPerformance);
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks));
            currPlot = currPlot+1;
            
            ylim([0 1]);
            title('Performance over blocks','fontWeight','bold');
            xlabel('Block number','fontWeight','bold');
            ylabel('Percent Correct ((hit+CR)/totalTrials', 'fontWeight','bold');
            %set(gca,'XTickLabel',SessionTypes)
            subplot(plotRows,plotCols,currPlot);
            plot(basicPropertiesToPlot{h,1}.segmentdPrime,'ok-','LineWidth',5)
            hold on
                %    line([1 max(xlim)],dprimeChance,'k--','LineWidth',2)
                plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
            % SessionTypes = {'Cntr' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %SessionTypes = {'Block 1' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %NumTicks = length(SessionTypes);
            NumTicks = length(basicPropertiesToPlot{h,1}.segmentdPrime);
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks));
            currPlot = currPlot+1;
            ylimit = ylim;

            if ylimit(1)<0;
                yLimitSet = ylimit(1);
            else
                yLimitSet = 0;
            end
            
            ylim([yLimitSet ylimit(2)]);
            set(gca, 'Ylim',[yLimitSet max(tempTrimdPrime)])
            title('Performance over blocks','fontWeight','bold');
            xlabel('Block number','fontWeight','bold');
            ylabel('d Prime', 'fontWeight','bold');
            text(1,ylimit(1)+0.1,basicPropertiesToPlot{h,1}.namedata)
            %set(gca,'XTickLabel',SessionTypes)
             baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
             saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\segmented',baseFileName),'jpeg');
numFigs = numFigs+1;
 currPlot=1;
  if (plotON==1)
        fV = figure;clf
        set(fV,'Position',positionGraph7);
    else
        figure('Visible','off');clf;
    end
        else
            indSeg(h,1)=0;
        end
    end
   % saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\segmented',tempTitle),'jpeg');
    
   
    if (plotON==1)
        if currPlot==1
            close(fV)
        end
    end
end
%
% XTickLabel = get(gca,'XTickLabel')
% hxLabel = get(gca,'XLabel')
% set(hxLabel,'Units','data')
% xLabelPosition = get(hxLabel,'Position')
% XTick = get(gca,'XTick')
% y = xLabelPosition(2)
%  y = repmat(y,length(XTick),1)
%  fs = get(gca,'fontsize')
%  hText = text(XTick,y,XTickLabel,'fontsize',fs)
%  set(hText,'Rotation',90,'HorizontalAlignment','right')


end
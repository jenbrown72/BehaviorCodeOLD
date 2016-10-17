function [basicPropertiesToPlot, possibleAngles,info] = JB_basicBehaviorProperties(plotON,cleanUp)

% JB_basicBehaviorProperties once in the Mouse ID.data file loads DATA.mat
% and run a basic analysis on each session
%

%   [] = JB_basicBehaviorProperties(plotON,cleanUp)
%   plotON = 1 (plots output), = 0 (no plot generated)
%   cleanUp = 1 (excluded sessions/data based on set criteria), = 0,
%   included every session and data point
%   [basicPropertiesToPlot] returned one session per day
%   [possibleAngles] returned vector of angles
%   [info] returns basic info from each session


% Examples:
%   [basicPropertiesToPlot, possibleAngles,info] = JB_basicBehaviorProperties(1,0);
%
%   [basicPropertiesToPlot] = JB_basicBehaviorProperties(0,0);
%

positionGraph1 = [1555 293 343 702];

if (exist('plotON','var') && plotON==1)
    f = figure;clf
    set(f,'Position',positionGraph1);
else
    plotON=0;
    figure('Visible','off');clf;
end

if (exist('cleanUp','var') && cleanUp==1)
else
    cleanUp=0;
end

% analysis stop early criteria
minNoTrialsTotal = 200; % no of trials required to have been completed for this criteria to be used
stopExThreshold = 32; % no of consequtive misses and correct rejections to say the mouse wasn't engaged and stop analysis (4xS8)

stepsPerResolution = 360*2; %current code looks for up and down changes
wheelDistancePerTurn = 47.75;
cmPerStep = wheelDistancePerTurn/stepsPerResolution;

load('DATA.mat');
disp(' ')
disp(DATA.mouseID);
possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];
basicProperties = cell(100,1);

%define which txt file column the different digital inputs are read in from
%processing
encoder0Pos = 1;
targetArea = 2;
reward = 3;
licks = 4;
LEDstate = 5;
count = 6;
millis = 7;
angle = 11;
catchTrial = 12;
trialType = 13;
timeDiff = 14;
encoderStepDiff = 15;
licksR = 16;
rewardR = 17;

for i=1:length(DATA.allFiles);
            clear startTrial endTrial 
            
            %delete any concaternated rows that contain all 1s
                    segIdx = find(strcmp('Segments = ',  DATA.allFiles{i}.metaData));
                    if isempty(segIdx) 
                        DATA.allFiles{i}.rawData(all(DATA.allFiles{i}.rawData==1,2),:)=[];
                    end
                %clean up and delete when catchTrial stayed on permanantely (big in
    %early Sep 2016)
    
    
 
    [idxCatchStart, ~] = find(diff(DATA.allFiles{i}.rawData(:,catchTrial))>0);
    [idxCatchStop, ~] = find(diff(DATA.allFiles{i}.rawData(:,catchTrial))<0);
    
    if (length(idxCatchStart))==(length(idxCatchStop))
    elseif isempty(idxCatchStart)
             continue
        else

        DATA.allFiles{i}.rawData(idxCatchStart(end):end,:)=[];
    end
    
   %clean up and delete when angle is 0 (start)
    noAngle = find(DATA.allFiles{i}.rawData(:,angle)==0);
    for ik = 1:length(noAngle)
        DATA.allFiles{i}.rawData(noAngle(ik),:)=nan;
    end
    DATA.allFiles{i}.rawData = DATA.allFiles{i}.rawData(~any(isnan(DATA.allFiles{i}.rawData),2),:);
    
     %clean up and delete if last trial wasn't completed
    [startTrial, ~] = find(diff(DATA.allFiles{i}.rawData(:,targetArea))>0);
    [endTrial, ~] = find(diff(DATA.allFiles{i}.rawData(:,targetArea))<0);
    
    if ~isempty(startTrial)
            if (length(startTrial))==(length(endTrial))
    else
        DATA.allFiles{i}.rawData(startTrial(end):end,:)=[];
        startTrial(end)=[];
            end
            
    else
        continue
    end
    
    if isempty( DATA.allFiles{i}.rawData)
        continue
    else
        % Delete any blips where the clock went backwards
        negTime = find(diff(DATA.allFiles{i}.rawData(:,millis))<=0)+1;
        for ik = 1:length(negTime)
            DATA.allFiles{i}.rawData(negTime(ik),:)=nan;
        end
        DATA.allFiles{i}.rawData = DATA.allFiles{i}.rawData(~any(isnan(DATA.allFiles{i}.rawData),2),:);
        
        %Extract basic informaiton from the data file
        basicProperties{i,1}.date = DATA.allFiles{i}.date;
        basicProperties{i,1}.datenum = DATA.allFiles{i}.dateFromFile;
        basicProperties{i,1}.sessionDuration = ((DATA.allFiles{i}.rawData(end,millis)-DATA.allFiles{i}.rawData(1,millis))/1000)/60; % Calculate duration of session in minutes
        basicProperties{i,1}.fileLength = length(DATA.allFiles{i}.rawData);
        basicProperties{i,1}.noTrials = sum(diff(DATA.allFiles{i}.rawData(:,count))>0); % DATA.allFiles{i}.rawData(end,count)-DATA.allFiles{i}.rawData(1,count);
        basicProperties{i,1}.trimType = NaN;
        basicProperties{i,1}.pairsDprime = nan(length(possibleAngles),1);
        basicProperties{i,1}.pairsDprimeSTIM = nan(length(possibleAngles),1);
        basicProperties{i,1}.pairsDprimeNoSTIM = nan(length(possibleAngles),1);
        
        tempName = DATA.allFiles{i}.name(strfind(DATA.allFiles{i}.name,'_201'):strfind(DATA.allFiles{i}.name,'_Box'));
        tempLocation = strfind(tempName,'_');
        for hh = 1:length(tempLocation)-1;
            if(tempLocation(hh+1)-tempLocation(hh))==2;
                tempName = strcat( tempName(1:tempLocation(hh)),'0',tempName(tempLocation(hh)+1:end));
                tempLocation = tempLocation+1;
            end
        end
        
        toDelete = strfind(tempName,'_');
        tempName(toDelete)='.';
        basicProperties{i,1}.namedata = tempName;
        
        %get running velocity
        threshold = 50000;
        aboveThreshold = find(DATA.allFiles{i}.rawData(:,encoder0Pos)>threshold);
        for j = 1:length(aboveThreshold)
            DATA.allFiles{i}.rawData(aboveThreshold(j),:)=nan;
        end
        basicProperties{i,1}.totalStepsPerMin = nansum(diff(DATA.allFiles{i}.rawData(:,encoder0Pos))>0)/basicProperties{i,1}.sessionDuration; % Calculate running velocity
        basicProperties{i,1}.avgVelocity = (nansum(DATA.allFiles{i}.rawData(:,encoderStepDiff))*cmPerStep)/(basicProperties{i,1}.sessionDuration*60); % Calculate running velocity in cm per sec
        
        %Extract information from the metaData file
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
        
        %find out if this session has positive or neg GO stim
        orgGOIdx = find(strcmp('orgGOangleTab = ', metaData));
        
        if isempty(orgGOIdx)
            basicProperties{i,1}.orgGOstim = 1;
        else
            orgGOstim = metaData{orgGOIdx,2};
            if ischar(orgGOstim)==1
                basicProperties{i,1}.orgGOstim = str2num(orgGOstim);
            else
                basicProperties{i,1}.orgGOstim = orgGOstim;
            end
        end
        
        
        %find out if this session was a trimming session
        trimIdx = find(strcmp('Trimming = ', metaData));
        if ~isempty(trimIdx)
            tempTrim = metaData{trimIdx,2};
            basicProperties{i,1}.trimming = tempTrim;
            trimPattern = {basicProperties{i,1}.trimming};
            
            countTrim = strfind(trimPattern,'C');
            if length(countTrim{1})==4; %this is a row data)
                basicProperties{i,1}.trimType = 'Row';
            end
        else
            basicProperties{i,1}.trimming = [];
        end
        
        %find out if this session was a 2AFC
        AFCIdx = find(strcmp('2ACF = ', metaData));
        if isempty(AFCIdx)
            basicProperties{i,1}.afc = 0;
        else
            afcStatus = metaData{AFCIdx,2};
            if ischar(afcStatus)==1
                basicProperties{i,1}.afc = str2num(afcStatus);
            else
                basicProperties{i,1}.afc = afcStatus;
            end
        end
        
        %find out if this session had segments
        segIdx = find(strcmp('Segments = ', metaData));
        if ~isempty(segIdx)
            tempSeg = metaData{segIdx,2};
            tempSegNum = str2num(tempSeg);
            if (tempSegNum(1)==1) && (tempSegNum(2)==1)
                tempData = DATA.allFiles{i}.rawData;
                tempData = tempData(~any(isnan(tempData),2),:);
                index = find(all(tempData==1,2));
                index = [index; length(tempData)]; %add a one to find the start of the first segment
                %             segments(1,1) = 1; %add first trial
                for kk=1:length(index)-1;
                    if ((tempData(index(kk)-1,6)==(tempData(index(kk)+1,6))))
                        countCummulative = 0;
                        diffCount = tempData(index(kk)-1,6);
                    else
                        countCummulative = tempData(index(kk)-1,6);
                        diffCount = tempData(index(kk)-1,6);
                    end
                    tempData(index(kk)+1:index(kk+1)-1,6) = tempData(index(kk)+1:index(kk+1)-1,6)+countCummulative;
                    basicProperties{i,1}.segments(kk,2) = diffCount;
                    basicProperties{i,1}.segments(kk+1,1) = diffCount+1;
                end
                basicProperties{i,1}.segments(length(basicProperties{i,1}.segments),2) = tempData(length(tempData),6); %add last trial
            else
                basicProperties{i,1}.segments = str2num(tempSeg);
            end
        else
            basicProperties{i,1}.segments = [];
        end
        
        %%New
        
        %get index for start of each trial
        basicProperties{i,1}.temptrialTypesCatch = NaN;
        basicProperties{i,1}.performanceCatch=NaN;
        basicProperties{i,1}.dPrimeCatch = NaN;
        tempAngleTypes=NaN;
        tempCatchTrials = NaN;
        tempOptogenetics=NaN;
        
        
        for jj = 1:length(startTrial)
        tempAngleTypes(jj,1) = DATA.allFiles{i}.rawData(startTrial(jj,1)+2,angle);
        tempCatchTrials(jj,1) = DATA.allFiles{i}.rawData(startTrial(jj,1),catchTrial);
        tempOptogenetics(jj,1) = DATA.allFiles{i}.rawData(startTrial(jj,1),LEDstate);
        end
        endTrialCatch = [];
        %first target area index of catch trials
        
        %see if mouse licked
        clear temptrialTypes;
        clear temptrialTypesCatch;
        
        
        if ~isempty(endTrial)
            for hh=1:length(endTrial)
                lick = find((DATA.allFiles{i}.rawData(startTrial(hh,1):endTrial(hh,1),licks))==1);
                if isempty(lick);
                    licksTrial(hh,1)=0;
                else
                    licksTrial(hh,1)=1;
                end
                if ((tempAngleTypes(hh,1)<270) &&  (basicProperties{i,1}.orgGOstim==1)) %origional Go trial
                    if (licksTrial(hh,1)==1) %HIT
                        temptrialTypes(hh,1) = 1;
                    elseif (licksTrial(hh,1)==0) %MISS
                        temptrialTypes(hh,1) = 3;
                    end
                elseif ((tempAngleTypes(hh,1)>270) &&  (basicProperties{i,1}.orgGOstim==0)) %switched NoGo trial
                    if (licksTrial(hh,1)==1) %False Alarm
                        temptrialTypes(hh,1) = 2;
                    elseif (licksTrial(hh,1)==0) %Correct Rejection
                        temptrialTypes(hh,1) = 4;
                    end
                elseif ((tempAngleTypes(hh,1)>270) &&  (basicProperties{i,1}.orgGOstim==1))%origional NoGo trial
                    if (licksTrial(hh,1)==1) %False Alarm
                        temptrialTypes(hh,1) = 2;
                    elseif (licksTrial(hh,1)==0) %Correct Rejection
                        temptrialTypes(hh,1) = 4;
                    end
                elseif ((tempAngleTypes(hh,1)>270) &&  (basicProperties{i,1}.orgGOstim==0))%switched Go trial
                    if (licksTrial(hh,1)==1) %HIT
                        temptrialTypes(hh,1) = 1;
                    elseif (licksTrial(hh,1)==0) %MISS
                        temptrialTypes(hh,1) = 3;
                    end
                end
            end
        end
        
        basicProperties{i,1}.temptrialTypes = temptrialTypes(find(tempCatchTrials==0));
        basicProperties{i,1}.temptrialTypesCatch = temptrialTypes(find(tempCatchTrials==1));
        
        tempTypesTrials = temptrialTypes(find(tempCatchTrials==0));
        tempTypesCatch = temptrialTypes(find(tempCatchTrials==1));
        
        tempAnglesTrials = tempAngleTypes(find(tempCatchTrials==0));
        tempAnglesCatch = tempAngleTypes(find(tempCatchTrials==1));
        
        
        tempOptogeneticsTrials = tempOptogenetics(find(tempCatchTrials==0));
        tempOptogeneticsCatch = tempOptogenetics(find(tempCatchTrials==1));
        
        %calculate performance on catch trials
        
        [Hit, FA, Miss,CR, performance, ~, dprime] = JB_countTrialTypes(tempTypesCatch);
        basicProperties{i,1}.performanceCatch = performance;
        basicProperties{i,1}.dPrimeCatch = dprime;

        %Look at raw data to get stats on performance
        %     [startTrial, ~] = find(diff(DATA.allFiles{i}.rawData(:,trialType))>0);
        %     temptrialTypes = DATA.allFiles{i}.rawData(startTrial+2,trialType);
        
        %Calculate Block Performance
        blockBinSize = 8;
        %  basicProperties{i,1}.temptrialTypes = temptrialTypes;
        %     blockPerformance = temptrialTypes;
        %     blockPerformance(blockPerformance==4)=1;
        %     blockPerformance(blockPerformance==1)=1;
        %     blockPerformance(blockPerformance==2)=2;
        %     blockPerformance(blockPerformance==3)=2;
        %     binspacing = 1:blockBinSize:length(blockPerformance);
        %       basicProperties{i,1}.blockPerformancePlot=[];
        %     for s = 1:length(binspacing)-1;
        %          basicProperties{i,1}.blockPerformancePlot(s) = length(find(blockPerformance(binspacing(s):binspacing(s+1)-1)==1))/blockBinSize;
        %     end
        
        %plot the session performance
        if ~isempty(tempTypesTrials)
            if (plotON==1)
                figure(f);clf
                tempTitle = basicProperties{i,1}.namedata;
                tempTitle(findstr(tempTitle,'_'))=[];
                set(gcf,'name',tempTitle,'numbertitle','off');
            end
            
            for kk=1:length(tempTypesTrials)
                if tempTypesTrials(kk)==1; % Hit
                    line([0 1], [kk kk], 'Color', 'g','LineWidth',2);
                    hold on;
                elseif tempTypesTrials(kk)==3; % No Lick on GO stimulus = Miss
                    line([0 1], [kk kk], 'Color', 'r','LineWidth',2);
                    hold on;
                elseif tempTypesTrials(kk)==4; % No Lick on Distractor stimulus = Correct Rejection
                    line([2 3], [kk kk], 'Color', 'g','LineWidth',2);
                    hold on;
                elseif tempTypesTrials(kk)==2; % Lick on Distractor stimulus = False Alarm
                    line([2 3], [kk kk], 'Color', 'r','LineWidth',2);
                    hold on;
                end
            end
            
            xlim([-0.5 3.5])
            %   ylim([1 length(tempTypesTrials)])
            set(gca,'YDir','reverse');
            
            %find how many consequtive trials with no lick there were
            ff = (tempTypesTrials(minNoTrialsTotal:end,1)==3) | (tempTypesTrials(minNoTrialsTotal:end,1)==4); %miss and correct rejection
            tempCum = 1;
            if ~isempty(ff)
                cumCount = nan(length(ff),1);
                for v = 1:length(ff)-1
                    if (ff(v) && ff(v+1)) ==1;
                        tempCum = tempCum+1;
                    else
                        tempCum = 1;
                    end
                    cumCount(v,1) = tempCum;
                end
                idxStop = (find(cumCount>stopExThreshold))+(minNoTrialsTotal);
                if ~isempty(idxStop) && (cleanUp == 1)
                    stopIdx = idxStop(1)-stopExThreshold;
                    if stopIdx>0
                        X = ['Date ', basicProperties{i,1}.date(1:11), ' Stopping Session Analysis Early By', num2str(length(startTrial)- stopIdx), 'Trials', 'On trial No.', num2str(stopIdx),'Out Of',num2str(length(startTrial))] ;
                        line(xlim, [stopIdx stopIdx], 'Color', 'k','LineWidth',2);
                        startTrial(stopIdx:end)=[];
                        disp(X)
                    end
                end
            end
            if (plotON==1)
                waitforbuttonpress;
            end
        end
 
        %if this was a segmented session, calculate performance in each session
        if ~isempty( basicProperties{i,1}.segments);
            sessionDivisions = basicProperties{i,1}.segments;
            for j = 1:length(sessionDivisions)
                if sessionDivisions(j,1)==0;
                    idxStart = 1;
                else
                    idxStart = sessionDivisions(j,1);
                end
                idxEnd = sessionDivisions(j,2);
                tempSegmentData = tempTypesTrials(idxStart:idxEnd,:);
                [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(tempSegmentData);
                basicProperties{i,1}.segmentPerformance(j,1) =  performance;
                basicProperties{i,1}.segmentdPrime(j,1) = dprime;
            end
        else
            basicProperties{i,1}.segmentPerformance = [];
            basicProperties{i,1}.segmentdPrime = [];
        end
        
        %Calculare performance for each angle presented and overall
        %performance
        
        basicProperties{i,1}.noTrialsSTIM = nan(length(possibleAngles),1);
        basicProperties{i,1}.noTrialsNoSTIM = nan(length(possibleAngles),1);
        
        for kk=1:length(possibleAngles)
            stimCount = 1;
            nostimCount = 1;
            idxAngles = find(tempAnglesTrials==possibleAngles(kk));
            performanceAngleTemp = tempTypesTrials(idxAngles);
            [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTemp);
            basicProperties{i,1}.performance{kk} =  performance;
            basicProperties{i,1}.probLicking{kk} =  licking;
            basicProperties{i,1}.trialTypesAngle{:,kk} = performanceAngleTemp;
            basicProperties{i,1}.performanceSTIM{kk} = nan;
            basicProperties{i,1}.probLickingSTIM{kk} = nan;
            %         basicProperties{i,1}.trialTypesAngleSTIM{kk} = nan;
            basicProperties{i,1}.performanceNoSTIM{kk} = nan;
            basicProperties{i,1}.probLickingNoSTIM{kk} = nan;
            %         basicProperties{i,1}.trialTypesAngleNoSTIM{kk} = nan;
            
            if ~isempty(idxAngles) %if this angle was used
                performanceAngleTempSTIM = [];
                performanceAngleTempNoSTIM = [];
                if basicProperties{i,1}.optogenetics==1 %if optogenetics was used calculate performance with and without stim
                    for ff = 1:length(idxAngles)
                        if (tempOptogeneticsTrials(idxAngles(ff))==1)
                            performanceAngleTempSTIM(stimCount,1) = tempTypesTrials(idxAngles(ff));
                            stimCount = stimCount+1;
                        else
                            performanceAngleTempNoSTIM(nostimCount,1) = tempTypesTrials(idxAngles(ff));
                            nostimCount = nostimCount+1;
                        end
                    end
                    [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTempSTIM);
                    basicProperties{i,1}.trialTypesAngleSTIM{:,kk} = performanceAngleTempSTIM;
                    basicProperties{i,1}.performanceSTIM{kk} = performance;
                    basicProperties{i,1}.probLickingSTIM{kk} =  licking;
                    basicProperties{i,1}.noTrialsSTIM(kk) =  length(performanceAngleTempSTIM);
                    
                    [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTempNoSTIM);
                    basicProperties{i,1}.trialTypesAngleNoSTIM{:,kk} = performanceAngleTempNoSTIM;
                    basicProperties{i,1}.performanceNoSTIM{kk} =  performance;
                    basicProperties{i,1}.probLickingNoSTIM{kk} =  licking;
                    basicProperties{i,1}.noTrialsNoSTIM(kk) =  length(performanceAngleTempNoSTIM);
                end
                
            else
                
                basicProperties{i,1}.trialTypesAngleSTIM{:,kk} = [];
                basicProperties{i,1}.trialTypesAngleNoSTIM{:,kk} = [];
            end
        end
        

        %Calculate how many of each trial type there were
        [Hit, FA, Miss,CR, performance, ~, dprime] = JB_countTrialTypes(tempTypesTrials);
        basicProperties{i,1}.HIT = Hit; %HitTrialCount
        basicProperties{i,1}.FA = FA; % FATrialCount
        basicProperties{i,1}.MISS = Miss; %MissTrialCount
        basicProperties{i,1}.CR = CR; %CRTrialCount
        basicProperties{i,1}.sessionperformance = performance;
        basicProperties{i,1}.dprime=dprime;
        
        %now look at optogenetic sessions
        if (basicProperties{i,1}.optogenetics==1)
            basicProperties{i,1}.HITStim = sum((tempTypesTrials==1)&(tempOptogeneticsTrials==1))+1; %HitTrialCount
            basicProperties{i,1}.FAStim = sum((tempTypesTrials==2)&(tempOptogeneticsTrials==1))+1; % FATrialCount
            basicProperties{i,1}.MISSStim = sum((tempTypesTrials==3)&(tempOptogeneticsTrials==1))+1; %MissTrialCount
            basicProperties{i,1}.CRStim = sum((tempTypesTrials==4)&(tempOptogeneticsTrials==1))+1; %CRTrialCount
            basicProperties{i,1}.HITNoStim = sum((tempTypesTrials==1)&(tempOptogeneticsTrials==0))+1; %HitTrialCount
            basicProperties{i,1}.FANoStim = sum((tempTypesTrials==2)&(tempOptogeneticsTrials==0))+1; % FATrialCount
            basicProperties{i,1}.MISSNoStim = sum((tempTypesTrials==3)&(tempOptogeneticsTrials==0))+1; %MissTrialCount
            basicProperties{i,1}.CRNoStim = sum((tempTypesTrials==4)&(tempOptogeneticsTrials==0))+1; %CRTrialCount
            basicProperties{i,1}.dprimeSTIM=norminv(( basicProperties{i,1}.HITStim/( basicProperties{i,1}.HITStim+basicProperties{i,1}.MISSStim)),0,1)-norminv((basicProperties{i,1}.FAStim/(basicProperties{i,1}.FAStim+basicProperties{i,1}.CRStim)),0,1);
            basicProperties{i,1}.dprimeNoSTIM=norminv(( basicProperties{i,1}.HITNoStim/( basicProperties{i,1}.HITNoStim+basicProperties{i,1}.MISSNoStim)),0,1)-norminv((basicProperties{i,1}.FANoStim/(basicProperties{i,1}.FANoStim+basicProperties{i,1}.CRNoStim)),0,1);
            basicProperties{i,1}.sessionperformanceSTIM = (basicProperties{i,1}.HITStim+basicProperties{i,1}.CRStim)/(sum(tempOptogeneticsTrials==1)+4);
            basicProperties{i,1}.sessionperformanceNoSTIM = (basicProperties{i,1}.HITNoStim+basicProperties{i,1}.CRNoStim)/(sum(tempOptogeneticsTrials==0)+4);
        end
        
        %import as metaDATA cell array
        for j=1:length(metaData)
            if(strcmp('Auto reward = ', metaData{j,1}));
                autoSession = metaData{j,2};
            end
        end
        
        match = strcmp('Orientation Selected = ', metaData(:,1));
        clear ind
        ind = nan(length(match),1);
        for kk=1:length(match)
            if match(kk,1)==1;
                if ischar(metaData{kk,2})
                    ind(kk,1) = str2num(metaData{kk,2});
                else
                    ind(kk,1) =(metaData{kk,2});
                end
            end
        end
        tempAngles = ind(ind>0);
        
        %Define session type
        tempsessionType = strcat('S',num2str(length(tempAngles)));
        basicProperties{i,1}.stimNumber = length(tempAngles);
        if autoSession=='1'
            tempsessionType = strcat(tempsessionType,'auto');
            basicProperties{i,1}.stimNumber = 0;
        end
        basicProperties{i,1}.mouseID = DATA.mouseID;
        basicProperties{i,1}.sessionType = tempsessionType;
        
        %Detect Lick Latency
        targetStartIdx = find(diff(DATA.allFiles{i}.rawData(:,targetArea))>0);
        targetStopIdx = find(diff(DATA.allFiles{i}.rawData(:,targetArea))<0);
        
        % make sure the first start index is greater than the  - if not, delete
        % the first stop startTrial.
        
        % added
        if (isempty(targetStartIdx) || isempty(targetStopIdx))
            
            targetStartIdx = [];
            targetStopIdx = [];
        end
        
        % added
        
        if ~isempty(targetStopIdx)
            if (targetStopIdx(1)<targetStartIdx(1))
                targetStopIdx(1) = [];
            end
            
            if (length(targetStartIdx)<2)
                basicProperties{i,1}.frequencyLicks = nan;
                basicProperties{i,1}.firstLickLatency = nan;
            else
                
                for kk = 1:length(targetStartIdx)-1
                    tempDiffLicks = diff(DATA.allFiles{i}.rawData(targetStartIdx(kk):targetStopIdx(kk),licks));
                    tempAngleLick = DATA.allFiles{i}.rawData(targetStartIdx(kk),angle);
                    tempStimLick = DATA.allFiles{i}.rawData(targetStartIdx(kk)+10,LEDstate);
                    basicProperties{i,1}.frequencyLicks(kk,:) = [sum(tempDiffLicks>0) tempAngleLick tempStimLick];
                    if ~isempty(find(tempDiffLicks>0,1));
                        firstLickidx = (find((tempDiffLicks)>0,1)+targetStartIdx(kk));
                        basicProperties{i,1}.firstLickLatency(kk,:) = [DATA.allFiles{i}.rawData(firstLickidx,millis)-DATA.allFiles{i}.rawData(targetStartIdx(kk),millis) tempAngleLick tempStimLick];
                    else
                        basicProperties{i,1}.firstLickLatency(kk,:) = [nan tempAngleLick tempStimLick];
                    end
                end
                
                %devide up into latency to each angle.
                
                uniqueAngles = unique(basicProperties{i,1}.firstLickLatency(:,2));
                delA = ~ismember(uniqueAngles,possibleAngles);
                uniqueAngles(delA)=[];
                
                for gg = 1:length(uniqueAngles)
                    %  basicProperties{i,1}.angleLickLatency{gg} = nan(length(basicProperties{i,1}.firstLickLatency),2);
                    idxAngleSTIM = (basicProperties{i,1}.firstLickLatency(:,2)==uniqueAngles(gg)) & (basicProperties{i,1}.firstLickLatency(:,3)==1);    %stimulated tri
                    idxAngleNoSTIM = (basicProperties{i,1}.firstLickLatency(:,2)==uniqueAngles(gg)) & (basicProperties{i,1}.firstLickLatency(:,3)==0);
                    basicProperties{i,1}.angleLickLatency{gg,1} = basicProperties{i,1}.firstLickLatency(idxAngleNoSTIM,1);
                    basicProperties{i,1}.angleLickLatency{gg,2} = basicProperties{i,1}.firstLickLatency(idxAngleSTIM,1);
                    basicProperties{i,1}.angleLick(gg,1) = uniqueAngles(gg);
                end
            end
            
            
        else
            basicProperties{i,1}.frequencyLicks = nan;
            basicProperties{i,1}.firstLickLatency = nan;
        end
        
        %         %Find angle pairs
        
        for d = 1:(floor(length(possibleAngles)/2));
            pair1 = basicProperties{i,1}.trialTypesAngle{1,(1+d-1)};
            pair2 = basicProperties{i,1}.trialTypesAngle{1,end-(d-1)};
            if (~isempty(pair1) && ~isempty(pair2))
                %  if (~isnan(pair1) || ~isnan(pair2))
                %find which was the go and nogo
                if any(pair1==1) || any(pair1==3) %pair 1 go stim
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
                else
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
                end
                basicProperties{i,1}.pairsDprime(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
            else
                basicProperties{i,1}.pairsDprime(d,1) =  NaN;
            end
            pair1 = basicProperties{i,1}.trialTypesAngleSTIM{1,(1+d-1)};
            pair2 = basicProperties{i,1}.trialTypesAngleSTIM{1,end-(d-1)};
            %          if (~isnan(pair1) || ~isnan(pair2))
            if (~isempty(pair1) && ~isempty(pair2))
                %find which was the go and nogo
                if any(pair1==1) || any(pair1==3) %pair 1 go stim
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
                else
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
                end
                basicProperties{i,1}.pairsDprimeSTIM(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
            else
                basicProperties{i,1}.pairsDprimeSTIM(d,1) =  NaN;
            end
            pair1 = basicProperties{i,1}.trialTypesAngleNoSTIM{1,(1+d-1)};
            pair2 = basicProperties{i,1}.trialTypesAngleNoSTIM{1,end-(d-1)};
            %         if (~isnan(pair1) || ~isnan(pair2))
            if (~isempty(pair1) && ~isempty(pair2))
                %find which was the go and nogo
                if any(pair1==1) || any(pair1==3) %pair 1 go stim
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
                else
                    [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                    [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
                end
                basicProperties{i,1}.pairsDprimeNoSTIM(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
            else
                basicProperties{i,1}.pairsDprimeNoSTIM(d,1) =  NaN;
            end
            basicProperties{i,1}.pairsDiff(d,1) = possibleAngles(end-(d-1))-possibleAngles(d);
        end
        
%         basicProperties{i,1}.pairsDprimeNoSTIM(isnan(basicProperties{i,1}.pairsDprimeNoSTIM(:,1)),:)=[];
%         basicProperties{i,1}.pairsDprimeSTIM(isnan(basicProperties{i,1}.pairsDprimeSTIM(:,1)),:)=[];
%         basicProperties{i,1}.pairsDprime(isnan(basicProperties{i,1}.pairsDprime(:,1)),:)=[];
        
    end
end

%delete empty cells
basicProperties(any(cellfun(@isempty,basicProperties),2),:)=[];

%Now go along and extra only a single session type from each day, e.g. if
%both S2 and S8 was run.. we look to see which file is the longest, and
%then use that
basicPropertiesToPlot = cell(100,1);
days = cell(length(basicProperties),1);
for k = 1:length(basicProperties),
    days{k,1} = basicProperties{k,1}.datenum(1:8);
end
plotNumber = 1;
[~,~,startTrial] = unique(days(:,1));
startTrial = sort(startTrial);
unique_idx = accumarray(startTrial(:),(1:length(startTrial))',[],@(x) {sort(x)});

for k = 1:length(unique_idx)
    tempRunVal=nan(length(unique_idx),1);
    if length(unique_idx{k})==1
        if (basicProperties{unique_idx{k}}.noTrials>minNoTrialsTotal || (cleanUp == 0))
            basicPropertiesToPlot{plotNumber,1} = basicProperties{unique_idx{k}};
            X = ['Date ', basicPropertiesToPlot{plotNumber,1}.namedata, ' ', basicPropertiesToPlot{plotNumber,1}.sessionType, ' Performance ', num2str( basicPropertiesToPlot{plotNumber,1}.sessionperformance), ' dPrime ', num2str(basicPropertiesToPlot{plotNumber,1}.dprime), ' #trials ',num2str(basicPropertiesToPlot{plotNumber,1}.noTrials),' Duration ' ,num2str(basicPropertiesToPlot{plotNumber,1}.sessionDuration)];
            disp(X)
            info(k,:) = {basicPropertiesToPlot{plotNumber,1}.namedata ,num2str(basicPropertiesToPlot{plotNumber,1}.dprime)};
            plotNumber = plotNumber+1;
            
        else
            X = ['EXCLUDED: Not enough Trials', basicProperties{k,1}.date(1:11),' #trials ',num2str(basicProperties{k,1}.noTrials)]
            disp(X)
        end
        
    elseif length(unique_idx{k})>1
        t = unique_idx{k};
        for kk = 1:length(t)
            tempRunVal(kk,1) = basicProperties{t(kk)}.noTrials;
        end
        [~,di] = max(tempRunVal); %find the data set where the mouse ran the most
        
        if (basicProperties{t(di)}.noTrials>minNoTrialsTotal || (cleanUp == 0))
            basicPropertiesToPlot{plotNumber,1} = basicProperties{t(di)};
            X = ['Date ', basicPropertiesToPlot{plotNumber,1}.namedata, ' ', basicPropertiesToPlot{plotNumber,1}.sessionType, ' Performance ', num2str( basicPropertiesToPlot{plotNumber,1}.sessionperformance), ' dPrime ', num2str(basicPropertiesToPlot{plotNumber,1}.dprime), ' #trials ',num2str(basicPropertiesToPlot{plotNumber,1}.noTrials),' Duration ' ,num2str(basicPropertiesToPlot{plotNumber,1}.sessionDuration)];
            disp(X)
            plotNumber = plotNumber+1;
        else
            X = ['EXCLUDED: Not enough Trials', ' ','Date ', basicProperties{t(di)}.date(1:11),' #trials ',num2str(basicProperties{t(di)}.noTrials)];
            disp(X)
        end
        
    end
end

%delete empty cells
basicPropertiesToPlot(any(cellfun(@isempty,basicPropertiesToPlot),2),:)=[];

% basicPropertiesToPlot = basicPropertiesToPlot';
%%

JB_plotPerformance(basicPropertiesToPlot,1)
% JB_plotOptogenetics(basicPropertiesToPlot,1)
JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,1)
% %
% % % [DATAavg] = JB_plotSelectionPerformance(basicPropertiesToPlot,possibleAngles,plotON,0,1,1);
% % % %JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,plotON)
% JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,1)
% JB_plotTrimming(basicPropertiesToPlot, plotON)
% JB_plotSegments(basicPropertiesToPlot,plotON)
end
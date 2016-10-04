function [ data ] = JB_plotRunProperties(DATA,sessionNo)
%UNTITLED Summary of this function goes here
%%   Detailed explanation goes here

tempDATA = DATA.allFiles{1,sessionNo}.rawData;
metaData = DATA.allFiles{1,sessionNo}.metaData;
positionGraph1 = [5 558 1847 420];
plotColor = {'DodgerBlue','Red'};
outcome = {'Hit & Correct Rejection';'Miss & False Alarm'};


% find out if this was a Go <270 or Go>270
orgGOIdx = find(strcmp('orgGOangleTab = ', metaData));
if isempty(orgGOIdx)
    orgGOstim = str2num(metaData{orgGOIdx,2});
else
    orgGOstim = str2num(metaData{orgGOIdx,2});
end

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
trialType = 13;
sampleRateSpacing = 1; %1ms sampling

possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];

trialTypes = {'Hit';'FalseAlarm';'Miss';'CorrectRejection'};
trialTypesCode = [1 2 3 4];

threshold = 50000;
aboveThreshold = find(tempDATA(:,1)>threshold);

% if mouse ran backwards before starting the trial, the step number will be
% 65540 - delete these points

for i = 1:length(aboveThreshold)
    tempDATA(aboveThreshold(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

%Delete any blips in the recordings - missaligned txt input - find these by
%searching for non binary values in the columns where these are expected
[ind, ~] = find(tempDATA(:,[2,3,4,5,8])>1); %all these inputs should be binary
if ~isempty(ind) %if there was a blip, delete this row
    UniqInd = unique(ind);
    for n=1:length(UniqInd);
        tempDATA(UniqInd(n),:)=nan;
    end
    tempDATA = tempDATA(~any(isnan(tempDATA),2),:);
end

% Delete any blips where the clock went backwards
negTime = find(diff(tempDATA(:,7))<=0) ;
for i = 1:length(negTime)
    tempDATA(negTime(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

% Delete any blips where the clock didn't increase
negTime = find(diff(tempDATA(:,7))==0);
for i = 1:length(negTime)
    tempDATA(negTime(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

C = unique(tempDATA(:,angle));
data.angleOrientations = C(C>0);
tempangles = data.angleOrientations;
%% Divide up data into single revolutions/trials

revolutionsDATA = find((diff(tempDATA(:,1)))<-100); % find where encoder is reset back to 0 to mark the end of a trial/start of enw trial
for j=1:length(revolutionsDATA)-1; %Dont include first revolution as wont be properly callibrated.
    data.trial{j}.raw(:,:) = tempDATA(revolutionsDATA(j)+1:revolutionsDATA(j+1),:);
end

% Check Thats Detection of stimulus worked
toDelete = [];
p=1;
for j = 1:length(data.trial);
    [loc2,~] = find((data.trial{j}.raw(:,targetArea)>0)); %did the target enter the target region - should do on all trials, if not - delete
    if isempty(loc2) || loc2(1)==1;
        toDelete(p) = j;
        p=p+1;
    end
end

if ~isempty(toDelete)
    data.trial(toDelete)=[];
end

%%
binSize = 20;
PreStim=2000; %time stamps are in 1ms bins
PostStim=2000;
TimePointsMin = -PreStim:1:PostStim-1;

baselinePeriod = 500/binSize; %500ms for baseline

toPlot = nan(length(data.trial),4);
performance = nan(length(data.trial),2);
acceraltionProps = nan(length(data.trial),4);
inflectionPoint = nan(length(data.trial),1);

for j = 1: length(data.trial);
    % make sample times uniform and zeroed
    clear tempFiltered;
    SampleTimes = data.trial{j}.raw(:,millis)-data.trial{j}.raw(1,millis);
    
    %make each data point uniformly sampled
    for i = 1:size((data.trial{j}.raw),2);
        tempTrial = data.trial{j}.raw(:,i);
        filteredSampleTimes = 0:sampleRateSpacing:max(SampleTimes); %1msec spacing
        tempFiltered(:,i)=interp1(SampleTimes,tempTrial,filteredSampleTimes);
    end
    
    % make sure all values are integers
    data.trial{j}.filteredDATA = [tempFiltered(:,1) (round(tempFiltered(:,2:end)))];
    
    tempRun = diff(data.trial{j}.filteredDATA(:,1));
    startIdx = find(diff(data.trial{j}.filteredDATA(:,2))==1);
    addPreStim = [];
    addPostStim = [];
    if (startIdx-PreStim<1)
        addPreStim=NaN*ones(1,PreStim - (startIdx-1)); %%add nans to start of vector
        currPreStim =  PreStim-(PreStim - startIdx+1);
    else
        currPreStim = PreStim;
    end
    
    if length(tempRun)<startIdx+PostStim-1
        addPostStim=NaN*ones(1,startIdx+PostStim-1-length(tempRun)); %%add nans to start of vector
        currPostStim =  PostStim-((startIdx+PostStim-1)-length(tempRun));
    else
        currPostStim = PostStim;
    end
    
    tempVelocityStim = tempRun(startIdx-currPreStim:startIdx+currPostStim-1,1)';
    tempRunA = [addPreStim tempVelocityStim addPostStim];
    tempTimeA = TimePointsMin;
    
    [running_speed,time_bins] = JB_computeSpeed(tempRunA,sampleRateSpacing,tempTimeA,1,binSize);
    running_speed = smooth(running_speed);
    
    acceleration = diff(running_speed)./diff(time_bins');
    [gg] = interp1(time_bins(1:end-1),acceleration,time_bins(1:end-1));
    [pMin,peakMin] = min(gg(10:end-10));
    [pMax,peakMax] = max(gg(10:end-10));
    
    baselineVel = nanmean(running_speed(1:baselinePeriod,1));
    acceraltionProps(j,:) = [pMin time_bins(peakMin) pMax time_bins(peakMax)];
        run_speed = running_speed';
    gg_subset = gg(1,peakMin:end);
    time_subset = time_bins(1,peakMin:end-1);
        running_speed_subset = run_speed(1,peakMin:end-1);
          tallyC=1;
    clear idxIn;
    idxIn = NaN;
    
          for kk = 1:length(gg_subset)-1
                if ((gg_subset(kk)<0) && (gg_subset(kk+1)>0))
                    idxIn(tallyC,1) = kk;
                    tallyC = tallyC+1;
                end
            end
%     idxIn = dsearchn(gg_subset',0);
if ~isnan(idxIn)
    inflectionPoint(j,:) = time_subset(idxIn(1)); 
else
        inflectionPoint(j,:) = NaN;
end 
%     
%     idxIn = dsearchn(gg_subset',0);
%     inflectionPoint(j,:) = time_subset(idxIn);
    
    %smooth data and find slopes
    %low pass filter acceleration data
    smoothed_gg = smooth(gg);
    slope_gg = diff(smoothed_gg)/diff(time_bins(1:end-1)');
    %max deceleration
    [min_slope min_slopeIdx] = min(slope_gg(10:end-10,1));
   
    if ~isnan(idxIn)
    toPlot(j,1) =  time_subset(idxIn(1));
    else 
        toPlot(j,1) = nan;
    end
    toPlot(j,2) = time_bins(peakMin);
    toPlot(j,3) = time_bins(peakMax);
    toPlot(j,4) = baselineVel-(running_speed(peakMin));
    toPlot(j,5) = time_bins(min_slopeIdx);
    
    performance(j,1) = data.trial{j}.filteredDATA(1,angle);
    performance(j,2) = data.trial{j}.raw(find(data.trial{j}.raw(:,trialType)>0,1),trialType);
%     
%     figure(77);clf
%     subplot(3,1,1)
%     plot(time_bins,running_speed)
%     hold on
%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),running_speed(peakMin),'*g')
%     plot(time_bins(peakMax),running_speed(peakMax),'*r')
%     ylabel('Velocity')
%     
%     subplot(3,1,2)
%     plot(time_bins(1:end-1),gg)
%     hold on
%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),gg(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),gg(peakMin),'*g')
%     plot(time_bins(peakMax),gg(peakMax),'*r')
%     ylabel('Acceleration')
%     
%     subplot(3,1,3)
%     plot(time_bins(1:end-2),slope_gg)
%     hold on
%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),slope_gg(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),slope_gg(peakMin),'*g')
%     plot(time_bins(peakMax),slope_gg(peakMax),'*r')
%     ylabel('Slope')
%     
        
    figure(78);
    
    subplot(2,length(tempangles)/2,find(tempangles==data.trial{j}.filteredDATA(1,angle)))
    if ((performance(j,2)==1) || performance(j,2)==4);
    plot(time_bins,running_speed,'Color',rgb(plotColor{1}))
    else
        plot(time_bins,running_speed,'Color',rgb(plotColor{2}))
    end
    hold on

%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),running_speed(peakMin),'*g')
%     plot(time_bins(peakMax),running_speed(peakMax),'*r')
    ylabel('Velocity')
    xlab = num2str(data.trial{j}.filteredDATA(1,angle)-270);
 xlabel(strcat('angle',xlab))
%         subplot(length(tempangles),1,find(tempangles==data.trial{j}.filteredDATA(1,angle)))
%     if ((performance(j,2)==1) || performance(j,2)==4);
%     plot(time_bins(1:end-1),gg,'Color',rgb(plotColor{1}))
%     else
%         plot(time_bins(1:end-1),gg,'Color',rgb(plotColor{2}))
%     end
%     hold on
%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),gg(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),gg(peakMin),'*g')
%     plot(time_bins(peakMax),gg(peakMax),'*r')
%     ylabel('Acceleration')
%     
%             subplot(length(tempangles),1,find(tempangles==data.trial{j}.filteredDATA(1,angle)))
%     if ((performance(j,2)==1) || performance(j,2)==4);
%     plot(time_bins(1:end-2),slope_gg,'Color',rgb(plotColor{1}))
%     else
%         plot(time_bins(1:end-2),slope_gg,'Color',rgb(plotColor{2}))
%     end
%     plot(time_bins(1:end-2),slope_gg)
%     hold on
%     plot([inflectionPoint(j,:) inflectionPoint(j,:)],[min(ylim) max(ylim)],'b--','LineWidth',2)
%     plot(time_bins(min_slopeIdx),slope_gg(min_slopeIdx),'*m')
%     plot(time_bins(peakMin),slope_gg(peakMin),'*g')
%     plot(time_bins(peakMax),slope_gg(peakMax),'*r')
%     ylabel('Slope')
%     
    
    %         jjj = waitforbuttonpress;
end
    figure(78);
hLL = legend(outcome)
newPosition = [0 0 0.1 0.1];
set(hLL, 'Position', newPosition, 'Box', 'off')

%%
clear accelerationTemp meanData semData
% varLabels = {'inflectionPoint';'minAcceleration';'MaxAcceleration'}
varLabels = {'inflection Point';'Min Acceleration Time';'Max Acceleration Time';'peak Vel Change';'min slope'};
outcome = {'Hit & Correct Rejection';'Miss & False Alarm'};
outcomePoints = [1 4;2 3];%HIT1,FA2,MISS3,CR4,
% outcomePoints = [1 1;2 2];%HIT1,FA2,MISS3,CR4,
figure(103)

% tempangles(tempangles<270)=[];

accelerationTemp = nan(100,length(tempangles));
plotNo = 1;

for gi = 1:length(outcome)
    clear accelerationTemp
    for gh = 1:length(varLabels);
        for gg = 1:length(tempangles) 
            tally = 1;
            %only for NoGo trials
            idx =  find(performance(:,1)==(tempangles(gg)));
            for ggh = 1:length(idx)
                if ((performance(idx(ggh),2)==outcomePoints(gi,1)) || (performance(idx(ggh),2)==outcomePoints(gi,2)))
                    accelerationTemp(tally,gg) = toPlot(idx(ggh),gh);
                    tally = tally+1;
                end
            end
        end
        
        subplot(length(outcome),length(varLabels),plotNo)
        [meanData,~,~,semData] = JB_calBasicStats(accelerationTemp);
        plot(tempangles,meanData,'o-k','MarkerSize',3,'LineWidth',3)
        hold on
        errorbar(tempangles,meanData,semData,'k','Linestyle','none');
        labels = num2str(tempangles-270);
        
        set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','Ydir','reverse')
        view(-90,90)
        xlabel('angles');
        ylabel(varLabels(gh));
        title(outcome{gi});
        plotNo = plotNo+1;
    end
end

figure(114);clf
set(gcf,'Position',positionGraph1)
tempangles = data.angleOrientations;
% tempangles(tempangles<270)=[];
accelerationTemp = nan(100,length(tempangles));

for gi = 1:length(outcome)
      clear accelerationTemp
    plotNo=1;
    for gh = 1:length(varLabels);
        for gg = 1:length(tempangles)
            tally = 1;
            %only for NoGo trials
            idx =  find(performance(:,1)==(tempangles(gg)));
            for ggh = 1:length(idx)
                if ((performance(idx(ggh),2)==outcomePoints(gi,1)) || (performance(idx(ggh),2)==outcomePoints(gi,2)))
                    accelerationTemp(tally,gg) = toPlot(idx(ggh),gh);
                    tally = tally+1;
                end
            end
        end
        
        subplot(1,length(varLabels),plotNo)
        [meanData,~,~,semData] = JB_calBasicStats(accelerationTemp);
        hold on
        errorbar(tempangles,meanData,semData,'o-','Color',rgb(plotColor{gi}),'MarkerSize',3,'LineWidth',3);
        labels = num2str(tempangles-270);
        
        set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','Ydir','reverse')
        view(-90,90)
        xlabel('angles');
        ylimit = ylim;
        if ((gi==2) && (plotNo==length(varLabels)))
            if (orgGOstim==1)
                text(tempangles(1),ylimit(2)+(ylimit(2)/5),'GO');
                text(tempangles(end),ylimit(2)+(ylimit(2)/5),'NOGO');
            else
                text(tempangles(1),ylimit(2),'NOGO');
                text(tempangles(end),ylimit(2),'GO');
            end
        end
        ylabel(varLabels(gh));
        plotNo = plotNo+1;
    end
end
hLL = legend(outcome);
newPosition = [0 0 0.1 0.1];
set(hLL, 'Position', newPosition, 'Box', 'off')

figure(25)
plot(acceraltionProps(:,2),acceraltionProps(:,1),'*g')
hold on
plot(acceraltionProps(:,4),acceraltionProps(:,3),'*r')
plot(inflectionPoint,0,'ok')




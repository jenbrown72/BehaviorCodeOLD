function [ output_args ] = JB_trialPerformance(tempDATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load DATA.mat
positionGraph2 = [680 288 314 690];
positionGraph1 = [1501 71 384 902];
tempDATA = DATA.allFiles{1,45}.rawData;

midPoint=270;
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

%% Divide up data into single revolutions/trials

revolutionsDATA = find((diff(tempDATA(:,1)))<-100); % find where encoder is reset back to 0 to mark the end of a trial/start of enw trial

for j=1:length(revolutionsDATA)-1; %Dont include first revolution as wont be properly callibrated.
    data.trial{j,1}.raw(:,:) = tempDATA(revolutionsDATA(j)+1:revolutionsDATA(j+1),:);
end

% Check Thats Detection of stimulus worked
toDelete = [];
p=1;

for j = 1:length(data.trial);
    [loc2,~] = find((data.trial{j,1}.raw(:,2)>0)); %did the target enter the target region - should do on all trials, if not - delete
    if isempty(loc2) || loc2(1)==1;
        toDelete(p) = j;
        p=p+1;
    end
end

if ~isempty(toDelete)
    data.trial(toDelete)=[];
end

%% Define Target Time, Position, Correct (Hit, correct Rejection) or Incorrect (Miss, or False Alarm) data.trial

sampleRateSpacing = 2;
data.SessionPerformance.midPoint = midPoint;
data.SessionPerformance.MissCount = 0;
data.SessionPerformance.HitCount = 0;
data.SessionPerformance.CorrectRejectionCount = 0;
data.SessionPerformance.FalseAlarmCount = 0;
data.SessionPerformance.midPointLickCount = 0;
data.SessionPerformance.midPointNoLickCount = 0;

data.SessionPerformance.MissCountStim = 0;
data.SessionPerformance.MissCountNoStim = 0;
data.SessionPerformance.HitCountStim = 0;
data.SessionPerformance.HitCountNoStim = 0;
data.SessionPerformance.CorrectRejectionCountStim = 0;
data.SessionPerformance.CorrectRejectionCountNoStim = 0;
data.SessionPerformance.FalseAlarmCountStim = 0;
data.SessionPerformance.FalseAlarmCountNoStim = 0;
data.SessionPerformance.midPointLickCountStim = 0;
data.SessionPerformance.midPointLickCountNoStim = 0;
data.SessionPerformance.midPointNoLickCountStim = 0;
data.SessionPerformance.midPointNoLickCountNoStim = 0;

for j = 1: length(data.trial);
    
    % make sample times uniform and zeroed
    SampleTimes = data.trial{j,1}.raw(:,7)-data.trial{j,1}.raw(1,7);
    clear tempFiltered;
    
    for i = 1:size((data.trial{j,1}.raw),2);
        tempTrial = data.trial{j,1}.raw(:,i);
        filteredSampleTimes = 0:sampleRateSpacing:max(SampleTimes); %10msec spacing
        tempFiltered(:,i)=interp1(SampleTimes,tempTrial,filteredSampleTimes);
    end
    
    % make sure all values are integers
    data.trial{j,1}.filteredDATA = [tempFiltered(:,1) (round(tempFiltered(:,2:end)))];
    A = data.trial{j,1}.filteredDATA(:,4);
    % make sure licks are represented as just one tick mark
    B = diff(data.trial{j,1}.filteredDATA(:,4))==1; %%make sure licks are represented as just one tick mark
    B(numel(A)) = 0;
    
    data.trial{j,1}.filteredDATA(:,4)=B;
    data.trial{j,1}.timeStamps(:,1) = [0;(cumsum(diff(data.trial{j,1}.filteredDATA(:,7)))/1000)];
    data.trial{j,1}.GracePeriod = ((data.trial{j,1}.filteredDATA(1,5)));
    data.trial{j,1}.Performance.StimType = data.trial{j,1}.filteredDATA(1,11);
    data.trial{j,1}.Performance.TrialType = data.trial{j,1}.raw(find(data.trial{j,1}.raw(:,13)>0),13);
    
    [loc2,~] = find((data.trial{j,1}.filteredDATA(:,2)>0)); % Identify when stim entered TargetArea
    idx3 = find((data.trial{j,1}.filteredDATA(:,3)==1)); % Identify RewardTime if Hit data.trial
    
    % If Reward was given, update information on record properties
    if ~isempty(idx3); %
        data.trial{j,1}.Stimulus.RewardIdx = idx3(1);
    else
        data.trial{j,1}.Stimulus.RewardIdx = nan;
    end
    
    %Gather information about stimulus presentation
    data.trial{j,1}.Stimulus.TargetMin = data.trial{j,1}.filteredDATA(loc2(1),1); % stim entered TargetArea Step
    data.trial{j,1}.Stimulus.TargetMinIdx = loc2(1); % stim entered TargetArea Idx
    data.trial{j,1}.Stimulus.TargetMinTime = data.trial{j,1}.timeStamps(loc2(1),1);% stim entered TargetArea time
    
    data.trial{j,1}.Stimulus.TargetMax = data.trial{j,1}.filteredDATA(loc2(1),end); % timeout TargetArea Step
    data.trial{j,1}.Stimulus.TargetMaxIdx = loc2(end); % timeout TargetArea Idx
    data.trial{j,1}.Stimulus.TargetMaxTime = data.trial{j,1}.timeStamps(loc2(1),end);% timeout TargetArea time
    
    %Gather information about Answer Period
    data.trial{j,1}.Stimulus.AnswerPeriodDuration = data.trial{j,1}.filteredDATA(loc2(1),12);
    data.trial{j,1}.Stimulus.AnswerPeriodStartIdx = loc2(1);
    data.trial{j,1}.Stimulus.AnswerPeriodStartTime = data.trial{j,1}.timeStamps(loc2(1));
    data.trial{j,1}.Stimulus.AnswerPeriodEndIdx = loc2(end);
    data.trial{j,1}.Stimulus.AnswerPeriodEndTime = data.trial{j,1}.timeStamps(loc2(end));
    
    % Determine if animal licked during Answer Period
    
    %Look at time of first lick for all stim and all trial types cover both
    %the answer period and 1 second after this
    SampleLickTime = 2000/sampleRateSpacing;  % = 2seconds
    EndTime = data.trial{j,1}.Stimulus.AnswerPeriodStartIdx+SampleLickTime;
    if EndTime>length(data.trial{j,1}.filteredDATA)
        EndTime = length(data.trial{j,1}.filteredDATA);
    end
    
    % get lick times aroound the 2 seconds after answer period started
    LickIdx = (data.trial{j,1}.Stimulus.AnswerPeriodStartIdx+(find((data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx:EndTime,4)==1),1)));
    % get time of first lick
    data.trial{j,1}.Performance.LickTime = data.trial{j,1}.timeStamps(LickIdx)-data.trial{j,1}.timeStamps(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx);
    
    % if no lick within 2 seconds - give a value of nan
    if isempty(data.trial{j,1}.Performance.LickTime)
        data.trial{j,1}.Performance.LickTime=nan;
    end
    
    %if reward port opened find licks before this happened
    if ~isnan(data.trial{j,1}.Stimulus.RewardIdx) %If there was a reward given
        data.trial{j,1}.Performance.LickFrequency = sum(data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx:data.trial{j,1}.Stimulus.RewardIdx,4)==1)/(data.trial{j,1}.Stimulus.AnswerPeriodDuration/1000);
        tempFirstLickIdx = (data.trial{j,1}.Stimulus.AnswerPeriodStartIdx+(find((data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx:data.trial{j,1}.Stimulus.RewardIdx,4)==1),1)));
        tempLickLatency = data.trial{j,1}.timeStamps(tempFirstLickIdx)-data.trial{j,1}.timeStamps(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx);
        if ~isempty(tempLickLatency)
            data.trial{j,1}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j,1}.Performance.FirstLickLatency = nan;
        end
    else
        if length(data.trial{j,1}.filteredDATA)>data.trial{j,1}.Stimulus.AnswerPeriodEndIdx;%if reward port didn't open find out if any licks occured during whloe answer period
            maxTimeIdx = data.trial{j,1}.Stimulus.AnswerPeriodEndIdx;
        else
            maxTimeIdx = length(data.trial{j,1}.filteredDATA);
        end
        data.trial{j,1}.Performance.LickFrequency = sum(data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx:maxTimeIdx,4)==1)/(data.trial{j,1}.Stimulus.AnswerPeriodDuration/1000);
        tempFirstLickIdx = (data.trial{j,1}.Stimulus.AnswerPeriodStartIdx+(find((data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx:maxTimeIdx,4)==1),1)));
        tempLickLatency = data.trial{j,1}.timeStamps(tempFirstLickIdx)-data.trial{j,1}.timeStamps(data.trial{j,1}.Stimulus.AnswerPeriodStartIdx);
        if ~isempty(tempLickLatency)
            data.trial{j,1}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j,1}.Performance.FirstLickLatency = nan;
        end
    end
    
    if any(data.trial{j,1}.raw(:,5),1);
        data.trial{j,1}.Performance.StimTrial = 1;
    else
        data.trial{j,1}.Performance.StimTrial = 0;
    end
    
    % Update Trial Type Results
    if data.trial{j,1}.Performance.LickFrequency>0; % If there was licking..
        if (data.trial{j,1}.Performance.StimType < midPoint); % If this was the go stimulus - Hit trial
            data.SessionPerformance.HitCount = data.SessionPerformance.HitCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 1; % Anticipatory lick - Hit trial
            data.trial{j,1}.Performance.OutcomeCode = 1;
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial - test to see if column contains a 1 =  i.e. stim on
                data.SessionPerformance.HitCountStim = data.SessionPerformance.HitCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.HitCountNoStim = data.SessionPerformance.HitCountNoStim+1;
            end
        elseif (data.trial{j,1}.Performance.StimType > midPoint) % If this was a distractor stimulus
            data.SessionPerformance.FalseAlarmCount = data.SessionPerformance.FalseAlarmCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 1;% False Alarm - incorrect lick
            data.trial{j,1}.Performance.OutcomeCode = 2;
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial
                data.SessionPerformance.FalseAlarmCountStim = data.SessionPerformance.FalseAlarmCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.FalseAlarmCountNoStim = data.SessionPerformance.FalseAlarmCountNoStim+1;
            end
        elseif (data.trial{j,1}.Performance.StimType == midPoint) % If this was exactly in the middle
            data.SessionPerformance.midPointLickCount = data.SessionPerformance.midPointLickCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 1;% licked for 50% verticle
            data.trial{j,1}.Performance.OutcomeCode = 5; % 5 if licked, 6 if didn't
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial
                data.SessionPerformance.midPointLickCountStim = data.SessionPerformance.midPointLickCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.midPointLickCountNoStim = data.SessionPerformance.midPointLickCountNoStim+1;
            end
        end
    else
        if (data.trial{j,1}.Performance.StimType < midPoint); % If this was the no go stimulus
            data.SessionPerformance.MissCount = data.SessionPerformance.MissCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 0;  % No anticipatory lick - Miss trial
            data.trial{j,1}.Performance.OutcomeCode = 3;
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial
                data.SessionPerformance.MissCountStim = data.SessionPerformance.MissCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.MissCountNoStim = data.SessionPerformance.MissCountNoStim+1;
            end
        elseif (data.trial{j,1}.Performance.StimType > midPoint) % If this was a distractor stimulus
            data.SessionPerformance.CorrectRejectionCount = data.SessionPerformance.CorrectRejectionCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 0; % Correctly withheld lick
            data.trial{j,1}.Performance.OutcomeCode = 4;
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial
                data.SessionPerformance.CorrectRejectionCountStim = data.SessionPerformance.CorrectRejectionCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.CorrectRejectionCountNoStim = data.SessionPerformance.CorrectRejectionCountNoStim+1;
            end
        elseif (data.trial{j,1}.Performance.StimType == midPoint) % If this was exactly in the middle
            data.SessionPerformance.midPointNoLickCount = data.SessionPerformance.midPointNoLickCount+1;
            data.trial{j,1}.Performance.OutcomeLick = 0;% licked for 50% verticle
            data.trial{j,1}.Performance.OutcomeCode = 6; % 5 if licked, 6 if didn't
            if any(data.trial{j,1}.raw(:,5),1); %Photostimulation Trial
                data.SessionPerformance.midPointNoLickCountStim = data.SessionPerformance.midPointNoLickCountStim+1;
            else %No stimulation Trial
                data.SessionPerformance.midPointNoLickCountNoStim = data.SessionPerformance.midPointNoLickCountNoStim+1;
            end
        end
    end
    data.SessionPerformance.CumulativePerformance(j,1) = ((data.SessionPerformance.HitCount+data.SessionPerformance.CorrectRejectionCount)/j);
end

disp('________________Performance Lick Criteria_________________________');
disp(['HitCount ' num2str(data.SessionPerformance.HitCount)]);
disp(['MissCount ' num2str(data.SessionPerformance.MissCount)]);
disp(['CorrectRejectionCount ' num2str(data.SessionPerformance.CorrectRejectionCount)]);
disp(['FalseAlarmCount ' num2str(data.SessionPerformance.FalseAlarmCount)]);
disp(['midPointLickCount ' num2str(data.SessionPerformance.midPointLickCount)]);
disp(['midPointNoLickCount ' num2str(data.SessionPerformance.midPointNoLickCount)]);

if (data.SessionPerformance.HitCountStim+data.SessionPerformance.MissCountStim) % optogenetic sesison
    data.SessionPerformance.optogenetics = 1;
else
    data.SessionPerformance.optogenetics = 0;
end

if data.SessionPerformance.HitCountStim>0
    disp('________________Performance Stimulation_________________________');
    disp(['HitCountStim ' num2str(data.SessionPerformance.HitCountStim)]);
    disp(['HitCountNoStim ' num2str(data.SessionPerformance.HitCountNoStim)]);
    disp(['MissCountStim ' num2str(data.SessionPerformance.MissCountStim)]);
    disp(['MissCountNoStim ' num2str(data.SessionPerformance.MissCountNoStim)]);
    disp(['CorrectRejectionCountStim ' num2str(data.SessionPerformance.CorrectRejectionCountStim)]);
    disp(['CorrectRejectionCountNoStim ' num2str(data.SessionPerformance.CorrectRejectionCountNoStim)]);
    disp(['FalseAlarmCountStim ' num2str(data.SessionPerformance.FalseAlarmCountStim)]);
    disp(['FalseAlarmCountNoStim ' num2str(data.SessionPerformance.FalseAlarmCountNoStim)]);
    disp(['midPointLickCountStim ' num2str(data.SessionPerformance.midPointLickCountStim)]);
    disp(['midPointLickCountNoStim ' num2str(data.SessionPerformance.midPointLickCountNoStim)]);
    disp(['midPointNoLickCountStim ' num2str(data.SessionPerformance.midPointNoLickCountStim)]);
    disp(['midPointNoLickCountNoStim ' num2str(data.SessionPerformance.midPointNoLickCountNoStim)]);
end

data.SessionPerformance.totalSessions = sum(data.SessionPerformance.HitCount+ data.SessionPerformance.MissCount+data.SessionPerformance.CorrectRejectionCount+data.SessionPerformance.FalseAlarmCount);
data.CumulativePerformance(1:length(data.trial),1:5)=nan;

for j=1:length(data.trial);
    data.CumulativePerformance(j,1) = data.trial{j,1}.Performance.StimType;
    data.CumulativePerformance(j,2) = data.trial{j,1}.Performance.OutcomeLick;
    data.CumulativePerformance(j,3) = data.trial{j,1}.Performance.LickFrequency;
    data.CumulativePerformance(j,4) = data.trial{j,1}.Performance.FirstLickLatency;
    data.CumulativePerformance(j,5) = data.trial{j,1}.Performance.OutcomeCode;
    data.CumulativePerformance(j,6) = data.trial{j,1}.Performance.StimTrial;
end

% For all stim - find running tally for each orietnation
%identify all stimulate orientations used in the session
C = unique(data.CumulativePerformance(:,1));
data.SessionPerformance.orientations = C;

%create empty matrix for incoming data
if (length(C)>1)
    data.SessionPerformance.licklatency = nan(100,length(data.SessionPerformance.orientations)/2);
    data.SessionPerformance.licklatencyStim = nan(100,length(data.SessionPerformance.orientations)/2);
    data.SessionPerformance.licklatencyNoStim = nan(100,length(data.SessionPerformance.orientations)/2);
    data.SessionPerformance.licklatencyAngle = nan(1,length(data.SessionPerformance.orientations)/2);
end
data.SessionPerformance.lickTimeALL = nan(100,length(data.SessionPerformance.orientations));
data.SessionPerformance.lickTimeALLStim = nan(100,length(data.SessionPerformance.orientations));
data.SessionPerformance.lickTimeALLNoStim = nan(100,length(data.SessionPerformance.orientations));
data.SessionPerformance.lickTimeALLAngle = nan(100,length(data.SessionPerformance.orientations));
data.SessionPerformance.lickCountALL = nan(6,length(data.SessionPerformance.orientations));

for k = 1:length(C)
    tally = 1;
    tallyStim=1;
    tallyNoStim=1;
    emptyTally = 1;
    emptyTallyStim=1;
    emptyTallyNoStim=1;
    for j=1:length(data.trial); % Go Stimulus
        if data.trial{j,1}.Performance.StimType==C(k)
            data.SessionPerformance.lickTimeALLAngle(tally,k) = C(k);
            data.SessionPerformance.lickTimeALL(tally,k) = data.trial{j,1}.Performance.LickTime;
            tally = tally+1;
            if isnan(data.trial{j,1}.Performance.LickTime)
                emptyTally = emptyTally+1;
            end
            if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                data.SessionPerformance.lickTimeALLStim(tallyStim, k) = data.trial{j,1}.Performance.LickTime;
                tallyStim=tallyStim+1;
                if isnan(data.trial{j,1}.Performance.LickTime)
                    emptyTallyStim = emptyTallyStim+1;
                end
            else
                data.SessionPerformance.lickTimeALLNoStim(tallyNoStim, k) = data.trial{j,1}.Performance.LickTime;
                tallyNoStim=tallyNoStim+1;
                if isnan(data.trial{j,1}.Performance.LickTime)
                    emptyTallyNoStim = emptyTallyNoStim+1;
                end
            end
        end
    end
    data.SessionPerformance.lickCountALL(1,k) = tally;
    data.SessionPerformance.lickCountALL(2,k) = emptyTally;
    data.SessionPerformance.lickCountALL(3,k) = tallyStim;
    data.SessionPerformance.lickCountALL(4,k) = emptyTallyStim;
    data.SessionPerformance.lickCountALL(5,k) = tallyNoStim;
    data.SessionPerformance.lickCountALL(6,k) = emptyTallyNoStim;
    data.SessionPerformance.lickCountALLNorm(1,k) = emptyTally/tally;
    data.SessionPerformance.lickCountALLNorm(2,k) = emptyTallyStim/tallyStim;
    data.SessionPerformance.lickCountALLNorm(3,k) = emptyTallyNoStim/tallyNoStim;
end

for k = 1:length(C) % create empty field for stim
    
    if C(k)<midPoint;
        
        lickfieldname = ['Go_licked' int2str(C(k))];
        data.SessionPerformance.(lickfieldname) = [];
        noLickfieldname = ['Go_noLick' int2str(C(k))];
        data.SessionPerformance.(noLickfieldname) = [];
        noLick=0;
        Lick=0;
        
        %Create arrays for GO with Stimulation
        lickfieldnameStim = ['Go_licked_Stim' int2str(C(k))];
        data.SessionPerformance.(lickfieldnameStim) = [];
        noLickfieldnameStim = ['Go_noLick_Stim' int2str(C(k))];
        data.SessionPerformance.(noLickfieldnameStim) = [];
        noLickStim=0;
        LickStim=0;
        
        %Create arrays for GO without Stimulation
        lickfieldnameNoStim = ['Go_licked_NoStim' int2str(C(k))];
        data.SessionPerformance.(lickfieldnameNoStim) = [];
        noLickfieldnameNoStim = ['Go_noLick_NoStim' int2str(C(k))];
        data.SessionPerformance.(noLickfieldnameNoStim) = [];
        noLickNoStim=0;
        LickNoStim=0;
        
        data.SessionPerformance.correctTrial{k,1} = lickfieldname;
        data.SessionPerformance.incorrectTrial{k,1} = noLickfieldname;
        
        %Create arrays for GO with Stimulation
        data.SessionPerformance.correctTrialStim{k,1} = lickfieldnameStim;
        data.SessionPerformance.incorrectTrialStim{k,1} = noLickfieldnameStim;
        
        %Create arrays for GO without Stimulation
        data.SessionPerformance.correctTrialNoStim{k,1} = lickfieldnameNoStim;
        data.SessionPerformance.incorrectTrialNoStim{k,1} = noLickfieldnameNoStim;
        
        
        for j=1:length(data.trial); % Go Stimulus
            if data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==0;  %no lick
                noLick = noLick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    noLickStim = noLickStim+1;
                else % Not Stimulated
                    noLickNoStim = noLickNoStim+1;
                end
            elseif data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==1;
                data.SessionPerformance.licklatencyAngle(Lick+1,k) = C(k);
                data.SessionPerformance.licklatency(Lick+1,k) = data.trial{j,1}.Performance.FirstLickLatency;
                Lick = Lick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    LickStim = LickStim+1;
                    data.SessionPerformance.licklatencyStim(LickStim+1,k) = data.trial{j,1}.Performance.FirstLickLatency;
                else % Not Stimulated
                    LickNoStim = LickNoStim+1;
                    data.SessionPerformance.licklatencyNoStim(LickNoStim+1,k) = data.trial{j,1}.Performance.FirstLickLatency;
                end
            end
        end
    elseif C(k)>midPoint
        
        lickfieldname = ['NoGo_licked' int2str(C(k))];
        data.SessionPerformance.(lickfieldname) = [];
        noLickfieldname = ['NoGo_noLick' int2str(C(k))];
        data.SessionPerformance.(noLickfieldname) = [];
        noLick=0;
        Lick=0;
        
        %Create arrays for NOGO with Stimulation
        lickfieldnameStim = ['NoGo_licked_Stim' int2str(C(k))];
        data.SessionPerformance.(lickfieldnameStim) = [];
        noLickfieldnameStim = ['NoGo_noLick_Stim' int2str(C(k))];
        data.SessionPerformance.(noLickfieldnameStim) = [];
        noLickStim=0;
        LickStim=0;
        
        %Create arrays for NOGO without Stimulation
        lickfieldnameNoStim = ['NoGo_licked_NoStim' int2str(C(k))];
        data.SessionPerformance.(lickfieldnameNoStim) = [];
        noLickfieldnameNoStim = ['NoGo_noLick_NoStim' int2str(C(k))];
        data.SessionPerformance.(noLickfieldnameNoStim) = [];
        noLickNoStim=0;
        LickNoStim=0;
        
        data.SessionPerformance.correctTrial{k,1} = noLickfieldname;
        data.SessionPerformance.incorrectTrial{k,1} = lickfieldname;
        
        %Create arrays for NoGO with Stimulation
        data.SessionPerformance.correctTrialStim{k,1} = noLickfieldnameStim;
        data.SessionPerformance.incorrectTrialStim{k,1} = lickfieldnameStim;
        
        %Create arrays for NoGO without Stimulation
        data.SessionPerformance.correctTrialNoStim{k,1} = noLickfieldnameNoStim;
        data.SessionPerformance.incorrectTrialNoStim{k,1} = lickfieldnameNoStim;
        
        
        for j=1:length(data.trial);
            if data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==0;
                noLick = noLick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    noLickStim = noLickStim+1;
                else % Not Stimulated
                    noLickNoStim = noLickNoStim+1;
                end
            elseif data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==1;
                Lick = Lick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    LickStim = LickStim+1;
                else % Not Stimulated
                    LickNoStim = LickNoStim+1;
                end
            end
        end
        
    elseif C(k)==midPoint
        
        lickfieldname = ['Mid_licked' int2str(C(k))];
        data.SessionPerformance.(lickfieldname) = [];
        noLickfieldname = ['Mid_noLick' int2str(C(k))];
        data.SessionPerformance.(noLickfieldname) = [];
        noLick=0;
        Lick=0;
        
        data.SessionPerformance.correctTrial{k,1} = lickfieldname;
        data.SessionPerformance.incorrectTrial{k,1} = noLickfieldname;
        
        for j=1:length(data.trial);
            if data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==0;
                noLick = noLick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    noLickStim = noLickStim+1;
                else % Not Stimulated
                    noLickNoStim = noLickNoStim+1;
                end
            elseif data.trial{j,1}.Performance.StimType==C(k) && data.trial{j,1}.Performance.OutcomeLick==1;
                Lick = Lick+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulated
                    LickStim = LickStim+1;
                else % Not Stimulated
                    LickNoStim = LickNoStim+1;
                end
            end
        end
    end
    
    data.SessionPerformance.(lickfieldname) = Lick;
    data.SessionPerformance.(noLickfieldname) = noLick;
    data.SessionPerformance.(lickfieldnameStim) = LickStim;
    data.SessionPerformance.(noLickfieldnameStim) = noLickStim;
    data.SessionPerformance.(lickfieldnameNoStim) = LickNoStim;
    data.SessionPerformance.(noLickfieldnameNoStim) = noLickNoStim;
end

%calculate lick latency parameters for each angle on GO/Licked stim

% if length(C)>1
%
% data.SessionPerformance.Meanlicklatency(1,:) = nanmean(data.SessionPerformance.licklatency);
% data.SessionPerformance.Stdlicklatency(1,:) = nanstd(data.SessionPerformance.licklatency);
% data.SessionPerformance.Meanlicklatency(2,:) = nanmean(data.SessionPerformance.licklatencyNoStim);
% data.SessionPerformance.Stdlicklatency(2,:) = nanstd(data.SessionPerformance.licklatencyNoStim);
% data.SessionPerformance.Meanlicklatency(3,:) = nanmean(data.SessionPerformance.licklatencyStim);
% data.SessionPerformance.Stdlicklatency(3,:) = nanstd(data.SessionPerformance.licklatencyStim);
% end


%% Test Pad data.trial and align to reward onset TEST
PreStim=2000; %time stamps are in 10ms bins, so 100 = 1 second 200 = 2seconds
PostStim=1000;
fullSampleTime = PreStim+PostStim;
TimePointsMin = -PreStim:1:PostStim-1;
%varToAlign = {'TargetMinIdx'};
varToAlign = {'AnswerPeriodEndIdx'};
%varToAlign = {'AnswerPeriodStartIdx' 'AnswerPeriodEndIdx'};

data.AlignedDATA.TrialType.HIT = [];
data.AlignedDATA.TrialType.MISS = [];
data.AlignedDATA.TrialType.CorrectRejection = [];
data.AlignedDATA.TrialType.FalseAlarm = [];
data.AlignedDATA.TrialType.HITStim = [];
data.AlignedDATA.TrialType.MISSStim = [];
data.AlignedDATA.TrialType.CorrectRejectionStim = [];
data.AlignedDATA.TrialType.FalseAlarmStim = [];
data.AlignedDATA.TrialType.HITNoStim = [];
data.AlignedDATA.TrialType.MISSNoStim = [];
data.AlignedDATA.TrialType.CorrectRejectionNoStim = [];
data.AlignedDATA.TrialType.FalseAlarmNoStim = [];
data.AlignedDATAallDATAHoloStim = [];


for f = 1:length(varToAlign)
    
    if ~isempty(data.SessionPerformance.correctTrial)
        for k = 1:length(data.SessionPerformance.correctTrial) % create empty field for each stimulus orientation
            orientationfieldname{k,1} = ['AlignedDATAallOrientations' int2str(C(k))];
            data.(orientationfieldname{k,1}).(varToAlign{f}) = [];
        end
    end
    clear data.AlignedDATAallDATA.(varToAlign{f});
    %Counters for Performance outcome
    h1=1;
    m1=1;
    cr1=1;
    fa1=1;
    %Counters for Glabal Stimulus type
    d1=1;
    s1=1;
    mp1=1;
    mnl1=1;
    ml1=1;
    
    %Counters for Performance outcome with and without Stimulation
    hStim1 = 1;
    hNoStim1 = 1;
    mStim1 = 1;
    mNoStim1 = 1;
    crStim1 = 1;
    crNoStim1 = 1;
    faStim1 = 1;
    faNoStim1 = 1;
    
    %Counters for Trials with and without Stimulation
    holo1 = 1;
    noholo1 = 1;
    
    for j=1:length(data.trial);
        raster_matrix = [];
        velocityTemp = diff(data.trial{j,1}.filteredDATA(:,1));
        addPreStim = [];
        addPostStim = [];
        
        
        if (data.trial{j,1}.Stimulus.(varToAlign{f})-PreStim<1)
            addPreStim=NaN*ones(1,PreStim - (data.trial{j,1}.Stimulus.(varToAlign{f})-1)); %%add nans to start of vector
            currPreStim =  PreStim-(PreStim - data.trial{j,1}.Stimulus.(varToAlign{f})+1);
        else
            currPreStim = PreStim;
        end
        
        if length(velocityTemp)<data.trial{j,1}.Stimulus.(varToAlign{f})+PostStim-1
            addPostStim=NaN*ones(1,data.trial{j,1}.Stimulus.(varToAlign{f})+PostStim-1-length(velocityTemp)); %%add nans to start of vector
            currPostStim =  PostStim-((data.trial{j,1}.Stimulus.(varToAlign{f})+PostStim-1)-length(velocityTemp));
        else
            currPostStim = PostStim;
        end
        
        %Get timestamps from End of sample time (target time)
        %ALL data
        
        tempVelocityStim = velocityTemp(data.trial{j,1}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j,1}.Stimulus.(varToAlign{f})+currPostStim-1,1)';
        tempLickStim = (data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j,1}.Stimulus.(varToAlign{f})+currPostStim-1,4))';
        IdxTarget = (data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j,1}.Stimulus.(varToAlign{f})+currPostStim-1,2))';
        IdxReward = (data.trial{j,1}.filteredDATA(data.trial{j,1}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j,1}.Stimulus.(varToAlign{f})+currPostStim-1,3))';
        
        data.AlignedDATAallDATA.(varToAlign{f})(1:fullSampleTime,j,1) = [addPreStim tempVelocityStim addPostStim];
        data.AlignedDATAallDATA.(varToAlign{f})(1:fullSampleTime,j,2) = [addPreStim tempLickStim addPostStim];
        data.AlignedDATAallDATA.(varToAlign{f})(1,j,3) = data.trial{j,1}.Performance.OutcomeLick;
        data.AlignedDATAallDATA.(varToAlign{f})(1,j,4) = data.trial{j,1}.Performance.StimType;
        data.AlignedDATAallDATA.(varToAlign{f})(1:fullSampleTime,j,5) = [addPreStim IdxTarget addPostStim];
        data.AlignedDATAallDATA.(varToAlign{f})(1:fullSampleTime,j,6) = [addPreStim IdxReward addPostStim];
        data.AlignedDATAallDATA.(varToAlign{f})(1,j,7) = data.trial{j,1}.Performance.OutcomeCode;
    end
    
    for j=1:length(data.trial);
        
        %Separate Data output into Go Stim and Distractor and Hit Miss Correct Rejection False Alarm
        %Trials
        if (data.trial{j,1}.Performance.StimTrial ==1); %Stimulated Session
            data.AlignedDATAallDATAHoloStim.(varToAlign{f})(:,holo1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            holo1 = holo1+1;
        else
            data.AlignedDATAallDATAHoloNoStim.(varToAlign{f})(:,noholo1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            noholo1 = noholo1+1;
        end
        
        if (data.trial{j,1}.Performance.StimType < midPoint); %Go Stim
            data.AlignedDATAallDATAStim.(varToAlign{f})(:,s1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            s1=s1+1;
            if data.trial{j,1}.Performance.OutcomeLick==1; %HIT trial
                data.AlignedDATA.TrialType.HIT.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                h1 = h1+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.HITStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    hStim1 = hStim1+1;
                else
                    data.AlignedDATA.TrialType.HITNoStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    hNoStim1 = hNoStim1+1;
                end
            elseif data.trial{j,1}.Performance.OutcomeLick==0 %Miss trial
                data.AlignedDATA.TrialType.MISS.(varToAlign{f})(:,m1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                m1 = m1+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.MISSStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    mStim1 = mStim1+1;
                else
                    data.AlignedDATA.TrialType.MISSNoStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    mNoStim1 = mNoStim1+1;
                end
            end
        elseif (data.trial{j,1}.Performance.StimType > midPoint); %NoGo Stim
            data.AlignedDATAallDATADist.(varToAlign{f})(:,d1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            d1=d1+1;
            if data.trial{j,1}.Performance.OutcomeLick==0; %Correct Rejection Trial
                data.AlignedDATA.TrialType.CorrectRejection.(varToAlign{f})(:,cr1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                cr1 = cr1+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.CorrectRejectionStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    crStim1 = crStim1+1;
                else
                    data.AlignedDATA.TrialType.CorrectRejectionNoStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    crNoStim1 = crNoStim1+1;
                end
                
            elseif data.trial{j,1}.Performance.OutcomeLick==1; %false alarm trial
                data.AlignedDATA.TrialType.FalseAlarm.(varToAlign{f})(:,fa1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                fa1 = fa1+1;
                if data.trial{j,1}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.FalseAlarmStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    faStim1 = faStim1+1;
                else
                    data.AlignedDATA.TrialType.FalseAlarmNoStim.(varToAlign{f})(:,h1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                    faNoStim1 = faNoStim1+1;
                end
            end
            
        elseif (data.trial{j,1}.Performance.StimType == midPoint); %midpoint
            data.AlignedDATAallDATAmidPoint.(varToAlign{f})(:,mp1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            mp1=mp1+1;
            if data.trial{j,1}.Performance.OutcomeLick==0; %didn't lick
                data.AlignedDATA.TrialType.midPointNoLickCount.(varToAlign{f})(:,mnl1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                mnl1 = mnl1+1;
            elseif data.trial{j,1}.Performance.OutcomeLick==1; %false alarm trial
                data.AlignedDATA.TrialType.midPointLickCount.(varToAlign{f})(:,ml1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
                ml1 = ml1+1;
            end
        end
        
        
        for k = 1:length(C)
            if data.trial{j,1}.Performance.StimType==C(k); % distractor 1
                data.(orientationfieldname{k}).(varToAlign{f})(:,size((data.(orientationfieldname{k}).(varToAlign{f})),2)+1,:) = data.AlignedDATAallDATA.(varToAlign{f})(:,j,:);
            end
        end
    end
    
    data = orderfields(data);
    
    %Calculate Block Performance
    blockBinSize = 12;
    blockPerformance = data.AlignedDATAallDATA.(varToAlign{f})(1,:,7);
    blockPerformance(blockPerformance==4)=1;
    blockPerformance(blockPerformance==1)=1;
    blockPerformance(blockPerformance==2)=2;
    blockPerformance(blockPerformance==3)=2;
    binspacing = 1:blockBinSize:size((blockPerformance),2);
    
    for s = 1:length(binspacing)-1;
        data.SessionPerformance.PerformanceBlock(s) = size((find(blockPerformance(1,binspacing(s):binspacing(s+1)-1)==1)),2)/blockBinSize;
    end
    
end

%%

fontSize = 16;

%varToPlot = {'AlignedDATAallDATA'};
%varToPlot = {'AlignedDATAallDATAStim'};
%varToPlot = {'AlignedDATAallDATAStim' 'AlignedDATAallDATADist'};
%varToPlot = {'AlignedDATAallDATA' 'AlignedDATAallDATAStim' 'AlignedDATAallDATADist'}

%varToPlot = {'AlignedDATAallDATAHoloStim'};
varToPlot = {'AlignedDATAallDATA' 'AlignedDATAallDATAHoloStim' 'AlignedDATAallDATAHoloNoStim'};

f=1;
TimeScaled= TimePointsMin;
currFig = nan;

for vp = 1:length(varToPlot);
    currFig(vp) = figure(vp+1);
    set(currFig(vp),'name',varToPlot{vp},'numbertitle','off');
    set(currFig(vp),'Position',positionGraph1)
    clf;
    %raster of licking
    clear raster_matrix;
    spike_matrix = data.(varToPlot{vp}).(varToAlign{f})(:,:,2)';
    
    for i=1:size((spike_matrix),1)
        raster_matrix(i,:) = spike_matrix(i,:).*i;
    end
    
    subplot(5,3,[1 2 4 5 7 8])
    %Velocity
    imagesc(TimeScaled,1,data.(varToPlot{vp}).(varToAlign{f})(:,:,1)')
    hold on
    plot(TimeScaled,raster_matrix,'k.','LineWidth',4);
    startLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    endLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    lickLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    
    for k = 1:size(data.(varToPlot{vp}).(varToAlign{f}),2);
        startLine(k,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,5)==1));
        endLine(k,1) = max(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,5)==1));
        if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,6)==1,1))
            lickLine(k,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,6)==1));
        else
            lickLine(k,1) = nan;
        end
    end
    
    for k = 1:size(data.(varToPlot{vp}).(varToAlign{f}),2);
        line([TimeScaled(startLine(k)) TimeScaled(startLine(k))], [k-1 k],'Color','g', 'LineWidth',2);
        line([TimeScaled(endLine(k)) TimeScaled(endLine(k))], [k-1 k],'Color','g', 'LineWidth',2);
        if (isnan(lickLine(k,1))==0) %if it doesn't equal nan
            line([TimeScaled(lickLine(k,1)) TimeScaled(lickLine(k,1))], [k-1 k],'Color','r', 'LineWidth',3);
        end
    end
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16);
    xlim([min(TimeScaled), max(TimeScaled)]);
    ylim([1 size(data.(varToPlot{vp}).(varToAlign{f}),2)])
    ylabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    
    subplot(5,3,[10 11]) %Plot Average Velocity
    AvgVelocity = nanmean(data.(varToPlot{vp}).(varToAlign{f})(:,:,1),2);
    enmmptyVel(vp,1:length(AvgVelocity)) = AvgVelocity;
    plot(TimeScaled,AvgVelocity,'k','LineWidth',2);
    set(gca,'FontSize',16);
    xlim([min(TimeScaled), max(TimeScaled)])
    ylabel('Velocity', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    %ylim([0 20]);
    
    subplot(5,3,[13 14]) %Average licking
    binSize = 200;
    t = TimeScaled;
    tpoints =linspace(min(t),max(t),binSize);
    y = nansum(data.(varToPlot{vp}).(varToAlign{f})(:,:,2),2);
    [~,bin]=histc(t,linspace(min(t),max(t),binSize));
    plot(tpoints,full(mean(sparse(1:length(t),bin,y))),'k','LineWidth',2)
    %     plot(TimeScaled,sumLicking,'k','LineWidth',2);
    set(gca,'FontSize',16);
    xlim([min(TimeScaled), max(TimeScaled)])
    ylabel('psth Licking', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    
    subplot(5,3,[3 6 9])
    for i=1:size(data.(varToPlot{vp}).(varToAlign{f}),2)
        if data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)<=midPoint; % Correct Lick on GO stimulus = Hit
            line([0 1], [i i], 'Color', 'g','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)<=midPoint; % No Lick on GO stimulus = Miss
            line([0 1], [i i], 'Color', 'r','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)>=midPoint; % No Lick on Distractor stimulus = Correct Rejection
            line([2 3], [i i], 'Color', 'g','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)>=midPoint; % Lick on Distractor stimulus = False Alarm
            line([2 3], [i i], 'Color', 'r','LineWidth',2);
            hold on;
        end
    end
    
    xlim([-0.5 3.5])
    ylim([1 size(data.(varToPlot{vp}).(varToAlign{f}),2)])
    set(gca,'YDir','reverse');
    
    if length(varToPlot)>1
        figure(99)
        subplot(2,1,1) %Plot Average Velocity
        plot(TimeScaled,AvgVelocity,'LineWidth',2);
        hold on
        set(gca,'FontSize',16);
        xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('Velocity', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        legend(varToPlot)
        %ylim([0 20]);
        
        subplot(2,1,2) %Average licking
        plot(tpoints,full(mean(sparse(1:length(t),bin,y))),'LineWidth',2)
        hold on
        set(gca,'FontSize',16);
        xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('psth Licking', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        legend(varToPlot)
    end
    
    
end

figure(21)
set(gcf,'Position',positionGraph2)
subplot(2,1,1);
plot(data.SessionPerformance.CumulativePerformance,'*b');
ylabel('cumulative performance', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
ylim([0,1])

subplot(2,1,2);
plot(data.SessionPerformance.PerformanceBlock,'-b','LineWidth',2);
ylabel('block performance', 'FontSize', fontSize,'fontWeight','bold');
xlabel('block (12) number', 'FontSize', fontSize,'fontWeight','bold');
ylim([0,1]);
%Run group analysis
%
% [sessPerformance d1] = JB_sessionPerformanceSinglePlot(data)

end




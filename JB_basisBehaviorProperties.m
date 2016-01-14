
load('DATA.mat')
basisPropertiesID{1,1} = 'sessionType';
sessionType = 1;
basisPropertiesID{2,1} = 'date';
date = 2;
basisPropertiesID{3,1} = 'datenum';
datenum = 3;
basisPropertiesID{4,1} = 'sessionDuration';
sessionDuration = 4;
basisPropertiesID{5,1} = 'totalStepsPerMin';
totalStepsPerMin = 5;

basicProperties = {};
encoder0Pos = 1;
rawSessionTime = 7;

for i=1:length(DATA.allFiles);
    
    basicProperties{date,i} = DATA.allFiles{i}.date;
    basicProperties{datenum,i} = DATA.allFiles{i}.datenum;
    
    tempDATA = DATA.allFiles{i}.rawData;
    basicProperties{sessionDuration,i} = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
    threshold = 50000;
    aboveThreshold = find(tempDATA(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        tempDATA(aboveThreshold(j),:)=nan;
    end
    
    totalSteps = nansum(diff(tempDATA(:,1))>0);
    basicProperties{totalStepsPerMin,i} = totalSteps/basicProperties{sessionDuration,i}; % Calculate running velocity
    
    %import as metaDATA cell array
    n=1;
    
    metaData = DATA.allFiles{i}.metaData;
    for j=1:length(metaData)
        
        match(j,1) = strcmp('Orientation Selected = ', metaData{j,1});
        
        if (match(j,1)==1) && (metaData{j,2}>0)
            tempAngles(n,1) = metaData{j,2};
            n=n+1;
        end
        
        %AutoReward
        if(strcmp('Auto reward = ', metaData{j,1}));
            autoSession = metaData{j,2};
        end
        
    end
    
    %Define session type
    tempsessionType = strcat('S',num2str(length(tempAngles)));
    
    if autoSession==1
        tempsessionType = strcat(tempsessionType,'auto');
    end
    
    basicProperties{sessionType,i} = tempsessionType;
    
end



for k = 1:length(basicProperties)
    
    days{k} = basicProperties{2,k}(1:11)
    
end

plotNumber = 1;

[un idx_last idx] = unique(days(1,:));
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});


for k = 1:length(unique_idx)
    
    if length(unique_idx{k})==1
        
        basicPropertiesToPlot(:,plotNumber) = basicProperties(:,(unique_idx{k}));
        
        plotNumber = plotNumber+1;
        
    elseif length(unique_idx{k})>1
        
        t = unique_idx{k};
        
        for kk = 1:length(t)
            
            tempRunVal(kk,1) = basicProperties{4,t(kk)};
            
        end
        
        [id di] = max(tempRunVal); %find the data set where the mouse ran the most
        
        basicPropertiesToPlot(:,plotNumber) = basicProperties(:,t(di));
        
        plotNumber = plotNumber+1;
        
    end
    
end


figure(1);clf
subplot(3,1,1)

numPoints = 1:1:length(basicPropertiesToPlot);
for j = 1:length(basicPropertiesToPlot);
    plot(numPoints(j),basicPropertiesToPlot{5,j},'or','MarkerSize', 10,'MarkerFaceColor','r')
    hold on
end

ylabel('totalStepsPerMin');
xlabel('Session Number');

subplot(3,1,2)

numPoints = 1:1:length(basicPropertiesToPlot);
for j = 1:length(basicPropertiesToPlot);
    plot(numPoints(j),basicPropertiesToPlot{4,j},'or','MarkerSize', 10,'MarkerFaceColor','r')
    hold on
end

ylabel('sessionDuration');
xlabel('Session Number');

%plot primary session type per day
subplot(3,1,3)
SessionTypes = {'S1auto' ; 'S1'; 'S2'; 'S6'; 'S12'};

for i=1:length(basicPropertiesToPlot)
    if strcmp('S1auto', basicPropertiesToPlot{1,i})
        colorCode(i,1) = 1; 
    elseif strcmp('S1',basicPropertiesToPlot{1,i})
        colorCode(i,1) = 2;
    elseif strcmp('S2',basicPropertiesToPlot{1,i})
        colorCode(i,1) = 3;
    elseif strcmp('S6',basicPropertiesToPlot{1,i})
        colorCode(i,1) = 4;
    elseif strcmp('S12',basicPropertiesToPlot{1,i})
        colorCode(i,1) = 5;
    end
end

ColorCodeColors = [1.0;0.8;0.6;0.4;0.2;0];

numPoints = 1:1:length(basicPropertiesToPlot);
for j = 1:length(basicPropertiesToPlot);
    plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
    hold on
end

NumTicks = max(colorCode);
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',SessionTypes(1:max(colorCode)))
%set(gca,'YTickLabel',{'S1Auto','S1', 'S2'})
ylabel('Session Type');
xlabel('Session Number');






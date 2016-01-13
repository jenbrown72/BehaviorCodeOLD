function [sessPerformance d1] = JB_sessionPerformanceSinglePlot(data)
%%
%load('data(i).mat')
fontSize = 16; %font Size for plotting

% for i = 1:size((data(i)),2)
%       for i = 3;

i=1;

o=1;
%Adding a 1 to each as to calculate d' each need to be an integer
hit = data(i).SessionPerformance.HitCount+1;
miss = data(i).SessionPerformance.MissCount+1;
FA = data(i).SessionPerformance.FalseAlarmCount+1;
CR = data(i).SessionPerformance.CorrectRejectionCount+1;

hitStim = data(i).SessionPerformance.HitCountStim+1;
missStim = data(i).SessionPerformance.MissCountStim+1;
FAStim = data(i).SessionPerformance.FalseAlarmCountStim+1;
CRStim = data(i).SessionPerformance.CorrectRejectionCountStim+1;

hitNoStim = data(i).SessionPerformance.HitCountNoStim+1;
missNoStim = data(i).SessionPerformance.MissCountNoStim+1;
FANoStim = data(i).SessionPerformance.FalseAlarmCountNoStim+1;
CRNoStim = data(i).SessionPerformance.CorrectRejectionCountNoStim+1;

sessPerformance(i,1) = (hit+CR)/(data(i).SessionPerformance.totalSessions);
sessPerformanceStim(i,1) = (hitStim+CRStim)/(hitStim+CRStim+FAStim+missStim);
sessPerformanceNoStim(i,1) = (hitNoStim+CRNoStim)/(hitNoStim+CRNoStim+FANoStim+missNoStim);

%D' - discrimitability index
d1(i,1)=norminv((hit/(hit+miss)),0,1)-norminv((FA/(FA+CR)),0,1);
d1Stim(i,1)=norminv((hitStim/(hitStim+missStim)),0,1)-norminv((FAStim/(FAStim+CRStim)),0,1);
d1NoStim(i,1)=norminv((hitNoStim/(hitNoStim+missNoStim)),0,1)-norminv((FANoStim/(FANoStim+CRNoStim)),0,1);

%for NoGo stim
tempCorrect = flipud(data(i).SessionPerformance.correctTrial);
tempInCorrect = flipud(data(i).SessionPerformance.incorrectTrial);
tempAngles = flipud(data(i).SessionPerformance.orientations);
TF = strncmpi(tempCorrect,'No',2); %Correctly didn't lick

%For stimulation
tempCorrectStim = flipud(data(i).SessionPerformance.correctTrialStim);
tempInCorrectStim = flipud(data(i).SessionPerformance.incorrectTrialStim);

%For No stimulation
tempCorrectNoStim = flipud(data(i).SessionPerformance.correctTrialNoStim);
tempInCorrectNoStim = flipud(data(i).SessionPerformance.incorrectTrialNoStim);

for k = 1:length(TF); % for each stim calculate percent correct
    
    if TF(k)==1;
        
        CRDist = data(i).SessionPerformance.(tempCorrect{k,1});
        FADist = data(i).SessionPerformance.(tempInCorrect{k,1});
        sessPerformanceALL(i,o) = CRDist/(CRDist+FADist)*100; % calculate percent correct for each orientation
        sessProbLickALL(i,o) = FADist/(CRDist+FADist)*100; % calculate percent correct for each orientation
        
        CRDistStim = data(i).SessionPerformance.(tempCorrectStim{k,1});
        FADistStim = data(i).SessionPerformance.(tempInCorrectStim{k,1});
        sessPerformanceALLStim(i,o) = CRDistStim/(CRDistStim+FADistStim)*100; % calculate percent correct for each orientation
        sessProbLickALLStim(i,o) = FADistStim/(CRDistStim+FADistStim)*100; % calculate percent correct for each orientation
        
        CRDistNoStim = data(i).SessionPerformance.(tempCorrectNoStim{k,1});
        FADistNoStim = data(i).SessionPerformance.(tempInCorrectNoStim{k,1});
        sessPerformanceALLNoStim(i,o) = CRDistNoStim/(CRDistNoStim+FADistNoStim)*100; % calculate percent correct for each orientation
        sessProbLickALLNoStim(i,o) = FADistNoStim/(CRDistNoStim+FADistNoStim)*100; % calculate percent correct for each orientation
        
        orientationsToPlot(o) = tempAngles(k);
        o = o+1;
        
    end
end

%For Mid Point stim
%     tempCorrect = data(i).SessionPerformance.correctTrial;
%     tempInCorrect = data(i).SessionPerformance.incorrectTrial;
%         tempAngles = data(i).SessionPerformance.orientations;

TF = strncmpi(tempCorrect,'Mid',3);

for k = 1:length(TF); % for each stim calculate percent correct
    
    if TF(k)==1;
        
        HIT = data(i).SessionPerformance.(tempCorrect{k,1});
        MISS = data(i).SessionPerformance.(tempInCorrect{k,1});
        sessPerformanceALL(i,o) = HIT/(HIT+MISS)*100; % calculate percent correct for each orientation
        sessProbLickALL(i,o) = HIT/(HIT+MISS)*100; % calculate percent correct for each orientation
        
        HITStim = data(i).SessionPerformance.(tempCorrectStim{k,1});
        MISSStim = data(i).SessionPerformance.(tempInCorrectStim{k,1});
        sessPerformanceALLStim(i,o) = HITStim/(HITStim+MISSStim)*100; % calculate percent correct for each orientation
        sessProbLickALLStim(i,o) = HITStim/(HITStim+MISSStim)*100; % calculate percent correct for each orientation
        
        HITNoStim = data(i).SessionPerformance.(tempCorrectNoStim{k,1});
        MISSNoStim = data(i).SessionPerformance.(tempInCorrectNoStim{k,1});
        sessPerformanceALLNoStim(i,o) = HITNoStim/(HITNoStim+MISSNoStim)*100; % calculate percent correct for each orientation
        sessProbLickALLNoStim(i,o) = HITNoStim/(HITNoStim+MISSNoStim)*100; % calculate percent correct for each orientation
        
        orientationsToPlot(o) = tempAngles(k);
        o = o+1;
        
    end
end


%For Go stim
%     tempCorrect = flipud(data(i).SessionPerformance.correctTrial);
%     tempInCorrect = flipud(data(i).SessionPerformance.incorrectTrial);
%         tempAngles = flipud(data(i).SessionPerformance.orientations);

TF = strncmpi(tempCorrect,'Go',2);

for k = 1:length(TF); % for each stim calculate percent correct
    
    if TF(k)==1;
        
        HIT = data(i).SessionPerformance.(tempCorrect{k,1});
        MISS = data(i).SessionPerformance.(tempInCorrect{k,1});
        sessPerformanceALL(i,o) = HIT/(HIT+MISS)*100; % calculate percent correct for each orientation
        sessProbLickALL(i,o) = HIT/(HIT+MISS)*100; % calculate percent correct for each orientation
        
        HITStim = data(i).SessionPerformance.(tempCorrectStim{k,1});
        MISSStim = data(i).SessionPerformance.(tempInCorrectStim{k,1});
        sessPerformanceALLStim(i,o) = HITStim/(HITStim+MISSStim)*100; % calculate percent correct for each orientation
        sessProbLickALLStim(i,o) = HITStim/(HITStim+MISSStim)*100; % calculate percent correct for each orientation
        
        HITNoStim = data(i).SessionPerformance.(tempCorrectNoStim{k,1});
        MISSNoStim = data(i).SessionPerformance.(tempInCorrectNoStim{k,1});
        sessPerformanceALLNoStim(i,o) = HITNoStim/(HITNoStim+MISSNoStim)*100; % calculate percent correct for each orientation
        sessProbLickALLNoStim(i,o) = HITNoStim/(HITNoStim+MISSNoStim)*100; % calculate percent correct for each orientation
        
        orientationsToPlot(o) = tempAngles(k);
        
        o = o+1;
        
    end
end

%%find discrimitability index for each pair of angles
numofStim = size((tempAngles),1);

for j=1:(size((tempAngles),1)/2);
    
    %check that you are comparing the correct angles
    
    
    if (tempAngles(max(size(tempAngles)-(j-1)))-tempAngles(j))<2
        
        anglesDiff(1,j) = tempAngles(max(size(tempAngles)-(j-1)))-270;
        
        tempAngleHit = data(i).SessionPerformance.(tempCorrect{max(numofStim)-(j-1)})+1;
        tempAngleMiss = data(i).SessionPerformance.(tempInCorrect{max(numofStim)-(j-1)})+1;
        tempAngleFA = data(i).SessionPerformance.(tempInCorrect{(j)})+1;
        tempAngleCR = data(i).SessionPerformance.(tempCorrect{(j)})+1;
        
        d1PairedAngle(i,j)=norminv((tempAngleHit/(tempAngleHit+tempAngleMiss)),0,1)-norminv((tempAngleFA/(tempAngleFA+tempAngleCR)),0,1);
        
    else
        
        display('error in angle comparison');
        
    end
    
    
end

% end

%convert orientations into difference from middle (3200)
for i=1:length(orientationsToPlot)
    orientationsToPlot(i) = orientationsToPlot(i)-270;
end


[orientationsToPlotTwo,Id] = sort(orientationsToPlot);

%Sort sessPerformanceDist using Id

stepAnlges = 3200;
anglePerStep = 3200/360;
%anglesToPlot = (orientationsToPlot*-1)/anglePerStep;



figure(11)
clf
subplot(3,2,1) %plot of % performance
plot(sessPerformance,'o-b','LineWidth',2);
set(gca,'FontSize',16);
title('Performance over sessions', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
ylabel('Percent Correct ((hit+CR)/totalTrials', 'FontSize', fontSize,'fontWeight','bold');
ylim([0 1]);
%xlimit = xlim;
line([1 max(xlim)],[0.8 0.8],'Color','r','LineStyle','--','LineWidth',2)
colorstring = 'bmgyrgr';

hold on
subplot(3,2,3) %plot of d'

for i = 1:size((d1),2)
    plot(d1(:, i),'o-','LineWidth',2, 'Color', colorstring(i));
    hold on
end

set(gca,'FontSize',16);
title('Discrimitability Index', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
% ylimit = ylim;
% ylim([0 ylimit(2)]);
line([1 max(xlim)],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);


subplot(3,2,2)
for i = 1:size((sessPerformanceALL),1);
    tempSess = sessPerformanceALL(i,Id);
    %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
    plot(orientationsToPlotTwo,tempSess,'o-','LineWidth',2,'Color', [0 0 (1/i)]);
    
    hold on
    
    
end
%
set(gca,'FontSize',16);
%set(gca, 'Xdir', 'reverse')
title('Performance for all orientation stimuli', 'FontSize', fontSize,'fontWeight','bold');
xlabel('G0 - Orientation Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
ylabel('Performance', 'FontSize', fontSize,'fontWeight','bold');
ylim([0 100])
xlimit = xlim;
line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);


subplot(3,2,4)

for i = 1:size((sessProbLickALL),1);
    tempSess = sessProbLickALL(i,Id);
    %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
    plot(orientationsToPlotTwo,tempSess,'o-','LineWidth',2,'Color', [0 0 (1/i)]);
    
    hold on
    
    
end
%
set(gca,'FontSize',16);
%set(gca, 'Xdir', 'reverse')
title('Performance for all orientation stimuli', 'FontSize', fontSize,'fontWeight','bold');
xlabel('G0 - Orientation Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
ylim([0 100])
xlimit = xlim;
line([xlimit(1) xlimit(2)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);



%Plot Stim V No Stim
subplot(3,2,5)

%for i=1:size(size(d1PairedAngle),1);

plot(anglesDiff',d1PairedAngle(:,:)','bo-','LineWidth',2);

set(gca,'FontSize',16);
title('Discrimitability Index/Angle', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
% ylimit = ylim;
% ylim([0 ylimit(2)]);
xlim([(min(anglesDiff)--(max(anglesDiff))) 0])
line([(min(anglesDiff)--(max(anglesDiff))) 0],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);





if data(i).SessionPerformance.optogenetics==1
    
    figure(14)
    clf
    subplot(2,2,1) %plot of % performance
    plot(sessPerformanceStim,'o-r','LineWidth',2);
    hold on
    plot(sessPerformanceNoStim,'o-k','LineWidth',2);
    
    set(gca,'FontSize',16);
    title('Performance over sessions', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('Percent Correct ((hit+CR)/totalTrials', 'FontSize', fontSize,'fontWeight','bold');
    ylim([0 1]);
    line([1 max(xlim)],[0.8 0.8],'Color','r','LineStyle','--','LineWidth',2)
    
    
    
    colorstring = 'bmgyrgr';
    
    hold on
    subplot(2,2,3) %plot of d'
    
    for i = 1:size((d1),2)
        plot(d1Stim(:, i),'o-r','LineWidth',2);
        hold on
        plot(d1NoStim(:, i),'o-k','LineWidth',2);
        
    end
    
    set(gca,'FontSize',16);
    title('Discrimicability Index', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
    % ylimit = ylim;
    % ylim([0 ylimit(2)]);
    line([1 max(xlim)],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);
    
    
    subplot(2,2,2)
    for i = 1:size((sessPerformanceALLNoStim),1);
        tempSessNoStim = sessPerformanceALLNoStim(i,Id);
        %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
        plot(orientationsToPlotTwo,tempSessNoStim,'o-k','LineWidth',2);
        hold on
        
        tempSessStim = sessPerformanceALLStim(i,Id);
        %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
        plot(orientationsToPlotTwo,tempSessStim,'o-r','LineWidth',2);
        hold on
        
        
    end
    %
    set(gca,'FontSize',16);
    %set(gca, 'Xdir', 'reverse')
    title('Performance for all orientation stimuli', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('G0 - Orientation Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('Performance', 'FontSize', fontSize,'fontWeight','bold');
    ylim([0 100])
    line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    
    
    subplot(2,2,4)
    
    for i = 1:size((sessProbLickALLNoStim),1);
        tempSessNoStim = sessProbLickALLNoStim(i,Id);
        %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
        plot(orientationsToPlotTwo,tempSessNoStim,'o-k','LineWidth',2);
        
        hold on
        
        tempSessStim = sessProbLickALLStim(i,Id);
        %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
        plot(orientationsToPlotTwo,tempSessStim,'o-r','LineWidth',2);
        
        
        hold on
    end
    %
    set(gca,'FontSize',16);
    %set(gca, 'Xdir', 'reverse')
    title('Performance for all orientation stimuli', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('G0 - Orientation Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
    ylim([0 100])
    line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    
end

%For IGOR

nIGOR_data = [orientationsToPlotTwo' sessPerformanceALL' sessProbLickALL'];
nIGOR_data2 = [sessPerformance d1];
nIGOR_data = [orientationsToPlotTwo' sessPerformanceALLNoStim' sessPerformanceALLStim' sessProbLickALLNoStim' sessProbLickALLStim'];

nIGOR_sessPerf = [sessPerformanceNoStim';sessPerformanceStim';d1NoStim';d1Stim'];
nIGOR_sessDis = [d1NoStim';d1Stim'];

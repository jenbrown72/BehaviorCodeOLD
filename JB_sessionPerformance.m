function [sessPerformance d1] = JB_sessionPerformance(filename,fileSession)
    %%
    load('DATA.mat')
    clear data
    k=1;
    Tempdata = DATA.(filename);
    fontSize = 12; %font Size for plotting
    sessionsToAnalyse = [];

    % if nargin ==2;
    %     
    % sessionsToAnalyse = fileSession;
    % 
    % end

    for kk = 1:size((Tempdata),2);
        stored=0;

        if ~isempty(sessionsToAnalyse)
            for jj = 1:size((sessionsToAnalyse),2)

                if ((kk==sessionsToAnalyse(jj)) && (stored==0))

        data(1,k) = Tempdata{kk}.analysedDATA;
        stored = 1;
        k=k+1;
                end

            end

        else
            data(k) = Tempdata{kk}.analysedDATA;
        k=k+1;

        end

    end


    for i = 1:size((data),2)
        %       for i = 3;

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

    end

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


    subplot(3,2,5)

    plot(anglesDiff',d1PairedAngle(:,:)','bo-');

    set(gca,'FontSize',16);
    title('Discrimitability Index/Angle', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Angle Difference from 0', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
    % ylimit = ylim;
    xlim([(min(anglesDiff)--(max(anglesDiff))) 0])
    line([(min(anglesDiff)--(max(anglesDiff))) 0],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);



    subplot(3,2,6)
    avgAngleD1 = mean(d1PairedAngle);

    plot(anglesDiff',avgAngleD1,'bo-','LineWidth',2);


    set(gca,'FontSize',16);
    title('Avg Discrimitability Index/Angle', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Angle Difference from 0', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
    % ylimit = ylim;
    xlim([(min(anglesDiff)--(max(anglesDiff))) 0])
    line([(min(anglesDiff)--(max(anglesDiff))) 0],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);



    %Plot Stim V No Stim

    %check to see if these involved stim and non stim trials


    x = length(sessPerformanceStim);

    figure(14)
    clf
    subplot(2,x+1,1) %plot of % performance
    plot(sessPerformanceStim,'o-r','LineWidth',4);
    hold on
    plot(sessPerformanceNoStim,'o-k','LineWidth',4);

    set(gca,'FontSize',fontSize);
    %title('Performance', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('Percent Correct', 'FontSize', fontSize,'fontWeight','bold');
    ylim([0 1]);
    line([1 max(xlim)],[0.8 0.8],'Color','r','LineStyle','--','LineWidth',2)
    line([1 max(xlim)],[0.5 0.5],'Color',[0.4 0.4 0.4]','LineStyle','--','LineWidth',2)

    linethickness = 2;

    hold on
    subplot(2,x+1,x+2) %plot of d'

    for i = 1:size((d1),2)
        plot(d1Stim(:, i),'o-r','LineWidth',4);
        hold on
        plot(d1NoStim(:, i),'o-k','LineWidth',4);

    end

    set(gca,'FontSize',fontSize);
    %title('Discriminability Index', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Session number', 'FontSize', fontSize,'fontWeight','bold');
    ylabel('d''', 'FontSize', fontSize,'fontWeight','bold');
    % ylimit = ylim;
    % ylim([0 ylimit(2)]);
    line([1 max(xlim)],[1.75 1.75],'Color','r','LineStyle','--','LineWidth',2);
    ColorShade = 0:1/size((sessPerformanceALLNoStim),1):1;

    %ColorShade = fliplr(0:1/size((sessPerformanceALLNoStim),1):1);

    for kk = 1:x; %for the lenght of sessions to plot)
        
        perfromanceHandles = subplot(2,x+1,kk+1);
        tempSessNoStim = sessPerformanceALLNoStim(kk,Id);
        tempSessStim = sessPerformanceALLStim(kk,Id);
        
        plot(orientationsToPlotTwo,tempSessNoStim,'k','LineWidth',linethickness);
        hold on
        plot(orientationsToPlotTwo,tempSessStim,'r','LineWidth',linethickness);
        
        xlabel(perfromanceHandles,'G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
        ylabel(perfromanceHandles,'Performance', 'FontSize', fontSize,'fontWeight','bold');
        line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
        line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
        
        lickingHandles = subplot(2,x+1,kk+(x+2));
        tempSessNoStim = sessProbLickALLNoStim(kk,Id);
        tempSessStim = sessProbLickALLStim(kk,Id);
        
        plot(orientationsToPlotTwo,tempSessNoStim,'k','LineWidth',linethickness);
        hold on
        plot(orientationsToPlotTwo,tempSessStim,'r','LineWidth',linethickness);
        
        xlabel(lickingHandles,'G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
           ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
        line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
        line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
        
        
    end

%                 %    plot(orientationsToPlotTwo,sessPerformanceALLNoStim(1,Id),'Color',[0 0 1],'LineWidth',4);
% 
%     %
%     set(perfromanceHandles(1),'FontSize',16);
%     %set(gca, 'Xdir', 'reverse')
%     title(perfromanceHandles(1),'Non Stimulated Trials', 'FontSize', fontSize,'fontWeight','bold');
%     xlabel(perfromanceHandles(1),'G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
%     ylabel(perfromanceHandles(1),'Performance', 'FontSize', fontSize,'fontWeight','bold');
%     ylim(perfromanceHandles(1),[0 100])
%     line(perfromanceHandles(1),[-50 50],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
%     line(perfromanceHandles(1),[0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% 
%     
%     set(lickingHandles(1),'FontSize',16);
%     %set(gca, 'Xdir', 'reverse')
%     title('Non Stimulated Trials', 'FontSize', fontSize,'fontWeight','bold');
%     xlabel('G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
%     ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
%     ylim([0 100])
%     line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
%     line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
%     
    

% subplot(2,x,)
% for i = 1:size((sessPerformanceALLNoStim),1);
% 
%         tempSessStim = sessPerformanceALLStim(i,Id);
%     %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
%         plot(orientationsToPlotTwo,tempSessStim,'Color',[greyShade greyShade greyShade],'LineWidth',linethickness);
% 
%    % plot(orientationsToPlotTwo,tempSessStim,'o-r','LineWidth',2);
%     hold on
% end
% 
%             %    plot(orientationsToPlotTwo,sessPerformanceALLStim(1,Id),'Color',[1 0 0],'LineWidth',4);
% 
% %
% set(gca,'FontSize',16);
% %set(gca, 'Xdir', 'reverse')
% title('Stimulated Trials', 'FontSize', fontSize,'fontWeight','bold');
% xlabel('G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
% ylabel('Performance', 'FontSize', fontSize,'fontWeight','bold');
% ylim([0 100])
% line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% 
% 
% subplot(2,3,5)
% 
% for i = 1:size((sessProbLickALLNoStim),1);
%     tempSessNoStim = sessProbLickALLNoStim(i,Id);
%     %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
%     plot(orientationsToPlotTwo,tempSessNoStim,'Color',[greyShade greyShade greyShade],'LineWidth',linethickness);
%     
%     hold on
%     
% end
% 
%      %           plot(orientationsToPlotTwo,sessProbLickALLNoStim(1,Id),'Color',[0 0 1],'LineWidth',4);
% 
% 
% 
% %
% set(gca,'FontSize',16);
% %set(gca, 'Xdir', 'reverse')
% title('Non Stimulated Trials', 'FontSize', fontSize,'fontWeight','bold');
% xlabel('G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
% ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
% ylim([0 100])
% line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% 
% 
% subplot(2,3,6)
% 
% for i = 1:size((sessProbLickALLNoStim),1);
% 
%     
%     tempSessStim = sessProbLickALLStim(i,Id);
%     %  plot(anglesToPlot,sessPerformanceDist(i,:),'o-','LineWidth',2,'Color', colorstring(i));
%     plot(orientationsToPlotTwo,tempSessStim,'Color',[greyShade greyShade greyShade],'LineWidth',linethickness);
%     
%     
%     hold on
% end
% 
% % plot(orientationsToPlotTwo,sessProbLickALLStim(1,Id),'Color',[1 0 0],'LineWidth',4);
% %
% set(gca,'FontSize',16);
% %set(gca, 'Xdir', 'reverse')
% title('Stimulated Trials', 'FontSize', fontSize,'fontWeight','bold');
% xlabel('G0 - Stim Angle - NoGo', 'FontSize', fontSize,'fontWeight','bold');
% ylabel('Probability of Lick', 'FontSize', fontSize,'fontWeight','bold');
% ylim([0 100])
% line([min(xlim) max(xlim)],[50 50],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% line([0 0],[0 100],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
% 
% 
% 


%For IGOR

nIGOR_data = [orientationsToPlotTwo' sessPerformanceALL' sessProbLickALL'];
nIGOR_data2 = [sessPerformance d1];
nIGOR_data = [orientationsToPlotTwo' sessPerformanceALLNoStim' sessPerformanceALLStim' sessProbLickALLNoStim' sessProbLickALLStim'];

nIGOR_sessPerf = [sessPerformanceNoStim';sessPerformanceStim';d1NoStim';d1Stim'];
nIGOR_sessDis = [d1NoStim';d1Stim'];
fi
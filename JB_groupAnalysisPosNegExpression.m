posMouseIDs = {'JB023_NT';'JB025_CE';'JB023_L';'JB024_NT';'JB023_R'};
negMouseIDs = {'JB025_R';'JB025_RL';'JB026_CE';'JB026_L';'JB026_R';'JB025_L';'JB024_R'};
 %posMouseIDs = {'JB025_CE'};
% negMouseIDs = {'JB025_R'};
variableIDs = {'dprime';'lick bias';'performance'};

for g=1:2; %ipsi or contra
    plotNoPos=1;
    plotNoNeg=1;
    pairs = [1,2;3,4;5,6];
    plotNO=1;
    
    figure(11*g)
    if g==1
        set(gcf,'name','Halo +ve Contralateral Stimulation','numbertitle','off')
    else
        set(gcf,'name','Halo +ve Ipsilateral Stimulation','numbertitle','off')
    end
    
    
    figure(12*g)
    if g==1
        set(gcf,'name','Halo -ve Contralateral Stimulation','numbertitle','off')
    else
        set(gcf,'name','Halo -ve Ipsilateral Stimulation','numbertitle','off')
    end
    
    
    figure(111*g)
    if g==1
        set(gcf,'name','Avg Contralateral Stimulation','numbertitle','off')
    else
        set(gcf,'name','Avg Ipsilateral Stimulation','numbertitle','off')
    end
    
    n=1;
    clear posIdx;
    clear negIdx;
    for i=1:length(posMouseIDs);
        for ii = 1:length(mouseIDs)
            if(strcmp(posMouseIDs{i},mouseIDs{ii,g}))
                posIdx(n,1) = ii;
                n=n+1;
            end
        end
    end
    
    n=1;
    for i=1:length(negMouseIDs);
        for ii = 1:length(mouseIDs)
            if(strcmp(negMouseIDs{i},mouseIDs{ii,g}))
                negIdx(n,1) = ii;
                n=n+1;
            end
        end
    end
    
    
    for j = 1:length(pairs);
        clear avgDATA semDATA stim;
        figure(11*g)
        hold on
        for i = 1:length(posIdx);
            clear stim
            if (~isempty(DATA{posIdx(i),g}))
                stim(1,:) = DATA{posIdx(i),g}.data(:,pairs(j,1))';
                stim(2,:) = DATA{posIdx(i),g}.data(:,pairs(j,2))';
                subplot(length(pairs),length(posIdx),plotNoPos);
                plot(stim,'*-k')
                hold on
                if j==1
                    ylim([ 0 3])
                else
                    ylim([0 1]);
                end
                xlimit = [0.5 2.5];
                if i==1;
                    ylabel(strcat(variableIDs{j},' ' , 'pos'));
                end
                text(xlimit(1),0.2,(posMouseIDs{i}));
                set(gca,'XTick',(1:2));
                set(gca,'XTickLabel',['Cntr';'Stim']);
                xlimit = [0.5 2.5];
                xlim(xlimit);
                plotNoPos=plotNoPos+1;
                avgDATA(i,:)= mean(stim');
                semDATA(i,:)= std(stim')/length(stim');
            end
        end
        [~,p] = ttest(avgDATA(:,1),avgDATA(:,2),'Alpha',0.05);
        
        
        figure(111*g)
        subplot(length(pairs),2,plotNO)
        %plot(avgDATA','k-');
        hold on
        errorbar(avgDATA',semDATA','k-','LineWidth',2);
        ylabel(strcat(variableIDs{j},' ' , 'pos'));
        set(gca,'XTick',(1:2));
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        
        if j==1
            ylim([ 0 3]);
        else
            ylim([0 1]);
        end
        str = ['P = ', num2str(p)];
        text(xlimit(1),0.2,str);
        plotNO=plotNO+1;
        clear avgDATA semDATA stim;
        figure(12*g)
        
        for i = 1:length(negIdx)
            clear stim
            if (~isempty(DATA{negIdx(i),g}))
                stim(1,:) = DATA{negIdx(i),g}.data(:,pairs(j,1))';
                stim(2,:) = DATA{negIdx(i),g}.data(:,pairs(j,2))';
                subplot(length(pairs),length(negIdx),plotNoNeg);
                plot(stim,'*-k')
                hold on
                if j==1
                    ylim([ 0 3])
                    text(xlimit(1),3.2,(negMouseIDs{i}));
                else
                    ylim([0 1]);
                end
                xlimit = [0.5 2.5];
                if i==1;
                    ylabel(strcat(variableIDs{j},' ' , 'neg'));
                end
                
                set(gca,'XTick',(1:2));
                set(gca,'XTickLabel',['Cntr';'Stim']);
                xlim(xlimit)
                plotNoNeg=plotNoNeg+1;
                avgDATA(i,:)= mean(stim');
                semDATA(i,:)= std(stim')/length(stim');
            end
        end
        [h,p] = ttest(avgDATA(:,1),avgDATA(:,2),'Alpha',0.05);
        
        figure(111*g)
        subplot(length(pairs),2,plotNO)
        %plot(avgDATA','k-');
        hold on
        errorbar(avgDATA',semDATA','k-','LineWidth',2);
        ylabel(strcat(variableIDs{j},' ' , 'neg'));
        set(gca,'XTick',(1:2));
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        str = ['P = ', num2str(p)];
        text(xlimit(1),0.2,str);
        if j==1
            ylim([ 0 3])
        else
            ylim([0 1]);
        end
        plotNO=plotNO+1;
    end
end

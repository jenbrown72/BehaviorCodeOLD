posMouseIDs = {'JB034_L';'JB034_R';'JB035_L';'JB035_NT';'JB035_RL'};
mouseIDs = {'JB034_L';'JB034_R';'JB035_L';'JB035_NT';'JB035_RL'};
% negMouseIDs = {'JB025_R';'JB025_RL';'JB026_CE';'JB026_L';'JB026_R';'JB025_L';'JB024_R'};
%posMouseIDs = {'JB025_CE'};
% negMouseIDs = {'JB025_R'};
variableIDs = {'dprime';'lick bias';'performance'};
tempstimulus = {'S8';'S6';'S2';'S1'};

for g=1; %ipsi or contra
    plotNoPos=1;
    pairs = [1,2;3,4;5,6]; %row pairs.. control stim control stim control stim
    pairsStimulus = [1,2;3,4;5,6;7,8]; %row pairs.. control stim control stim control stim
    plotNO=1;
    
    figure(11*g)
    set(gcf,'name','PvAi32 +ve Contralateral Stimulation','numbertitle','off')
    
    figure(111*g)
        set(gcf,'name','Avg Contralateral Stimulation','numbertitle','off')
    
    n=1;
    clear posIdx;
    for i=1:length(posMouseIDs);
        for ii = 1:length(mouseIDs)
            if(strcmp(posMouseIDs{i},mouseIDs{ii,g}))
                posIdx(n,1) = ii;
                n=n+1;
            end
        end
    end
    
    for j = 1:length(pairs); %for each parameter (dprime, lick bias, performance)
        clear avgDATA semDATA stim;
        figure(11*g)
        hold on
        for i = 1:length(posIdx); %for each mouse
            clear stim
            if (~isempty(DATA{posIdx(i),g}))
                
                for kk = 1:length(pairsStimulus);
                    data = DATA{posIdx(i),g}.data(pairsStimulus(kk,:),pairs(j,1));
                    [meanData,~,~,semData] = JB_calBasicStats(data);
                    
                    stim(1,kk) = meanData;
                    stimSEM(1,kk) = semData;
                    
                    data = DATA{posIdx(i),g}.data(pairsStimulus(kk,:),pairs(j,2));
                    [meanData,stdData,countData,semData] = JB_calBasicStats(data);
                    
                    stim(2,kk) = meanData;
                    stimSEM(2,kk) = semData;

                end
                subplot(length(pairs),length(posIdx),plotNoPos);
                errorbar(stim(1,:),stimSEM(1,:),'o-k','MarkerSize',3,'LineWidth',2);
                hold on
                errorbar(stim(2,:),stimSEM(2,:),'o-b','MarkerSize',3,'LineWidth',2);
                if j==1
                    ylim([-1 4])
                else
                    ylim([0 1]);
                end
                
                if j>1
        line([0 5],[0.5 0.5],'LineStyle','--','Color','k')
                end
           
                if i==1;
                    ylabel(strcat(variableIDs{j},' ' , 'pos'));
                end
          
                text(xlimit(1),0.2,(posMouseIDs{i}));
                set(gca,'XTick',(1:4));
                set(gca,'XTickLabel',(tempstimulus)');
                     xlimit = [0.5 4.5];
                set(gca,'Xdir','reverse')
                plotNoPos=plotNoPos+1;
                %                 avgDATA(i,:)= mean(stim');
                %                 semDATA(i,:)= std(stim')/length(stim');
            end
        end

    end
end

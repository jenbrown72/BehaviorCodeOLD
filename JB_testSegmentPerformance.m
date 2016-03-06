
% for j=1:length(data.trial);
%     data.CumulativePerformance(j,1) = data.trial{j,1}.Performance.StimType;
%     data.CumulativePerformance(j,2) = data.trial{j,1}.Performance.OutcomeLick;
%     data.CumulativePerformance(j,3) = data.trial{j,1}.Performance.LickFrequency;
%     data.CumulativePerformance(j,4) = data.trial{j,1}.Performance.FirstLickLatency;
%     data.CumulativePerformance(j,5) = data.trial{j,1}.Performance.OutcomeCode;
%     data.CumulativePerformance(j,6) = data.trial{j,1}.Performance.StimTrial;
% end
fontSize = 16; %font Size for plotting

prompt = 'If this was a segmented session - input time points';
sessionDivisions = input(prompt);
%sessionDivisions = [1 136;137 237;238 337;338 459];
segmentPerformance = nan(length(sessionDivisions),2);

for j = 1:length(sessionDivisions)
    A = [data.CumulativePerformance(:,7)];
    
    if sessionDivisions(j,1)==1;
        idxStart = 1;
    else
        
        idxStart = find(A-sessionDivisions(j,1)< 0, 1,'last');
    end
    
    idxEnd = find(A-sessionDivisions(j,2)< 0, 1,'last');
    
    tempDATA = data.CumulativePerformance(idxStart:idxEnd,:);
    HITangle = sum(tempDATA(:,5)==1);
    FAangle = sum(tempDATA(:,5)==2);
    MISSangle = sum(tempDATA(:,5)==3);
    CRangle = sum(tempDATA(:,5)==4);
    segmentPerformance(j,1) =  (HITangle+CRangle)/(length(tempDATA));
end

figure;clf;
plot(segmentPerformance,'ok-')
hold on
line([1 max(xlim)],[0.5 0.5],'Color','r','LineStyle','--','LineWidth',2)
SessionTypes = {'Cntr' ; 'Stim out'; 'Cntr'; 'Stim out'};
NumTicks = length(SessionTypes);
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks));

ylim([0 1]);
title('Performance over blocks', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Block number', 'FontSize', fontSize,'fontWeight','bold');
ylabel('Percent Correct ((hit+CR)/totalTrials', 'FontSize', fontSize,'fontWeight','bold');
set(gca,'XTickLabel',SessionTypes)
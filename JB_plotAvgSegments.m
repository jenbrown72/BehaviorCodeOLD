function [ avg ] = JB_plotAvgSegments(performance,dprime)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


figure(1);clf;
currPlot=1;


varToPlot = {'performance';'dprime'};
threshold = [0.6,1];
axisLimY = [0 1;-1 3];
axisLimX = [0 4;0 4];
for k = 1:length(varToPlot);
dataTemp = eval(varToPlot{k});
figure(1)
subplot(2,2,currPlot)
tally1 = 1;
tally2 = 1;
for i=1:length(dataTemp);
    if (dataTemp(1,i)>threshold(k))
        plot(dataTemp(:,i),'ro-','MarkerSize',4)
        avg.abovethres.(varToPlot{k})(:,tally1) = (dataTemp(:,i));
        tally1=tally1+1;
        hold on
    else
        plot(dataTemp(:,i),'ko-','MarkerSize',4)
        avg.belowthres.(varToPlot{k})(:,tally2) = (dataTemp(:,i));
        tally2=tally2+1;
        hold on
    end
    ylim(axisLimY(k,:))
    xlim(axisLimX(k,:));
    ylabel((varToPlot{k}))
end


avg.(varToPlot{k})(:,1) = (mean(avg.abovethres.(varToPlot{k})'));
avg.(varToPlot{k})(:,2) = (mean(avg.belowthres.(varToPlot{k})'));
figure(1)
subplot(2,2,currPlot+1)

b = bar(avg.(varToPlot{k}));
b(1).FaceColor = 'red';
b(2).FaceColor = 'w';
ylim(axisLimY(k,:));
xlim(axisLimX(k,:));
 ylabel((varToPlot{k}));
currPlot = currPlot+2;
end






for i=1:length(data.trial);
    figure(1);clf
trialNo = i;

plot(data.trial{1,trialNo}.raw(:,2),'Color',rgb('Red'),'linewidth',3) %target position
hold on
plot(data.trial{1,trialNo}.raw(:,4),'Color',rgb('Lime'))%Lick L 
plot(data.trial{1,trialNo}.raw(:,3),'Color',rgb('DarkGreen'),'linewidth',3)%Lick L Reward
plot(data.trial{1,trialNo}.raw(:,16),'Color',rgb('Magenta')) %Lick R
plot(data.trial{1,trialNo}.raw(:,17),'Color',rgb('Indigo'),'linewidth',3)%Lick R reward

xlim([find(data.trial{1,trialNo}.raw(:,2)>0,1)-500 find(data.trial{1,trialNo}.raw(:,2)>0,1)+1000]) 
ylim([-0.1 1.1])
    xlimit = xlim;
  str = strcat('session Angle ',num2str(data.trial{1,trialNo}.Performance.StimType));
text(xlimit(1),0,str)
if (data.trial{1,trialNo}.Performance.StimType<270)
    text(xlimit(1),-0.1,'LEFT')
else
        text(xlimit(1),-0.1,'RIGHT')
end
legend('Target','L Lick','L Reward','R Lick','R Reward')
waitforbuttonpress;
end
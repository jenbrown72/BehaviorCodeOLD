
positionGraph1 = [  3   778   861   217];
varLabels = {'dprime';'lick bias';'performance'};
mouseIDs = {'JB023_NT';'JB024_R';'JB024_NT';'JB025_L';'JB025_R';'JB025_RL';'JB025_CE';'JB026_L';'JB026_R';'JB026_CE'}

for k = 1:length(mouseIDs)
figure(k)
 set(gcf,'Position',positionGraph1);
  set(gcf,'name',mouseIDs{k},'numbertitle','off')
  countPlot=1;
  countPlot2=1;
for i=1:3;

subplot(2,3,i)
tempDATA = DATA{k,1}.data(:,countPlot:countPlot+1);
  meandPrime = mean(tempDATA);
  contraSide(k,:) = meandPrime;
            SEMdPrime = std(tempDATA)/sqrt(length(tempDATA));
            [h,p] = ttest(tempDATA(:,1),tempDATA(:,2));
            if h==1
                sigPairs = [1,2];
                sigPairsPval = p;
            end
            
plot(tempDATA','k-')
           hold on
            errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
ylabel(varLabels{i});
 set(gca,'XTickLabel',['Cntr';'Stim']);
            set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',['Cntr';'Stim']);
            xlimit = [0.5 2.5];
            set(gca, 'Xlim',[xlimit(1) xlimit(2)])
            if i>1
                   set(gca, 'Ylim',[0 1])
            end
                str = ['P = ', num2str(p)];
            text(xlimit(1),0.2,str)
                
countPlot=countPlot+2;
hold on

subplot(2,3,i+3)
tempDATA = DATA{k,2}.data(:,countPlot2:countPlot2+1);
  meandPrime = mean(tempDATA);
  ipsiSide(k,:) = meandPrime;
            SEMdPrime = std(tempDATA)/sqrt(length(tempDATA));
            [h,p] = ttest(tempDATA(:,1),tempDATA(:,2));
            if h==1
                sigPairs = [1,2];
                sigPairsPval = p;
            end
            
plot(tempDATA','k-')
           hold on
            errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)

ylabel(varLabels{i});

 set(gca,'XTickLabel',['Cntr';'Stim']);
            set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',['Cntr';'Stim']);
            xlimit = [0.5 2.5];
            set(gca, 'Xlim',[xlimit(1) xlimit(2)])
                str = ['P = ', num2str(p)];
            text(xlimit(1),0.2,str)
countPlot2=countPlot2+2;
hold on
end
end


figure(11)
subplot(1,2,1)
plot(contraSide','ok-')
ylabel('contralateral stimulation')
      set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',['Cntr';'Stim']);
subplot(1,2,2)
plot(ipsiSide','ok-')
ylabel('ipsilateral stimulation')
           set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',['Cntr';'Stim']);


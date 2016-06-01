function [running_speed,time_bins] = JB_computeSpeed(AvgVelocity,sampleRateSpacing,TimeScaled,run,binSize)
 
% JB_computeSpeed uses data from wheel encoder to extract run speed
 % 
 
 %   [] = JB_computeSpeed(AvgVelocity,sampleRateSpacing,TimeScaled)
 %   AvgVelocity 
 %   sampleRateSpacing sample rate of the session
 %   TimeScaled time data
 

 % Examples:
 %   [running_speed,time_bins] = JB_computeSpeed(AvgVelocity,sampleRateSpacing,TimeScaled);
 %   
 
% default rotary encoder has 360 pulses per revolution
% default circumference is 2*pi*7.6 cm = 47.75cm/revolution or 18.84 inches/revolution
% 
% MotorA = velocityTemp;
% %speed is determined by MotorA right now
% runTicks = diff(MotorA);

if (nargin>4)
    binWidth = binSize;
else
   binWidth = 20; % in ms binSize
end

runTicks = AvgVelocity*1000;
binInPnts = sampleRateSpacing*binWidth;

binnedTicks=[];
nBins = length(runTicks) /binInPnts;
for i=1:nBins
    binnedTicks(i)  = sum(runTicks((i-1)*binInPnts+1:i*binInPnts));
    time_bins(i)  = sum(TimeScaled((i-1)*binInPnts+1:i*binInPnts));
end

if (run==1) 
circumference = 47.75; %cm at edge of wheel
angularDistance = (binnedTicks)/(360*2);
LinearDistance = angularDistance * circumference;
running_speed = LinearDistance / 1000; %run speed cm/s
else    
running_speed = binnedTicks / (binWidth/1); %run speed cm/s
end   
time_bins = time_bins/(binWidth/1);

end

